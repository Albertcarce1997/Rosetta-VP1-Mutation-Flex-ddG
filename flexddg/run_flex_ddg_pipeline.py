#!/usr/bin/env python3
"""Run Rosetta flex ddG analyses for VP1 haplotypes.

This script is repository-relative: after cloning the repository, run the
preparation scripts first, then run this file from either the repository root or
the flexddg directory.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import logging
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
from dataclasses import dataclass
from io import StringIO
from pathlib import Path

import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
DEFAULT_PREP_DIR = REPO_ROOT / "prep"
DEFAULT_XML_SCRIPT = SCRIPT_DIR / "flex_ddg_pipeline" / "ddG-backrub.xml"
DEFAULT_REFERENCE_FASTA = SCRIPT_DIR / "Sabin_2_VP1.fasta"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "flex_ddg_results"
FASTA_EXTENSIONS = {".fa", ".faa", ".fas", ".fasta"}

CRYSTAL_MIN_VP1_POSITION = 27
CRYSTAL_MISSING_PDB_RESIDUES = set(range(96, 102))
MUTATION_PATTERN = re.compile(r"^([A-Z])(\d+)([A-Z])$")


@dataclass(frozen=True)
class TargetConfig:
    name: str
    pdb_path: Path
    chains_to_move: str
    vp1_chain: str = "A"
    extra_res_fa: Path | None = None

    @property
    def interface_partners(self) -> str:
        partner_chains = self.chains_to_move.replace(",", "")
        return f"{self.vp1_chain}_{partner_chains}"


@dataclass
class RunConfig:
    haplotype_fasta: Path
    reference_fasta: Path
    xml_script: Path
    output_dir: Path
    rosetta_image: str
    rosetta_bin: str
    num_backrub_trials: int
    num_replicates: int
    max_minimization_iter: int
    abs_score_convergence_thresh: float
    backrub_trajectory_stride: int
    max_workers: int
    min_vp1_position: int
    missing_pdb_residues: set[int]
    redo_excel: bool
    force: bool
    compress_pdbs: bool
    dry_run: bool
    constant_seed: bool
    seed_offset: int


def resolve_cli_path(value: str | Path, base_dir: Path = REPO_ROOT) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = base_dir / path
    return path.resolve()


def is_reference_like_fasta(path: Path, reference_fasta: Path) -> bool:
    if path.resolve() == reference_fasta.resolve():
        return True
    stem = path.stem.lower()
    return stem in {"ref", "refs", "reference", "references"} or "reference" in stem or "sabin_2" in stem or "sabin2" in stem


def discover_haplotype_fasta(reference_fasta: Path, search_dir: Path = SCRIPT_DIR) -> Path:
    candidates = sorted(
        path.resolve()
        for path in search_dir.iterdir()
        if path.is_file() and path.suffix.lower() in FASTA_EXTENSIONS and not is_reference_like_fasta(path, reference_fasta)
    )

    if len(candidates) == 1:
        return candidates[0]

    if not candidates:
        print(
            f"No haplotype FASTA file was found in {search_dir}. "
            "Add one FASTA file next to the flex ddG script or pass it with --input.",
            file=sys.stderr,
        )
        sys.exit(1)

    print(
        "Multiple non-reference FASTA files were found. Choose one with --input:",
        file=sys.stderr,
    )
    for candidate in candidates:
        print(f"  {candidate}", file=sys.stderr)
    sys.exit(1)


def resolve_haplotype_fasta(input_fasta: str | None, reference_fasta: Path) -> Path:
    if input_fasta:
        return resolve_cli_path(input_fasta)
    return discover_haplotype_fasta(reference_fasta)


def repo_relative(path: Path) -> str:
    resolved_repo = REPO_ROOT.resolve()
    resolved_path = path.resolve()
    try:
        return resolved_path.relative_to(resolved_repo).as_posix()
    except ValueError as error:
        raise ValueError(f"Path must be inside the repository root: {resolved_path}") from error


def format_time(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.0f}s"
    if seconds < 3600:
        minutes = seconds // 60
        remaining_seconds = seconds % 60
        return f"{minutes:.0f}m {remaining_seconds:.0f}s"
    hours = seconds // 3600
    minutes = (seconds % 3600) // 60
    return f"{hours:.0f}h {minutes:.0f}m"


def configure_logging(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / "flex_ddg_pipeline.log"

    logger = logging.getLogger()
    logger.handlers.clear()
    logger.setLevel(logging.INFO)

    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)


def build_targets(prep_dir: Path, target_selection: str) -> list[TargetConfig]:
    pocket_dir = prep_dir / "pocket_factor_PV2"
    all_targets = {
        "receptor": TargetConfig(
            name="Receptor",
            pdb_path=prep_dir / "VP1-receptor_final.pdb",
            chains_to_move="R",
        ),
        "9H2": TargetConfig(
            name="Antibody_9H2",
            pdb_path=prep_dir / "VP1-9H2_final.pdb",
            chains_to_move="H,L",
        ),
        "10D2": TargetConfig(
            name="Antibody_10D2",
            pdb_path=prep_dir / "VP1-10D2_final.pdb",
            chains_to_move="H,L",
        ),
        "pocket_factor": TargetConfig(
            name="Pocket_factor",
            pdb_path=pocket_dir / "8E8Y_VP1_PLM_final_ligX.pdb",
            chains_to_move="X",
            extra_res_fa=pocket_dir / "PLM.params",
        ),
    }

    selection_map = {
        "all": ["receptor", "9H2", "10D2", "pocket_factor"],
        "protein_complexes": ["receptor", "9H2", "10D2"],
        "receptor": ["receptor"],
        "antibodies": ["9H2", "10D2"],
        "9H2": ["9H2"],
        "10D2": ["10D2"],
        "pocket_factor": ["pocket_factor"],
        "plm": ["pocket_factor"],
    }
    return [all_targets[target_key] for target_key in selection_map[target_selection]]


def validate_inputs(targets: list[TargetConfig], config: RunConfig) -> None:
    required_paths = [config.haplotype_fasta, config.reference_fasta, config.xml_script]
    for target in targets:
        required_paths.append(target.pdb_path)
        if target.extra_res_fa is not None:
            required_paths.append(target.extra_res_fa)

    missing_paths = [path for path in required_paths if not path.exists()]
    if missing_paths:
        logging.error("Required input files are missing:")
        for path in missing_paths:
            logging.error("  %s", path)
        logging.error("Run the preparation scripts first, or pass explicit paths with the CLI options.")
        sys.exit(1)


def check_external_tools(config: RunConfig, needs_mafft: bool) -> None:
    if needs_mafft and shutil.which("mafft") is None:
        logging.error("MAFFT was not found in PATH. Install the conda environment or add MAFFT to PATH.")
        sys.exit(1)

    if not config.dry_run and shutil.which("docker") is None:
        logging.error("Docker was not found in PATH. Docker is required to run the Rosetta container.")
        sys.exit(1)


def docker_user_args() -> list[str]:
    if os.name != "posix":
        return []
    try:
        return ["--user", f"{os.getuid()}:{os.getgid()}"]
    except AttributeError:
        return []


def build_rosetta_command(
    target: TargetConfig,
    resfile_path: Path,
    output_dir: Path,
    config: RunConfig,
    *,
    nstruct: int,
    suffix: str = "",
    seed: int | None = None,
) -> list[str]:
    cmd = [
        "docker",
        "run",
        "--rm",
        "-v",
        f"{REPO_ROOT.resolve()}:/data",
        "-w",
        "/data",
    ]
    cmd.extend(docker_user_args())
    cmd.extend(
        [
            config.rosetta_image,
            config.rosetta_bin,
            "-parser:protocol",
            repo_relative(config.xml_script),
            "-s",
            repo_relative(target.pdb_path),
        ]
    )

    if target.extra_res_fa is not None:
        cmd.extend(["-extra_res_fa", repo_relative(target.extra_res_fa)])

    scorefile_name = f"score{suffix}.sc" if suffix else "score.sc"
    cmd.extend(
        [
            "-parser:script_vars",
            f"chainstomove={target.chains_to_move}",
            f"interface_partners={target.interface_partners}",
            f"mutate_resfile_relpath={repo_relative(resfile_path)}",
            f"number_backrub_trials={config.num_backrub_trials}",
            f"max_minimization_iter={config.max_minimization_iter}",
            f"abs_score_convergence_thresh={config.abs_score_convergence_thresh}",
            f"backrub_trajectory_stride={config.backrub_trajectory_stride}",
            "-in:file:fullatom",
            "-ignore_unrecognized_res",
            "-ignore_zero_occupancy",
            "false",
            "-ex1",
            "-ex2",
            "-nstruct",
            str(nstruct),
            "-out:path:all",
            repo_relative(output_dir),
            "-out:file:scorefile",
            scorefile_name,
        ]
    )

    if suffix:
        cmd.extend(["-out:suffix", suffix])
    if config.constant_seed:
        cmd.extend(["-constant_seed", "-jran", str(seed if seed is not None else config.seed_offset)])

    return cmd


def run_command(cmd: list[str], description: str, config: RunConfig) -> bool:
    if config.dry_run:
        logging.info("DRY RUN %s", description)
        logging.info("%s", " ".join(cmd))
        return True

    logging.info("Running: %s", description)
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
    if result.returncode == 0:
        return True

    logging.error("Command failed for %s with exit code %s", description, result.returncode)
    if result.stdout:
        logging.error("STDOUT tail:\n%s", result.stdout[-4000:])
    if result.stderr:
        logging.error("STDERR tail:\n%s", result.stderr[-4000:])
    return False


def align_sequences(ref_record: SeqRecord, ht_record: SeqRecord, output_fasta: Path, temp_dir: Path) -> bool:
    temp_dir.mkdir(parents=True, exist_ok=True)
    temp_input_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False, dir=temp_dir, encoding="utf-8") as temp_input:
            SeqIO.write([ref_record, ht_record], temp_input, "fasta")
            temp_input_path = Path(temp_input.name)

        cmd = ["mafft", "--localpair", "--maxiterate", "1000", "--quiet", str(temp_input_path)]
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
        if result.returncode != 0:
            logging.error("MAFFT failed for %s: %s", ht_record.id, result.stderr)
            return False
        output_fasta.write_text(result.stdout, encoding="utf-8")
        return True
    except Exception as error:
        logging.error("Alignment failed for %s: %s", ht_record.id, error)
        return False
    finally:
        if temp_input_path is not None and temp_input_path.exists():
            temp_input_path.unlink()


def parse_mutations(alignment_file: Path, ref_id: str, ht_id: str) -> tuple[list[tuple[str, int, str]], str]:
    alignment = AlignIO.read(alignment_file, "fasta")
    ref_sequence = None
    ht_sequence = None

    for record in alignment:
        if record.id == ref_id:
            ref_sequence = str(record.seq)
        elif record.id == ht_id:
            ht_sequence = str(record.seq)

    if not ref_sequence or not ht_sequence:
        return [], "Error: sequences not found in alignment"

    start_gap_length = 0
    in_start_gap = True
    for amino_acid in ht_sequence:
        if amino_acid == "-" and in_start_gap:
            start_gap_length += 1
        else:
            in_start_gap = False

    end_gap_length = 0
    for amino_acid in reversed(ht_sequence):
        if amino_acid == "-":
            end_gap_length += 1
        else:
            break

    if start_gap_length > 20:
        return [], "Skipped: start gap > 20"
    if end_gap_length > 7:
        return [], "Skipped: end gap > 7 (affecting antigenic site 3)"

    mutations: list[tuple[str, int, str]] = []
    reference_position = 0
    alignment_length = len(ref_sequence)

    for alignment_index, (ref_amino_acid, ht_amino_acid) in enumerate(zip(ref_sequence, ht_sequence)):
        if ref_amino_acid != "-":
            reference_position += 1

        if alignment_index < start_gap_length:
            continue
        if alignment_index >= alignment_length - end_gap_length:
            continue

        if ref_amino_acid == "-":
            return [], f"Skipped: contains insertion after RefPos {reference_position}. Docking required."
        if ht_amino_acid == "-":
            return [], f"Skipped: contains deletion at RefPos {reference_position}. Docking required."
        if ref_amino_acid != ht_amino_acid:
            mutations.append((ref_amino_acid, reference_position, ht_amino_acid))

    return mutations, "OK"


def create_resfile(mutations: list[tuple[str, int, str]], chain: str, filename: Path, min_seq_pos: int) -> None:
    filename.parent.mkdir(parents=True, exist_ok=True)
    with filename.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("NATAA\n")
        handle.write("start\n")
        for wild_type, position, mutant in mutations:
            if position >= min_seq_pos:
                handle.write(f"{position} {chain} PIKAA {mutant}\n")


def parse_flex_ddg_output(output_dir: Path) -> dict[str, float | int | str] | None:
    score_files = sorted(output_dir.glob("*.sc"))
    ddg_files = sorted(output_dir.glob("*ddg_predictions.out"))
    result_files = ddg_files + score_files
    if not result_files:
        return None

    dataframes = []
    for result_file in result_files:
        try:
            with result_file.open("r", encoding="utf-8", errors="ignore") as handle:
                score_lines = [line for line in handle if line.startswith("SCORE:")]
            if not score_lines:
                logging.warning("No SCORE lines found in %s", result_file)
                continue
            dataframe = pd.read_csv(StringIO("".join(score_lines)), sep=r"\s+")
            if "SCORE:" in dataframe.columns:
                dataframe = dataframe.drop(columns=["SCORE:"])
            dataframes.append(dataframe)
        except Exception as error:
            logging.warning("Failed to parse score file %s: %s", result_file, error)

    if not dataframes:
        return None

    combined = pd.concat(dataframes, ignore_index=True)
    preferred_columns = ["ddG", "ddg", "dG_separated", "dG_bind", "dG_binding", "I_sc", "total_score", "score"]
    metric_column = next((column for column in preferred_columns if column in combined.columns), None)
    if metric_column is None:
        logging.error("No recognized metric column found in %s. Columns: %s", output_dir, list(combined.columns))
        return None

    metric_values = pd.to_numeric(combined[metric_column], errors="coerce").dropna()
    if metric_values.empty:
        logging.error("Metric column %s in %s did not contain numeric values", metric_column, output_dir)
        return None

    return {
        "mean": float(metric_values.mean()),
        "std": float(metric_values.std()),
        "sem": float(metric_values.sem()),
        "n": int(metric_values.shape[0]),
        "metric": metric_column,
    }


def compress_and_cleanup_pdbs(output_dir: Path, enabled: bool) -> None:
    if not enabled:
        return
    pdb_files = sorted(output_dir.glob("**/*.pdb"))
    if not pdb_files:
        return

    archive_path = output_dir / "structures.tar.gz"
    logging.info("Compressing %s PDB files into %s", len(pdb_files), archive_path)
    try:
        with tarfile.open(archive_path, "w:gz") as archive:
            for pdb_file in pdb_files:
                archive.add(pdb_file, arcname=pdb_file.relative_to(output_dir))
        for pdb_file in pdb_files:
            pdb_file.unlink()
    except Exception as error:
        logging.error("Failed to compress PDBs in %s: %s", output_dir, error)


def extract_haplotype_id(record: SeqRecord) -> str:
    match = re.search(r"HT-\d+", record.description)
    if match:
        return match.group(0)
    cleaned_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", record.id)
    return cleaned_id or "haplotype"


def mutation_to_string(mutation: tuple[str, int, str]) -> str:
    wild_type, position, mutant = mutation
    return f"{wild_type}{position}{mutant}"


def mutation_list_to_string(mutations: list[tuple[str, int, str]] | tuple[tuple[str, int, str], ...]) -> str:
    return ", ".join(mutation_to_string(mutation) for mutation in mutations)


def filter_mutations(mutations: list[tuple[str, int, str]], config: RunConfig) -> tuple[list[tuple[str, int, str]], list[tuple[str, int, str]]]:
    filtered = []
    skipped = []
    for mutation in mutations:
        _wild_type, position, _mutant = mutation
        if position < config.min_vp1_position or position in config.missing_pdb_residues:
            skipped.append(mutation)
        else:
            filtered.append(mutation)
    return filtered, skipped


def analyze_haplotype(ht_record: SeqRecord, ref_record: SeqRecord, config: RunConfig) -> tuple[str, list[tuple[str, int, str]], str]:
    ht_id = extract_haplotype_id(ht_record)
    temp_align_dir = config.output_dir / "temp_alignments"
    alignment_file = temp_align_dir / f"{ht_id}_alignment.fasta"

    if not align_sequences(ref_record, ht_record, alignment_file, temp_align_dir):
        return ht_id, [], "Alignment failed"

    mutations, status = parse_mutations(alignment_file, ref_record.id, ht_record.id)
    try:
        alignment_file.unlink()
    except OSError:
        pass

    if status != "OK":
        return ht_id, [], status

    valid_mutations, skipped_mutations = filter_mutations(mutations, config)
    if skipped_mutations:
        logging.warning("%s: skipped mutations outside modeled residues: %s", ht_id, mutation_list_to_string(skipped_mutations))
    if not valid_mutations:
        return ht_id, [], "No qualifying mutations"
    return ht_id, valid_mutations, "OK"


def has_existing_scores(output_dir: Path) -> bool:
    return any(output_dir.glob("*.sc")) or any(output_dir.glob("*ddg_predictions.out"))


def run_execution_tasks(tasks: list[tuple[list[str], str]], config: RunConfig, pool_executor: concurrent.futures.Executor | None = None) -> None:
    if not tasks:
        return

    if pool_executor is not None:
        futures = [pool_executor.submit(run_command, cmd, description, config) for cmd, description in tasks]
        for future in concurrent.futures.as_completed(futures):
            future.result()
        return

    max_workers = min(config.max_workers, len(tasks))
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(run_command, cmd, description, config) for cmd, description in tasks]
        for future in concurrent.futures.as_completed(futures):
            future.result()


def run_rosetta_on_mutations(
    mutations: list[tuple[str, int, str]],
    output_dir: Path,
    label: str,
    targets: list[TargetConfig],
    config: RunConfig,
    *,
    parallel_replicates: bool = False,
    pool_executor: concurrent.futures.Executor | None = None,
) -> list[dict[str, object]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    tasks: list[tuple[list[str], str]] = []

    for target in targets:
        target_dir = output_dir / target.name
        target_dir.mkdir(parents=True, exist_ok=True)

        if has_existing_scores(target_dir) and not config.force:
            logging.info("%s %s: score files found; skipping Rosetta run", label, target.name)
            continue

        resfile_path = target_dir / "mutations.resfile"
        create_resfile(mutations, target.vp1_chain, resfile_path, config.min_vp1_position)

        if parallel_replicates:
            for replicate_index in range(config.num_replicates):
                replicate_number = replicate_index + 1
                suffix = f"_rep{replicate_number:03d}"
                cmd = build_rosetta_command(
                    target,
                    resfile_path,
                    target_dir,
                    config,
                    nstruct=1,
                    suffix=suffix,
                    seed=config.seed_offset + replicate_index,
                )
                tasks.append((cmd, f"{label} {target.name} replicate {replicate_number}"))
        else:
            cmd = build_rosetta_command(
                target,
                resfile_path,
                target_dir,
                config,
                nstruct=config.num_replicates,
                seed=config.seed_offset,
            )
            tasks.append((cmd, f"{label} {target.name}"))

    run_execution_tasks(tasks, config, pool_executor=pool_executor)

    results: list[dict[str, object]] = []
    for target in targets:
        target_dir = output_dir / target.name
        parsed = parse_flex_ddg_output(target_dir)
        if parsed:
            results.append(
                {
                    "Target": target.name,
                    "Mutation": label,
                    "Mutations_List": mutation_list_to_string(mutations) if mutations else "None",
                    **parsed,
                }
            )

    compress_and_cleanup_pdbs(output_dir, enabled=config.compress_pdbs)
    return results


def save_results_excel(results: list[dict[str, object]], excel_path: Path) -> None:
    if not results:
        return
    excel_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(results).to_excel(excel_path, index=False)


def process_item(
    item_type: str,
    item_id: str,
    mutations: list[tuple[str, int, str]],
    out_parent: Path,
    targets: list[TargetConfig],
    config: RunConfig,
    *,
    parallel_replicates: bool = False,
    pool_executor: concurrent.futures.Executor | None = None,
) -> tuple[str, list[dict[str, object]]]:
    out_dir = out_parent / item_id
    excel_path = out_dir / f"{item_id}_results.xlsx"

    if excel_path.exists() and not config.redo_excel and not config.force:
        logging.info("%s %s: existing Excel found; skipping", item_type, item_id)
        return item_id, []

    if out_dir.exists() and config.force:
        logging.info("%s %s: force enabled; removing %s", item_type, item_id, out_dir)
        shutil.rmtree(out_dir)
    elif out_dir.exists() and not excel_path.exists() and not config.redo_excel:
        logging.info("%s %s: incomplete output; cleaning %s", item_type, item_id, out_dir)
        shutil.rmtree(out_dir)

    out_dir.mkdir(parents=True, exist_ok=True)
    results = run_rosetta_on_mutations(
        mutations,
        out_dir,
        item_id,
        targets,
        config,
        parallel_replicates=parallel_replicates,
        pool_executor=pool_executor,
    )
    save_results_excel(results, excel_path)
    return item_id, results


def run_global_wt_baseline(targets: list[TargetConfig], config: RunConfig, parallel_replicates: bool) -> pd.DataFrame | None:
    wt_dir = config.output_dir / "WT_Baseline"
    wt_excel_path = wt_dir / "WT_flex_ddg_results.xlsx"

    if wt_excel_path.exists() and not config.redo_excel and not config.force:
        logging.info("WT baseline: existing Excel found; skipping")
        try:
            return pd.read_excel(wt_excel_path)
        except Exception as error:
            logging.warning("WT baseline Excel could not be read: %s", error)

    if wt_dir.exists() and config.force:
        shutil.rmtree(wt_dir)

    results = run_rosetta_on_mutations([], wt_dir, "WT_Baseline", targets, config, parallel_replicates=parallel_replicates)
    if results:
        dataframe = pd.DataFrame(results)
        dataframe.insert(0, "Haplotype", "WT")
        dataframe.to_excel(wt_excel_path, index=False)
        return dataframe
    return None


def load_reference_record(path: Path) -> SeqRecord:
    try:
        return SeqIO.read(path, "fasta")
    except ValueError:
        return next(SeqIO.parse(path, "fasta"))


def load_haplotype_records(path: Path) -> list[SeqRecord]:
    return list(SeqIO.parse(path, "fasta"))


def load_mutation_map(map_file: Path) -> tuple[dict[str, tuple[tuple[str, int, str], ...]], dict[tuple[tuple[str, int, str], ...], str], set[tuple[str, int, str]]]:
    map_df = pd.read_excel(map_file)
    ht_mutation_map: dict[str, tuple[tuple[str, int, str], ...]] = {}
    unique_mutation_sets: dict[tuple[tuple[str, int, str], ...], str] = {}
    unique_single_mutations: set[tuple[str, int, str]] = set()

    for _row_index, row in map_df.iterrows():
        ht_id = str(row["Haplotype"])
        set_id = str(row["Mutation_Set_ID"])
        mutations_text = str(row.get("Mutations", ""))
        mutation_list = []
        if mutations_text and mutations_text.lower() != "nan":
            for mutation_text in [part.strip() for part in mutations_text.split(",") if part.strip()]:
                match = MUTATION_PATTERN.match(mutation_text)
                if not match:
                    logging.warning("Could not parse mutation in map file: %s", mutation_text)
                    continue
                mutation_tuple = (match.group(1), int(match.group(2)), match.group(3))
                mutation_list.append(mutation_tuple)
                unique_single_mutations.add(mutation_tuple)

        mutation_tuple = tuple(sorted(mutation_list, key=lambda mutation: mutation[1]))
        ht_mutation_map[ht_id] = mutation_tuple
        unique_mutation_sets[mutation_tuple] = set_id

    return ht_mutation_map, unique_mutation_sets, unique_single_mutations


def analyze_haplotypes(
    ht_records: list[SeqRecord],
    ref_record: SeqRecord,
    config: RunConfig,
) -> tuple[dict[str, tuple[tuple[str, int, str], ...]], dict[tuple[tuple[str, int, str], ...], str], set[tuple[str, int, str]]]:
    ht_mutation_map: dict[str, tuple[tuple[str, int, str], ...]] = {}
    unique_mutation_sets: dict[tuple[tuple[str, int, str], ...], str] = {}
    unique_single_mutations: set[tuple[str, int, str]] = set()
    skipped_docking: list[tuple[str, str]] = []

    set_id_counter = 1
    for ht_record in ht_records:
        ht_id, mutations, status = analyze_haplotype(ht_record, ref_record, config)
        if status != "OK":
            logging.warning("Skipping %s: %s", ht_id, status)
            if "Docking required" in status:
                skipped_docking.append((ht_id, status))
            continue

        mutation_tuple = tuple(sorted(mutations, key=lambda mutation: mutation[1]))
        ht_mutation_map[ht_id] = mutation_tuple
        if mutation_tuple not in unique_mutation_sets:
            unique_mutation_sets[mutation_tuple] = f"Set_{set_id_counter:04d}"
            set_id_counter += 1
        unique_single_mutations.update(mutation_tuple)

    if skipped_docking:
        skipped_file = config.output_dir / "skipped_haplotypes_docking_required.txt"
        with skipped_file.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("Haplotypes skipped because they require docking:\n")
            for ht_id, reason in skipped_docking:
                handle.write(f"{ht_id}: {reason}\n")

    mapping_rows = []
    for ht_id, mutation_tuple in ht_mutation_map.items():
        set_id = unique_mutation_sets[mutation_tuple]
        mapping_rows.append(
            {
                "Haplotype": ht_id,
                "Mutation_Set_ID": set_id,
                "Mutations": mutation_list_to_string(mutation_tuple),
                "Mutation_Count": len(mutation_tuple),
            }
        )
    pd.DataFrame(mapping_rows).to_excel(config.output_dir / "Haplotype_Mutation_Map.xlsx", index=False)

    return ht_mutation_map, unique_mutation_sets, unique_single_mutations


def aggregate_optimized_results(
    ht_mutation_map: dict[str, tuple[tuple[str, int, str], ...]],
    unique_mutation_sets: dict[tuple[tuple[str, int, str], ...], str],
    unique_single_mutations: set[tuple[str, int, str]],
    config: RunConfig,
) -> None:
    set_results_cache: dict[str, pd.DataFrame] = {}
    single_results_cache: dict[str, pd.DataFrame] = {}

    for set_id in unique_mutation_sets.values():
        result_path = config.output_dir / "all_mutations" / set_id / f"{set_id}_results.xlsx"
        if result_path.exists():
            set_results_cache[set_id] = pd.read_excel(result_path)

    for mutation in unique_single_mutations:
        mutation_id = mutation_to_string(mutation)
        result_path = config.output_dir / "single_mutations" / mutation_id / f"{mutation_id}_results.xlsx"
        if result_path.exists():
            single_results_cache[mutation_id] = pd.read_excel(result_path)

    final_frames: list[pd.DataFrame] = []
    for ht_id, mutation_tuple in ht_mutation_map.items():
        set_id = unique_mutation_sets[mutation_tuple]
        if set_id in set_results_cache:
            set_df = set_results_cache[set_id].copy()
            set_df.insert(0, "Haplotype", ht_id)
            set_df["Type"] = "All_Mutations"
            final_frames.append(set_df)

        for mutation in mutation_tuple:
            mutation_id = mutation_to_string(mutation)
            if mutation_id in single_results_cache:
                single_df = single_results_cache[mutation_id].copy()
                single_df.insert(0, "Haplotype", ht_id)
                single_df["Type"] = "Single_Mutation"
                final_frames.append(single_df)

    wt_excel = config.output_dir / "WT_Baseline" / "WT_flex_ddg_results.xlsx"
    if wt_excel.exists():
        wt_df = pd.read_excel(wt_excel)
        if "Haplotype" not in wt_df.columns:
            wt_df.insert(0, "Haplotype", "WT")
        wt_df["Type"] = "WT_Baseline"
        final_frames.append(wt_df)

    if not final_frames:
        logging.warning("No results found to aggregate")
        return

    final_df = pd.concat(final_frames, ignore_index=True)
    final_df.to_excel(config.output_dir / "Optimized_ALL_HT_results.xlsx", index=False)
    logging.info("Saved optimized aggregate results")


def monitor_futures(futures: list[concurrent.futures.Future], label: str) -> None:
    total_tasks = len(futures)
    completed_tasks = 0
    start_time = time.time()
    for future in concurrent.futures.as_completed(futures):
        completed_tasks += 1
        elapsed_time = time.time() - start_time
        average_time = elapsed_time / completed_tasks
        remaining_tasks = total_tasks - completed_tasks
        estimated_remaining = average_time * remaining_tasks
        logging.info(
            "%s progress: %s/%s (%.1f%%), elapsed %s, ETA %s",
            label,
            completed_tasks,
            total_tasks,
            (completed_tasks / total_tasks) * 100.0,
            format_time(elapsed_time),
            format_time(estimated_remaining),
        )
        future.result()


def run_optimized_pipeline(ht_records: list[SeqRecord], ref_record: SeqRecord, targets: list[TargetConfig], config: RunConfig, dedicated_cores: bool) -> None:
    map_file = config.output_dir / "Haplotype_Mutation_Map.xlsx"
    if map_file.exists() and not config.force:
        logging.info("Using existing mutation map: %s", map_file)
        ht_mutation_map, unique_mutation_sets, unique_single_mutations = load_mutation_map(map_file)
    else:
        logging.info("Analyzing haplotypes for mutations")
        ht_mutation_map, unique_mutation_sets, unique_single_mutations = analyze_haplotypes(ht_records, ref_record, config)

    logging.info(
        "Found %s unique mutation sets and %s unique single mutations",
        len(unique_mutation_sets),
        len(unique_single_mutations),
    )
    if not unique_mutation_sets and not unique_single_mutations:
        return

    all_mutations_dir = config.output_dir / "all_mutations"
    single_mutations_dir = config.output_dir / "single_mutations"

    if dedicated_cores:
        manager_workers = min(4, max(1, config.max_workers))
        logging.info("Running dedicated-core mode with %s Rosetta workers", config.max_workers)
        with concurrent.futures.ThreadPoolExecutor(max_workers=config.max_workers) as rosetta_pool:
            with concurrent.futures.ThreadPoolExecutor(max_workers=manager_workers) as manager_pool:
                futures = []
                for mutation_tuple, set_id in unique_mutation_sets.items():
                    futures.append(
                        manager_pool.submit(
                            process_item,
                            "SET",
                            set_id,
                            list(mutation_tuple),
                            all_mutations_dir,
                            targets,
                            config,
                            parallel_replicates=True,
                            pool_executor=rosetta_pool,
                        )
                    )
                for mutation in unique_single_mutations:
                    mutation_id = mutation_to_string(mutation)
                    futures.append(
                        manager_pool.submit(
                            process_item,
                            "SINGLE",
                            mutation_id,
                            [mutation],
                            single_mutations_dir,
                            targets,
                            config,
                            parallel_replicates=True,
                            pool_executor=rosetta_pool,
                        )
                    )
                monitor_futures(futures, "Item")
    else:
        logging.info("Running optimized parallel-item mode with %s workers", config.max_workers)
        with concurrent.futures.ThreadPoolExecutor(max_workers=config.max_workers) as executor:
            futures = []
            for mutation_tuple, set_id in unique_mutation_sets.items():
                futures.append(executor.submit(process_item, "SET", set_id, list(mutation_tuple), all_mutations_dir, targets, config))
            for mutation in unique_single_mutations:
                mutation_id = mutation_to_string(mutation)
                futures.append(executor.submit(process_item, "SINGLE", mutation_id, [mutation], single_mutations_dir, targets, config))
            monitor_futures(futures, "Item")

    aggregate_optimized_results(ht_mutation_map, unique_mutation_sets, unique_single_mutations, config)


def run_legacy_pipeline(ht_records: list[SeqRecord], ref_record: SeqRecord, targets: list[TargetConfig], config: RunConfig, serial: bool) -> None:
    def run_one_haplotype(ht_record: SeqRecord) -> tuple[str, str]:
        ht_id, mutations, status = analyze_haplotype(ht_record, ref_record, config)
        ht_dir = config.output_dir / ht_id
        ht_dir.mkdir(parents=True, exist_ok=True)
        if status != "OK":
            (ht_dir / "status.txt").write_text(status, encoding="utf-8")
            return ht_id, status

        (ht_dir / "mutations.txt").write_text("\n".join(mutation_to_string(mutation) for mutation in mutations), encoding="utf-8")
        results = []
        results.extend(run_rosetta_on_mutations(mutations, ht_dir / "All_Mutations", "ALL", targets, config))
        for mutation in mutations:
            mutation_id = mutation_to_string(mutation)
            results.extend(run_rosetta_on_mutations([mutation], ht_dir / "Single_Mutations" / mutation_id, mutation_id, targets, config))
        if results:
            dataframe = pd.DataFrame(results)
            dataframe.insert(0, "Haplotype", ht_id)
            dataframe.to_excel(ht_dir / f"{ht_id}_flex_ddg_results.xlsx", index=False)
        return ht_id, "OK"

    if serial:
        for ht_record in ht_records:
            run_one_haplotype(ht_record)
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=config.max_workers) as executor:
            futures = [executor.submit(run_one_haplotype, ht_record) for ht_record in ht_records]
            monitor_futures(futures, "Haplotype")

    result_frames = []
    wt_excel = config.output_dir / "WT_Baseline" / "WT_flex_ddg_results.xlsx"
    if wt_excel.exists():
        result_frames.append(pd.read_excel(wt_excel))
    for excel_path in config.output_dir.glob("HT-*/*_flex_ddg_results.xlsx"):
        result_frames.append(pd.read_excel(excel_path))
    if result_frames:
        pd.concat(result_frames, ignore_index=True).to_excel(config.output_dir / "ALL_HT_flex_ddg_results.xlsx", index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Rosetta flex ddG for VP1 haplotypes.")
    parser.add_argument(
        "--input",
        "--haplotypes",
        dest="input_fasta",
        default=None,
        help="FASTA file containing VP1 haplotype sequences. If omitted, auto-detects one non-reference FASTA in flexddg/.",
    )
    parser.add_argument("--reference", default=str(DEFAULT_REFERENCE_FASTA), help="Sabin 2 VP1 reference FASTA.")
    parser.add_argument("--prep-dir", default=str(DEFAULT_PREP_DIR), help="Directory containing prepared target PDB files.")
    parser.add_argument("--xml-script", default=str(DEFAULT_XML_SCRIPT), help="RosettaScripts XML file for flex ddG.")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR), help="Output directory for Rosetta results.")
    parser.add_argument(
        "--targets",
        choices=["all", "protein_complexes", "receptor", "antibodies", "9H2", "10D2", "pocket_factor", "plm"],
        default="all",
        help="Targets to process.",
    )
    parser.add_argument("--rosetta-image", default="rosettacommons/rosetta:ml-408", help="Docker image containing Rosetta.")
    parser.add_argument("--rosetta-bin", default="/usr/local/bin/rosetta_scripts", help="rosetta_scripts path inside the Docker image.")
    parser.add_argument("--num-backrub-trials", type=int, default=1000, help="Backrub trials per trajectory.")
    parser.add_argument("--num-replicates", type=int, default=35, help="Independent Rosetta trajectories per mutation set.")
    parser.add_argument("--max-minimization-iter", type=int, default=5000, help="Maximum minimization iterations.")
    parser.add_argument("--abs-score-convergence-thresh", type=float, default=200.0, help="Rosetta score convergence threshold.")
    parser.add_argument("--backrub-trajectory-stride", type=int, default=1000, help="Backrub trajectory stride.")
    parser.add_argument("--max-workers", type=int, default=min(os.cpu_count() or 1, 64), help="Maximum concurrent Rosetta jobs/items.")
    parser.add_argument("--legacy", action="store_true", help="Process each haplotype independently instead of deduplicating mutation sets.")
    parser.add_argument("--serial", action="store_true", help="In legacy mode, process haplotypes serially.")
    parser.add_argument("--dedicated-cores", action="store_true", help="Run each replicate as a separate Rosetta job in optimized mode.")
    parser.add_argument("--redo-excel", action="store_true", help="Regenerate Excel summaries from existing score files when possible.")
    parser.add_argument("--force", action="store_true", help="Remove item output directories before rerunning Rosetta.")
    parser.add_argument("--keep-pdbs", action="store_true", help="Keep Rosetta output PDB files instead of archiving them into structures.tar.gz.")
    parser.add_argument("--dry-run", action="store_true", help="Print Rosetta commands without running Docker.")
    parser.add_argument("--constant-seed", action="store_true", help="Pass Rosetta -constant_seed and -jran for deterministic reruns.")
    parser.add_argument("--seed-offset", type=int, default=1111, help="Base random seed used with --constant-seed.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = resolve_cli_path(args.output_dir)
    reference_fasta = resolve_cli_path(args.reference)
    haplotype_fasta = resolve_haplotype_fasta(args.input_fasta, reference_fasta)
    config = RunConfig(
        haplotype_fasta=haplotype_fasta,
        reference_fasta=reference_fasta,
        xml_script=resolve_cli_path(args.xml_script),
        output_dir=output_dir,
        rosetta_image=args.rosetta_image,
        rosetta_bin=args.rosetta_bin,
        num_backrub_trials=args.num_backrub_trials,
        num_replicates=args.num_replicates,
        max_minimization_iter=args.max_minimization_iter,
        abs_score_convergence_thresh=args.abs_score_convergence_thresh,
        backrub_trajectory_stride=args.backrub_trajectory_stride,
        max_workers=max(1, args.max_workers),
        min_vp1_position=CRYSTAL_MIN_VP1_POSITION,
        missing_pdb_residues=set(CRYSTAL_MISSING_PDB_RESIDUES),
        redo_excel=args.redo_excel,
        force=args.force,
        compress_pdbs=not args.keep_pdbs,
        dry_run=args.dry_run,
        constant_seed=args.constant_seed,
        seed_offset=args.seed_offset,
    )

    configure_logging(config.output_dir)
    prep_dir = resolve_cli_path(args.prep_dir)
    targets = build_targets(prep_dir, args.targets)

    logging.info("Repository root: %s", REPO_ROOT)
    logging.info("Selected targets: %s", ", ".join(target.name for target in targets))
    logging.info("Haplotype FASTA: %s", config.haplotype_fasta)
    logging.info("Reference FASTA: %s", config.reference_fasta)
    validate_inputs(targets, config)

    map_file_exists = (config.output_dir / "Haplotype_Mutation_Map.xlsx").exists()
    needs_mafft = args.legacy or not map_file_exists or config.force
    check_external_tools(config, needs_mafft=needs_mafft)

    ref_record = load_reference_record(config.reference_fasta)
    ht_records = load_haplotype_records(config.haplotype_fasta)
    logging.info("Found %s haplotypes", len(ht_records))

    run_global_wt_baseline(targets, config, parallel_replicates=args.dedicated_cores and not args.legacy)
    if args.legacy:
        run_legacy_pipeline(ht_records, ref_record, targets, config, serial=args.serial)
    else:
        run_optimized_pipeline(ht_records, ref_record, targets, config, dedicated_cores=args.dedicated_cores)


if __name__ == "__main__":
    main()