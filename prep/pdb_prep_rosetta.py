import glob
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from statistics import mean

from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.Polypeptide import protein_letters_3to1

try:
    from Bio import pairwise2
except ImportError:
    print("Error: Bio.pairwise2 is required. Install biopython with pairwise2 support.")
    sys.exit(1)


ROSETTA_IMAGE = "rosettacommons/rosetta:ml-408"
MOUNT_POINT = "/data"
WORKDIR = os.path.dirname(os.path.abspath(__file__))

MUTATION_FILE = "Sabin_rev_mutations.txt"
REFERENCE_FASTA = "Sabin_2_VP1.fasta"
NSTRUCT = 20
CPU_THREADS = max(1, os.cpu_count() or 1)

TARGETS = {
    "3epf": {
        "candidates": ["3epf.pdb", "3EPF.pdb"],
        "apply_mutations": True,
        "keep_chains": ["A", "R"],
        "interface_partner_chains": ["R"],
        "final_output": "VP1-receptor_final.pdb",
    },
    "8e8y": {
        "candidates": ["8e8y.pdb", "8E8Y.pdb", "8e8s.pdb", "8E8S.pdb"],
        "apply_mutations": True,
        "keep_chains": ["A", "H", "L"],
        "interface_partner_chains": ["H", "L"],
        "final_output": "VP1-9H2_final.pdb",
    },
    "9ocl": {
        "candidates": ["9ocl.pdb", "9OCL.pdb"],
        "apply_mutations": True,
        "keep_chains": ["A", "H", "L"],
        "interface_partner_chains": ["H", "L"],
        "final_output": "VP1-10D2_final.pdb",
    },
}

CHAIN_MAP = {"1": "A", "2": "B", "3": "C", "4": "D"}
VP1_CHAINS = ["A"]


class CleanSelect(Select):
    def __init__(self, keep_hetero_resnames=None, keep_chains=None):
        self.keep_hetero_resnames = set(keep_hetero_resnames or [])
        self.keep_chains = set(keep_chains) if keep_chains else None

    def accept_chain(self, chain):
        if self.keep_chains is None:
            return True
        return chain.id in self.keep_chains

    def accept_residue(self, residue):
        hetflag = residue.id[0]
        if hetflag == " ":
            return True
        if hetflag == "W":
            return False
        return residue.get_resname().strip() in self.keep_hetero_resnames

    def accept_atom(self, atom):
        altloc = atom.get_altloc()
        if altloc in (" ", "A"):
            return True
        return False


def resolve_input_file(candidates):
    for name in candidates:
        if os.path.exists(name):
            return name
    return None


def safe_three_to_one(resname):
    return protein_letters_3to1.get(resname.upper(), "X")


def parse_fasta_sequence(path):
    seq = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq.append(line)
    return "".join(seq)


def parse_mutations(path):
    mutations = []
    pattern = re.compile(r"^([A-Z])(\d+)([A-Z])$")
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            match = pattern.match(line)
            if not match:
                print(f"Skipping invalid mutation line: {line}")
                continue
            old_aa = match.group(1)
            seq_pos = int(match.group(2))
            new_aa = match.group(3)
            mutations.append((old_aa, seq_pos, new_aa))
    return mutations


def get_chain_sequence_and_residues(pdb_file, chain_id):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_file)
    model = next(structure.get_models())
    if chain_id not in model:
        raise KeyError(f"Chain {chain_id} not found in {pdb_file}")
    chain = model[chain_id]

    sequence = []
    residues = []
    for residue in chain:
        if residue.id[0] != " ":
            continue
        sequence.append(safe_three_to_one(residue.get_resname()))
        residues.append(residue)
    return "".join(sequence), residues


def map_mutations_to_chain_a(clean_pdb, reference_fasta, seq_mutations):
    ref_seq = parse_fasta_sequence(reference_fasta)
    pdb_seq, pdb_residues = get_chain_sequence_and_residues(clean_pdb, "A")

    if not ref_seq or not pdb_seq:
        raise RuntimeError(f"Empty sequence for mapping in {clean_pdb}")

    aln = pairwise2.align.globalms(
        ref_seq,
        pdb_seq,
        2.0,
        -1.0,
        -10.0,
        -0.5,
        one_alignment_only=True,
    )[0]
    aln_ref = aln.seqA
    aln_pdb = aln.seqB

    ref_index = -1
    pdb_index = -1
    ref_to_pdb = {}

    for char_ref, char_pdb in zip(aln_ref, aln_pdb):
        if char_ref != "-":
            ref_index += 1
        if char_pdb != "-":
            pdb_index += 1
        if char_ref != "-" and char_pdb != "-":
            ref_to_pdb[ref_index] = pdb_residues[pdb_index]

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("direct", clean_pdb)
    chain_a = next(structure.get_models())["A"]

    mapped = []
    for old_aa, seq_pos, new_aa in seq_mutations:
        ref0 = seq_pos - 1
        residue = ref_to_pdb.get(ref0)
        if residue is None:
            direct_key = (" ", seq_pos, " ")
            if direct_key in chain_a:
                residue = chain_a[direct_key]
                print(f"Fallback direct-number mapping for position {seq_pos} in {clean_pdb}")
            else:
                print(f"Warning: sequence position {seq_pos} could not be mapped to chain A in {clean_pdb}")
                continue
        pdb_one = safe_three_to_one(residue.get_resname())
        if pdb_one != old_aa:
            print(
                f"Warning: expected {old_aa}{seq_pos} but mapped residue is {pdb_one}{residue.id[1]}{residue.id[2].strip()}"
            )
        pdb_resid = f"{residue.id[1]}{residue.id[2].strip()}"
        mapped.append((pdb_resid, new_aa, residue))
        print(f"Mapped {old_aa}{seq_pos}{new_aa} -> A {pdb_resid} ({residue.get_resname()})")

    return mapped


def generate_resfile(mapped_mutations, out_resfile, chain_id="A"):
    seen = set()
    with open(out_resfile, "w", encoding="utf-8") as handle:
        handle.write("NATAA\n")
        handle.write("start\n")
        for pdb_resid, new_aa, _ in mapped_mutations:
            key = (pdb_resid, chain_id)
            if key in seen:
                continue
            handle.write(f"{pdb_resid} {chain_id} PIKAA {new_aa}\n")
            seen.add(key)
    print(f"Wrote {out_resfile} with {len(seen)} unique mutation entries")


def clean_and_prepare_pdb(in_pdb, out_pdb, keep_hetero=None, keep_chains=None):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", in_pdb)
    model = next(structure.get_models())

    for chain in list(model):
        if chain.id in CHAIN_MAP:
            chain.id = CHAIN_MAP[chain.id]

    io = PDBIO()
    io.set_structure(structure)
    io.save(out_pdb, CleanSelect(keep_hetero_resnames=keep_hetero or [], keep_chains=keep_chains))
    print(f"[OK] cleaned PDB: {out_pdb}")


def write_relax_xml(xml_path):
    xml = """<ROSETTASCRIPTS>
    <SCOREFXNS>
        <ScoreFunction name="ref15_cart" weights="ref2015_cart">
            <Reweight scoretype="coordinate_constraint" weight="1.0"/>
            <Reweight scoretype="cart_bonded" weight="0.5"/>
            <Reweight scoretype="pro_close" weight="0.0"/>
        </ScoreFunction>
    </SCOREFXNS>
    <TASKOPERATIONS>
        <ReadResfile name="rrf" filename="%%resfile%%"/>
        <InitializeFromCommandline name="init"/>
        <IncludeCurrent name="current"/>
    </TASKOPERATIONS>
    <MOVERS>
        <FastRelax name="relax" scorefxn="ref15_cart" task_operations="rrf,init,current" cartesian="1" bondangle="1" disable_design="false"/>
    </MOVERS>
    <PROTOCOLS>
        <Add mover="relax"/>
    </PROTOCOLS>
</ROSETTASCRIPTS>
"""
    with open(xml_path, "w", encoding="utf-8") as handle:
        handle.write(xml)


def run_docker(command_args, env_vars=None):
    user_args = []
    if os.name == "posix":
        try:
            user_args = ["--user", f"{os.getuid()}:{os.getgid()}"]
        except AttributeError:
            pass

    env_args = []
    for key, value in (env_vars or {}).items():
        env_args.extend(["-e", f"{key}={value}"])

    cmd = [
        "docker",
        "run",
        "--rm",
        "-v",
        f"{WORKDIR}:{MOUNT_POINT}",
        "-w",
        MOUNT_POINT,
    ] + env_args + user_args + [ROSETTA_IMAGE] + command_args

    print("Running:", " ".join(command_args))
    result = subprocess.run(cmd)
    return result.returncode == 0


def find_docker_binary(pattern):
    cmd = ["docker", "run", "--rm", ROSETTA_IMAGE, "find", "/usr/local/bin", "/usr/bin", "-name", pattern]
    try:
        output = subprocess.check_output(cmd, text=True).strip().splitlines()
    except subprocess.CalledProcessError:
        return None
    for line in output:
        line = line.strip()
        if line:
            return line
    return None


def parse_scorefile(scorefile_path):
    if not os.path.exists(scorefile_path):
        return []
    rows = []
    headers = None
    with open(scorefile_path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.strip()
            if not line.startswith("SCORE:"):
                continue
            tokens = line.split()
            if len(tokens) < 3:
                continue
            if tokens[1] == "score":
                headers = tokens[1:]
                continue
            if headers is None:
                continue
            values = tokens[1:]
            if len(values) != len(headers):
                continue
            row = dict(zip(headers, values))
            rows.append(row)
    return rows


def collect_relaxed_models(tag):
    return sorted(glob.glob(f"{tag}_clean_relaxed_*.pdb"))


def ca_map(structure):
    mapping = {}
    model = next(structure.get_models())
    for chain in model:
        for residue in chain:
            if residue.id[0] != " ":
                continue
            if "CA" in residue:
                key = (chain.id, residue.id[1], residue.id[2].strip())
                mapping[key] = residue["CA"].coord
    return mapping


def ca_rmsd_between(pdb_a, pdb_b):
    parser = PDBParser(QUIET=True)
    s1 = parser.get_structure("a", pdb_a)
    s2 = parser.get_structure("b", pdb_b)
    m1 = ca_map(s1)
    m2 = ca_map(s2)
    common = set(m1.keys()) & set(m2.keys())
    if not common:
        return None
    sq = 0.0
    n = 0
    for key in common:
        d = m1[key] - m2[key]
        sq += float(d.dot(d))
        n += 1
    return (sq / n) ** 0.5


def pick_centroid_model(model_paths):
    if len(model_paths) <= 1:
        return model_paths[0] if model_paths else None
    mean_rmsd = {}
    for i, path_i in enumerate(model_paths):
        values = []
        for j, path_j in enumerate(model_paths):
            if i == j:
                continue
            rmsd = ca_rmsd_between(path_i, path_j)
            if rmsd is not None:
                values.append(rmsd)
        mean_rmsd[path_i] = mean(values) if values else float("inf")
    return min(mean_rmsd, key=mean_rmsd.get)


def compute_interface_contacts(pdb_path, group_a, group_b, cutoff=5.0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("x", pdb_path)
    model = next(structure.get_models())

    atoms_a = []
    atoms_b = []

    for chain in model:
        if chain.id in group_a:
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                if "CB" in residue:
                    atoms_a.append(residue["CB"].coord)
                elif "CA" in residue:
                    atoms_a.append(residue["CA"].coord)
        elif chain.id in group_b:
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                if "CB" in residue:
                    atoms_b.append(residue["CB"].coord)
                elif "CA" in residue:
                    atoms_b.append(residue["CA"].coord)

    if not atoms_a or not atoms_b:
        return 0

    cutoff_sq = cutoff * cutoff
    contacts = 0
    for a in atoms_a:
        for b in atoms_b:
            diff = a - b
            if float(diff.dot(diff)) <= cutoff_sq:
                contacts += 1
    return contacts


def missing_mutations_in_model(pdb_path, mapped_mutations, chain_id="A"):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("mutchk", pdb_path)
    model = next(structure.get_models())
    if chain_id not in model:
        return [("CHAIN", chain_id, "missing")]
    chain = model[chain_id]

    missing = []
    for _, expected_new_aa, residue in mapped_mutations:
        resseq = residue.id[1]
        icode = residue.id[2].strip()
        lookup = (" ", resseq, icode if icode else " ")
        if lookup not in chain:
            missing.append((f"{resseq}{icode}", expected_new_aa, "missing"))
            continue
        observed = safe_three_to_one(chain[lookup].get_resname())
        if observed != expected_new_aa:
            missing.append((f"{resseq}{icode}", expected_new_aa, observed))
    return missing


def run_rosetta_with_threads(cmd, score_tag, cpu_threads):
    threading_env = {"OMP_NUM_THREADS": str(cpu_threads)}
    threaded_cmd = cmd + [
        "-multithreading:total_threads",
        str(cpu_threads),
        "-multithreading:interaction_graph_threads",
        str(max(1, cpu_threads // 2)),
    ]

    print(f"Attempting threaded Rosetta run for {score_tag} with {cpu_threads} CPU threads")
    ok = run_docker(threaded_cmd, env_vars=threading_env)
    if ok:
        return True

    print(f"Threaded Rosetta flags failed for {score_tag}; retrying without multithreading flags")
    return run_docker(cmd, env_vars=threading_env)


def cleanup_target_outputs(tag, scorefile, final_output):
    for old_model in glob.glob(f"{tag}_clean_relaxed_*.pdb"):
        try:
            os.remove(old_model)
        except OSError:
            pass
    for path in [scorefile, final_output]:
        if os.path.exists(path):
            try:
                os.remove(path)
            except OSError:
                pass


def process_target(tag, config, source_pdb, rosetta_scripts, seq_mutations, xml_path, cpu_threads):
    clean_pdb = f"{tag}_clean.pdb"
    resfile = f"{tag}.resfile"
    scorefile = f"{tag}_relax.sc"
    final_output = config.get("final_output", f"{tag}_final.pdb")

    print(f"\n=== {tag} ({source_pdb}) ===")
    cleanup_target_outputs(tag, scorefile, final_output)

    clean_and_prepare_pdb(
        source_pdb,
        clean_pdb,
        keep_hetero=["ZN", "MG", "CA", "NA", "K", "CL"],
        keep_chains=config["keep_chains"],
    )

    mapped = []
    if config["apply_mutations"]:
        mapped = map_mutations_to_chain_a(clean_pdb, REFERENCE_FASTA, seq_mutations)
    generate_resfile(mapped, resfile, chain_id="A")

    cmd = [
        rosetta_scripts,
        "-parser:protocol",
        xml_path,
        "-parser:script_vars",
        f"resfile={resfile}",
        "-s",
        clean_pdb,
        "-relax:constrain_relax_to_start_coords",
        "-relax:coord_constrain_sidechains",
        "-relax:cartesian",
        "-relax:respect_resfile",
        "true",
        "-nstruct",
        str(NSTRUCT),
        "-out:suffix",
        "_relaxed",
        "-out:file:scorefile",
        scorefile,
        "-packing:repack_only",
        "false",
        "-use_input_sc",
        "-ignore_zero_occupancy",
        "false",
        "-ex1",
        "-ex2aro",
    ]

    ok = run_rosetta_with_threads(cmd, tag, cpu_threads)
    if not ok:
        return tag, False, "Rosetta failed"

    models = collect_relaxed_models(tag)
    if not models:
        return tag, False, "No relaxed models found"

    scores = parse_scorefile(scorefile)
    best_by_score = None
    best_by_score_pdb = None
    if scores and "description" in scores[0] and "total_score" in scores[0]:
        try:
            row = min(scores, key=lambda r: float(r["total_score"]))
            best_by_score = row["description"]
            best_by_score_pdb = best_by_score if best_by_score.endswith(".pdb") else f"{best_by_score}.pdb"
        except ValueError:
            best_by_score = None
            best_by_score_pdb = None

    candidate_models = list(models)
    if mapped:
        valid_models = []
        for model_path in models:
            missing = missing_mutations_in_model(model_path, mapped, chain_id="A")
            if not missing:
                valid_models.append(model_path)
        if not valid_models:
            return tag, False, "No relaxed models contain all requested mutations"
        if len(valid_models) < len(models):
            print(f"Mutation filter kept {len(valid_models)}/{len(models)} models for {tag}")
        candidate_models = valid_models

    if best_by_score_pdb and best_by_score_pdb in candidate_models:
        final_model = best_by_score_pdb
    else:
        final_model = pick_centroid_model(candidate_models) or candidate_models[0]

    centroid_model = pick_centroid_model(candidate_models)
    first_model = candidate_models[0]
    first_rmsd = ca_rmsd_between(clean_pdb, first_model)

    print(f"Generated models ({tag}): {len(models)}")
    if best_by_score:
        print(f"Lowest-energy model ({tag}): {best_by_score}")
    if centroid_model:
        print(f"Centroid model ({tag}): {centroid_model}")
    print(f"Selected final model ({tag}): {final_model}")
    if first_rmsd is not None:
        print(f"RMSD(clean vs first relaxed, {tag}): {first_rmsd:.3f} A")

    contacts = compute_interface_contacts(
        final_model,
        group_a=VP1_CHAINS,
        group_b=config["interface_partner_chains"],
        cutoff=5.0,
    )
    print(f"Interface contacts (VP1 vs partner, {tag}) at 5.0 A: {contacts}")

    shutil.copyfile(final_model, final_output)
    print(f"Exported final complex ({tag}): {final_output}")
    return tag, True, final_output


def main():
    if not os.path.exists(MUTATION_FILE):
        print(f"Missing mutation file: {MUTATION_FILE}")
        sys.exit(1)
    if not os.path.exists(REFERENCE_FASTA):
        print(f"Missing reference FASTA: {REFERENCE_FASTA}")
        sys.exit(1)

    rosetta_scripts = find_docker_binary("rosetta_scripts.*linuxgccrelease")
    if not rosetta_scripts:
        print("Could not find rosetta_scripts binary in Docker image.")
        sys.exit(1)

    seq_mutations = parse_mutations(MUTATION_FILE)
    xml_path = "relax_mutations.xml"
    write_relax_xml(xml_path)

    print("Preparing VP1-partner complexes (without VP2/VP3/VP4) and running constrained Cartesian relax...")
    print(f"rosetta_scripts: {rosetta_scripts}")
    available_targets = []
    for tag, config in TARGETS.items():
        source_pdb = resolve_input_file(config["candidates"])
        if not source_pdb:
            print(f"Skipping {tag}: none of {config['candidates']} found")
            continue
        available_targets.append((tag, config, source_pdb))

    if not available_targets:
        print("No valid input PDB files found.")
        sys.exit(1)

    max_workers = min(3, len(available_targets))
    threads_per_target = max(1, CPU_THREADS // max_workers)
    print(f"Configured CPU threads: total={CPU_THREADS}, parallel_targets={max_workers}, per_target={threads_per_target}")

    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_tag = {}
        for tag, config, source_pdb in available_targets:
            future = executor.submit(
                process_target,
                tag,
                config,
                source_pdb,
                rosetta_scripts,
                seq_mutations,
                xml_path,
                threads_per_target,
            )
            future_to_tag[future] = tag

        for future in as_completed(future_to_tag):
            tag = future_to_tag[future]
            try:
                result_tag, ok, message = future.result()
            except Exception as exc:
                print(f"[{tag}] Failed with exception: {exc}")
                results.append((tag, False, str(exc)))
                continue
            status = "OK" if ok else "FAILED"
            print(f"[{result_tag}] {status}: {message}")
            results.append((result_tag, ok, message))

    failed = [tag for tag, ok, _ in results if not ok]
    if failed:
        print(f"\nCompleted with failures for: {', '.join(sorted(failed))}")
    else:
        print("\nAll targets completed successfully.")

    print("\nDone. VP1-partner complexes prepared for downstream flex ddG.")


if __name__ == "__main__":
    os.chdir(WORKDIR)
    main()
