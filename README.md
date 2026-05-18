# Rosetta VP1 Mutation Flex ddG Pipeline

This repository contains the scripts used to prepare poliovirus VP1 complexes and evaluate the effect of VP1 haplotype mutations with a Rosetta flex ddG workflow. The analysis compares mutated haplotypes against the Sabin 2 VP1 reference in four structural contexts:

- VP1 bound to the poliovirus receptor CD155
- VP1 bound to antibody 9H2
- VP1 bound to antibody 10D2
- VP1 retaining the pocket factor palmitic acid, PLM

The pipeline has two stages. First, the preparation scripts extract the relevant VP1-containing complexes from the source PDB files, introduce the Sabin 2 reversion mutations where needed, retain the desired binding partner, and perform constrained Rosetta relaxation. The pocket-factor input `8E8Y.pdb` is already Sabin 2, so that branch keeps PLM and relaxes the structure without introducing mutations. Second, the flex ddG script aligns haplotype sequences to the Sabin 2 reference, converts sequence differences to Rosetta resfiles, runs Rosetta backrub/flex ddG calculations, and summarizes the score files.

## Repository Layout

```text
.
├── rosetta_env.yml
├── prep/
│   ├── pdb_prep_rosetta.py
│   ├── relax_mutations.xml
│   ├── Sabin_2_VP1.fasta
│   ├── Sabin_rev_mutations.txt
│   ├── 3epf.pdb
│   ├── 8E8Y.pdb
│   ├── 9OCL.pdb
│   └── pocket_factor_PV2/
│       ├── pdb_prep_rosetta.py
│       ├── cif_to_mol2.py
│       ├── prepare_target.bat
│       ├── PLM.cif
│       ├── PLM.mol2
│       ├── PLM.params
│       └── 8E8Y.pdb
└── flexddg/
    ├── run_flex_ddg_pipeline.py
    ├── Sabin_2_VP1.fasta
    ├── Top_sequences_by_HT_20seq_new.fasta
    └── flex_ddg_pipeline/
        └── ddG-backrub.xml
```

The main flex ddG entry point is `flexddg/run_flex_ddg_pipeline.py`. Use `--targets pocket_factor` to run only the PLM pocket-factor analysis.

## Requirements

The Python environment is intentionally small:

- Python 3.10 or 3.11
- Biopython
- pandas
- openpyxl
- MAFFT

Rosetta itself is run through Docker using the `rosettacommons/rosetta:ml-408` image. You must have Docker installed, running, and able to pull or access that image.

Create the conda environment with:

```bash
mamba env create -f rosetta_env.yml
mamba activate rosetta-vp1-flexddg
```

On Windows, the easiest route is usually WSL2 plus Docker Desktop integration, because the environment uses Bioconda packages such as MAFFT. Native Windows can work if MAFFT is available on `PATH` and Docker can mount the repository drive.

Check the external tools before launching long jobs:

```bash
mafft --version
docker run --rm rosettacommons/rosetta:ml-408 /usr/local/bin/rosetta_scripts -help
```

If the Rosetta binary path differs in your image, pass it with `--rosetta-bin` when running the flex ddG script.

## Stage 1: Prepare Structures

The preparation scripts are now script-relative, so they can be launched from the repository root or from their own folders.

### Receptor and Antibody Complexes

Run:

```bash
python prep/pdb_prep_rosetta.py
```

Inputs:

- `prep/3epf.pdb`: VP1-receptor source structure
- `prep/8E8Y.pdb`: VP1-9H2 source structure
- `prep/9OCL.pdb`: VP1-10D2 source structure
- `prep/Sabin_2_VP1.fasta`: Sabin 2 VP1 reference sequence
- `prep/Sabin_rev_mutations.txt`: mutations used to generate the Sabin 2 version

Main operations:

- Keeps only the VP1 chain and the relevant partner chains
- Normalizes numeric chain identifiers to alphabetic chain IDs when needed
- Maps mutation positions from the Sabin 2 VP1 FASTA to PDB chain A using sequence alignment
- Writes Rosetta resfiles for the Sabin reversion mutations
- Runs constrained Cartesian relaxation in the Rosetta Docker image
- Chooses a final relaxed model using score and centroid/RMSD checks

Output:

- `prep/VP1-receptor_final.pdb`
- `prep/VP1-9H2_final.pdb`
- `prep/VP1-10D2_final.pdb`

### Pocket-Factor Complex

Run:

```bash
python prep/pocket_factor_PV2/pdb_prep_rosetta.py
```

Inputs:

- `prep/pocket_factor_PV2/8E8Y.pdb`: Sabin 2 VP1 source structure containing PLM
- `prep/pocket_factor_PV2/PLM.params`: Rosetta ligand params for palmitic acid

Main operations:

- Keeps VP1 chain A and the PLM hetero residue
- Skips mutation introduction because this `8E8Y.pdb` input is already Sabin 2
- Passes `PLM.params` to Rosetta with `-extra_res_fa`
- Relaxes the VP1-PLM structure
- Writes a final PDB and a second copy with PLM assigned to chain `X` for the downstream interface calculation

Expected outputs:

- `prep/pocket_factor_PV2/8E8Y_VP1_PLM_final.pdb`
- `prep/pocket_factor_PV2/8E8Y_VP1_PLM_final_ligX.pdb`

The repository includes `PLM.params`. If you need to regenerate it, run `prep/pocket_factor_PV2/prepare_target.bat` on Windows or adapt the same commands for Bash. That helper downloads `8E8Y.pdb` and `PLM.cif`, converts the CIF to MOL2, and runs Rosetta `molfile_to_params.py` inside Docker.

## Stage 2: Run Flex ddG

After the prepared PDB files exist, run the unified flex ddG script:

```bash
python flexddg/run_flex_ddg_pipeline.py --input flexddg/your_haplotypes.fasta --targets all --max-workers 8
```

Default inputs:

- Haplotypes: auto-detected from a single non-reference FASTA file in `flexddg/`, or set explicitly with `--input`
- Reference: `flexddg/Sabin_2_VP1.fasta`
- Prepared structures: files generated under `prep/`
- RosettaScripts XML: `flexddg/flex_ddg_pipeline/ddG-backrub.xml`
- Output directory: `flex_ddg_results/`

If `--input` is omitted, the script scans `flexddg/` for FASTA files with extensions `.fa`, `.faa`, `.fas`, or `.fasta`, excludes the configured reference FASTA and reference-like filenames, and uses the remaining file only when exactly one candidate is present. If multiple non-reference FASTA files are present, pass the desired haplotype file with `--input`. The older `--haplotypes` flag is still accepted as an alias for `--input`.

Target choices:

- `--targets all`: receptor, both antibodies, and pocket factor
- `--targets protein_complexes`: receptor plus both antibodies
- `--targets receptor`
- `--targets antibodies`
- `--targets 9H2`
- `--targets 10D2`
- `--targets pocket_factor`

Common runtime options:

- `--max-workers` controls concurrent Rosetta jobs or mutation-set workers. Choose a value that matches the available CPUs and memory
- `--input` selects the haplotype FASTA file. Use this whenever more than one non-reference FASTA file is present
- `--dedicated-cores` runs each replicate as a separate Rosetta job. This can use large machines efficiently, but it launches many Docker containers
- `--legacy` processes each haplotype independently. The default optimized mode deduplicates identical mutation sets and single mutations before running Rosetta
- `--redo-excel` regenerates Excel summaries from existing score files when possible
- `--force` removes item output directories and reruns Rosetta
- `--dry-run` prints the Docker/Rosetta commands without executing them
- `--keep-pdbs` leaves Rosetta output PDB files uncompressed. By default they are archived into `structures.tar.gz` under each output item

The manuscript-style protocol parameters are encoded as script defaults. Use `python flexddg/run_flex_ddg_pipeline.py --help` to inspect all advanced options, and change protocol-level settings only when intentionally running a different analysis.

## Mutation Handling

The flex ddG script aligns every haplotype to the Sabin 2 VP1 reference using MAFFT. It then records substitutions as `(wild-type residue, VP1 reference position, mutant residue)` tuples.

By default, the script uses the crystal-structure residue coverage rules:

- Mutations before VP1 position 27 are ignored
- Positions 96 to 101 are skipped because they are missing in the crystal structures
- Haplotypes requiring insertions or real deletions are skipped and reported as requiring docking/remodeling

For each valid haplotype, optimized mode creates:

- One mutation set representing all valid mutations in that haplotype
- One single-mutation job for each individual mutation

Identical mutation sets shared by several haplotypes are run only once, then mapped back to each haplotype in the final aggregate table.

## Outputs

The default output directory is `flex_ddg_results/`.

Key files and folders:

- `flex_ddg_results/flex_ddg_pipeline.log`: runtime log
- `flex_ddg_results/Haplotype_Mutation_Map.xlsx`: haplotype to mutation-set mapping
- `flex_ddg_results/WT_Baseline/WT_flex_ddg_results.xlsx`: Sabin 2 baseline results for each target
- `flex_ddg_results/all_mutations/Set_*/`: Rosetta outputs and summaries for deduplicated full haplotype mutation sets
- `flex_ddg_results/single_mutations/<mutation>/`: Rosetta outputs and summaries for individual mutations
- `flex_ddg_results/Optimized_ALL_HT_results.xlsx`: final aggregated table mapping set and single-mutation results back to haplotypes
- `flex_ddg_results/skipped_haplotypes_docking_required.txt`: haplotypes skipped because the alignment implied insertions or deletions requiring structural remodeling

The summary parser looks for Rosetta score columns in this order:

1. `ddG`
2. `ddg`
3. `dG_separated`
4. `dG_bind` or `dG_binding`
5. `I_sc`
6. `total_score`
7. `score`

The selected column is reported in the `metric` field. Interpret the sign and magnitude according to the Rosetta metric actually present in your score files.

## Reproducing the Full Workflow

From a fresh clone:

```bash
mamba env create -f rosetta_env.yml
mamba activate rosetta-vp1-flexddg

python prep/pdb_prep_rosetta.py
python prep/pocket_factor_PV2/pdb_prep_rosetta.py

python flexddg/run_flex_ddg_pipeline.py --input flexddg/your_haplotypes.fasta --targets all --max-workers 8
```

For a quick command check without launching Rosetta:

```bash
python flexddg/run_flex_ddg_pipeline.py --input flexddg/your_haplotypes.fasta --targets pocket_factor --dry-run
```

## Troubleshooting

If the flex ddG script reports missing prepared PDB files, run the preparation stage first. The default flex ddG paths point to the expected outputs in `prep/` and `prep/pocket_factor_PV2/`.

If Docker cannot mount the repository on Windows, ensure Docker Desktop has access to the drive containing the repository or run from WSL2.

If Rosetta reports an unknown residue for PLM, confirm that `prep/pocket_factor_PV2/PLM.params` exists and that the pocket-factor target is being run through the unified script, which passes `-extra_res_fa` automatically.

If MAFFT is missing, recreate the conda environment or install MAFFT separately and ensure it is available on `PATH`.
