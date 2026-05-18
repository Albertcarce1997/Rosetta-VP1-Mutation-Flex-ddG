@echo off
cd /d "%~dp0"

REM Download 8E8Y.pdb
curl -o 8E8Y.pdb https://files.rcsb.org/download/8E8Y.pdb

REM Download PLM.cif
curl -o PLM.cif https://files.rcsb.org/ligands/view/PLM.cif

REM Run Python script to convert CIF to MOL2
python cif_to_mol2.py

REM Use Rosetta docker to generate params file
echo Running Rosetta molfile_to_params inside Docker...
docker run --rm -v "%CD%:/data" -w /data rosettacommons/rosetta:ml-408 ^
    python3 /rosetta/main/source/scripts/python/public/molfile_to_params.py ^
    -n PLM -p PLM --conformers-in-one-file PLM_from_cif.mol2

echo Done! You can now run pdb_prep_rosetta.py
