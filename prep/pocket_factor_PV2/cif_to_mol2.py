import shlex
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
cif_path = SCRIPT_DIR / 'PLM.cif'
mol2_path = SCRIPT_DIR / 'PLM_from_cif.mol2'

lines = cif_path.read_text(encoding='utf-8').splitlines()

def parse_loop(target_prefix):
    i = 0
    while i < len(lines):
        if lines[i].strip() != 'loop_':
            i += 1
            continue
        i += 1
        headers = []
        while i < len(lines) and lines[i].startswith('_'):
            headers.append(lines[i].strip())
            i += 1
        if not headers:
            continue
        if not any(h.startswith(target_prefix) for h in headers):
            while i < len(lines) and lines[i].strip() != '#':
                i += 1
            i += 1
            continue
        rows = []
        while i < len(lines):
            s = lines[i].strip()
            if not s or s == '#':
                if s == '#':
                    i += 1
                break
            if s.startswith('loop_') or s.startswith('_'):
                break
            rows.append(shlex.split(lines[i]))
            i += 1
        return headers, rows
    raise RuntimeError(f'Loop not found for prefix {target_prefix}')

atom_headers, atom_rows = parse_loop('_chem_comp_atom.')
bond_headers, bond_rows = parse_loop('_chem_comp_bond.')

atom_cols = {h: idx for idx, h in enumerate(atom_headers)}
bond_cols = {h: idx for idx, h in enumerate(bond_headers)}

atoms = []
for row in atom_rows:
    atom_id = row[atom_cols['_chem_comp_atom.atom_id']]
    elem = row[atom_cols['_chem_comp_atom.type_symbol']]
    x = float(row[atom_cols['_chem_comp_atom.model_Cartn_x']])
    y = float(row[atom_cols['_chem_comp_atom.model_Cartn_y']])
    z = float(row[atom_cols['_chem_comp_atom.model_Cartn_z']])
    atoms.append((atom_id, elem, x, y, z))

atom_index = {name: i + 1 for i, (name, *_rest) in enumerate(atoms)}

bonds = []
for row in bond_rows:
    a1 = row[bond_cols['_chem_comp_bond.atom_id_1']]
    a2 = row[bond_cols['_chem_comp_bond.atom_id_2']]
    order = row[bond_cols['_chem_comp_bond.value_order']]
    aromatic = row[bond_cols['_chem_comp_bond.pdbx_aromatic_flag']]
    if aromatic.upper() == 'Y':
        btype = 'ar'
    else:
        btype = {'SING': '1', 'DOUB': '2', 'TRIP': '3'}.get(order.upper(), '1')
    bonds.append((a1, a2, btype))

with mol2_path.open('w', encoding='utf-8', newline='\n') as out:
    out.write('@<TRIPOS>MOLECULE\n')
    out.write('PLM\n')
    out.write(f'{len(atoms):6d} {len(bonds):6d} {1:6d}\n')
    out.write('SMALL\n')
    out.write('USER_CHARGES\n\n')
    out.write('Converted from PLM.cif\n')
    out.write('@<TRIPOS>ATOM\n')
    for i, (name, elem, x, y, z) in enumerate(atoms, start=1):
        out.write(f'{i:8d} {name:<8s} {x:9.4f} {y:9.4f} {z:9.4f} {elem:<5s} {1:5d} {"PLM":<8s} {0.0:9.4f}\n')
    out.write('@<TRIPOS>BOND\n')
    for i, (a1, a2, btype) in enumerate(bonds, start=1):
        out.write(f'{i:8d} {atom_index[a1]:8d} {atom_index[a2]:8d} {btype}\n')

print(f'Wrote {mol2_path}')