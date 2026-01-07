import json
from pathlib import Path

nb_path = Path('notebooks/kaggle/md_simulation_264THM_PPARG.ipynb')
with open(nb_path, 'r', encoding='utf-8') as f:
    nb = json.load(f)

# Find and fix the ITP extraction cell
for i, cell in enumerate(nb['cells']):
    if cell['cell_type'] == 'code':
        source = ''.join(cell['source'])
        
        # Fix the ligand.itp extraction - separate atomtypes and moleculetype
        if 'Extract ligand ITP from generated top file' in source:
            new_source = '''%%bash
set -e
cd /kaggle/working/264THM_PPARG/topol

# Extract ligand ITP from generated top file
python << 'EOF'
with open('ligand.top', 'r') as f:
    lines = f.readlines()

# Find atomtypes and moleculetype sections
in_atomtypes = False
in_moleculetype = False
atomtypes_lines = []
moleculetype_lines = []

for line in lines:
    if '[ atomtypes ]' in line:
        in_atomtypes = True
        in_moleculetype = False
        atomtypes_lines.append(line)
        continue
    if '[ moleculetype ]' in line:
        in_atomtypes = False
        in_moleculetype = True
    if in_atomtypes:
        atomtypes_lines.append(line)
    if in_moleculetype:
        moleculetype_lines.append(line)

# Write atomtypes to separate file
with open('ligand_atomtypes.itp', 'w') as f:
    f.write('; Ligand atomtypes from AmberTools\\n')
    f.writelines(atomtypes_lines)

# Write moleculetype to ligand.itp (without atomtypes)
with open('ligand.itp', 'w') as f:
    f.write('; Ligand topology from AmberTools\\n')
    f.writelines(moleculetype_lines)

print('ligand_atomtypes.itp and ligand.itp created!')
EOF'''
            cell['source'] = [new_source]
            print(f'Fixed ITP extraction cell {i}')

        # Fix the topology update to include atomtypes in correct position
        if "insert_pos = topol.find('[ system ]')" in source:
            new_source = '''os.chdir(WORK_DIR / 'topol')

with open('protein.gro', 'r') as f:
    protein_lines = f.readlines()
with open('ligand.gro', 'r') as f:
    ligand_lines = f.readlines()

protein_atoms = protein_lines[2:-1]
ligand_atoms = ligand_lines[2:-1]
box = protein_lines[-1]
total_atoms = len(protein_atoms) + len(ligand_atoms)

with open('complex.gro', 'w') as f:
    f.write(f"{CONFIG['complex_name']} complex\\n")
    f.write(f' {total_atoms}\\n')
    f.writelines(protein_atoms)
    f.writelines(ligand_atoms)
    f.write(box)

print(f'Complex: {total_atoms} atoms')

# Update topology - atomtypes must come AFTER forcefield but BEFORE water
with open('topol.top', 'r') as f:
    topol = f.read()

# Find position after forcefield include and before position restraints
ff_include_end = topol.find('#include "amber99sb-ildn.ff/tip3p.itp"')
if ff_include_end == -1:
    # Try alternative pattern
    ff_include_end = topol.find('#include "./amber99sb-ildn.ff')
if ff_include_end == -1:
    # Fallback: find after all forcefield includes
    import re
    matches = list(re.finditer(r'#include.*amber99sb.*\\.itp.*\\n', topol))
    if matches:
        ff_include_end = matches[-1].end()

# Insert atomtypes after forcefield
if ff_include_end > 0:
    topol = topol[:ff_include_end] + '\\n#include "ligand_atomtypes.itp"\\n' + topol[ff_include_end:]

# Add ligand.itp before [ system ]
insert_pos = topol.find('[ system ]')
if insert_pos > 0:
    topol = topol[:insert_pos] + '#include "ligand.itp"\\n\\n' + topol[insert_pos:]

# Add LIG to molecules
topol += '\\nLIG     1\\n'

with open('topol.top', 'w') as f:
    f.write(topol)

print('Topology updated with correct atomtypes order!')
os.chdir(WORK_DIR)'''
            cell['source'] = [new_source]
            print(f'Fixed topology update cell {i}')

with open(nb_path, 'w', encoding='utf-8') as f:
    json.dump(nb, f, indent=4)

print('Notebook fixed!')
