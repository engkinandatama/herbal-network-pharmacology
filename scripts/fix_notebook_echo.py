import json
from pathlib import Path

nb_path = Path('notebooks/kaggle/md_simulation_264THM_PPARG.ipynb')
with open(nb_path, 'r', encoding='utf-8') as f:
    content = f.read()

# Fix all echo pipe patterns to use printf for more reliable piping
replacements = [
    ("echo '1' | gmx pdb2gmx", "printf '1\\n' | gmx pdb2gmx"),
    ("echo 'SOL' | gmx genion", "printf 'SOL\\n' | gmx genion"),
    ("echo '4 4' | gmx rms", "printf '4\\n4\\n' | gmx rms"),
    ("echo '4' | gmx rmsf", "printf '4\\n' | gmx rmsf"),
    ("echo '1' | gmx gyrate", "printf '1\\n' | gmx gyrate"),
]

for old, new in replacements:
    if old in content:
        content = content.replace(old, new)
        print(f"Fixed: {old}")

with open(nb_path, 'w', encoding='utf-8') as f:
    f.write(content)

print("All echo pipe patterns fixed!")
