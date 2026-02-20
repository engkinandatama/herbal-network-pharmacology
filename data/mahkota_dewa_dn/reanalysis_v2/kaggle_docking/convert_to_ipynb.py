"""Convert moldock_vina.py to .ipynb notebook"""
import json, re

py_file = "data/mahkota_dewa_dn/reanalysis_v2/kaggle_docking/moldock_vina.py"
ipynb_file = "data/mahkota_dewa_dn/reanalysis_v2/kaggle_docking/moldock_vina.ipynb"

with open(py_file, encoding="utf-8") as f:
    content = f.read()

# Split on # %% markers
parts = re.split(r'^# %%(.*)$', content, flags=re.MULTILINE)

cells = []

# First part is the docstring header -> markdown cell
header = parts[0].strip()
if header.startswith('"""') and header.endswith('"""'):
    header_text = header.strip('"""').strip()
    cells.append({
        "cell_type": "markdown",
        "metadata": {},
        "source": header_text.splitlines(keepends=True)
    })

# Process remaining parts (marker_suffix, code) pairs
i = 1
while i < len(parts):
    marker = parts[i].strip() if i < len(parts) else ""
    code = parts[i+1] if i+1 < len(parts) else ""
    i += 2
    
    if marker.startswith("[markdown]"):
        # Markdown cell - extract the # comment lines
        lines = code.strip().splitlines(keepends=True)
        md_lines = []
        for line in lines:
            if line.startswith("# "):
                md_lines.append(line[2:])
            elif line.strip() == "#":
                md_lines.append("\n")
            else:
                md_lines.append(line)
        cells.append({
            "cell_type": "markdown",
            "metadata": {},
            "source": md_lines
        })
    else:
        # Code cell
        code = code.strip()
        if code:
            cells.append({
                "cell_type": "code",
                "metadata": {},
                "source": [line + "\n" for line in code.splitlines()],
                "outputs": [],
                "execution_count": None
            })

notebook = {
    "nbformat": 4,
    "nbformat_minor": 4,
    "metadata": {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3"
        },
        "language_info": {
            "name": "python",
            "version": "3.10.0"
        }
    },
    "cells": cells
}

with open(ipynb_file, "w", encoding="utf-8") as f:
    json.dump(notebook, f, indent=1)

print(f"Created {ipynb_file}")
print(f"Total cells: {len(cells)}")
for i, c in enumerate(cells):
    ctype = c["cell_type"]
    preview = "".join(c["source"])[:60].replace("\n", " ")
    print(f"  [{i}] {ctype:>10s}: {preview}...")
