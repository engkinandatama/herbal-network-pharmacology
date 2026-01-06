# Network Pharmacology Toolkit

A modular Python toolkit for network pharmacology analysis of herbal medicines. Designed for reproducible research with config-based workflows.

## Features

- **Compound Collection**: Fetch compounds from KNApSAcK, PubChem, and literature
- **Target Prediction**: Predict drug targets using SwissTargetPrediction
- **Disease Genes**: Retrieve disease-associated genes from GeneCards, DisGeNET, OMIM
- **Network Analysis**: Build and analyze PPI networks using STRING database
- **Enrichment Analysis**: KEGG, GO, and Reactome pathway enrichment
- **ADMET Prediction**: Drug-likeness evaluation (Lipinski, Veber rules)
- **Visualization**: Publication-ready figures

## Installation

```bash
# Clone the repository
git clone <repository-url>
cd herbal-network-pharmacology

# Create virtual environment
python -m venv venv
venv\Scripts\activate  # Windows
# source venv/bin/activate  # Linux/Mac

# Install dependencies
pip install -r requirements.txt
```

## Quick Start

### 1. Create Project Config

Copy and customize the config file:

```bash
cp config/default.yaml config/projects/my_project.yaml
```

Edit `my_project.yaml` with your plant and disease information.

### 2. Run Analysis

**Using CLI:**

```bash
# Show project info
python -m src.cli --config config/projects/mahkota_dewa_dn.yaml info

# Run individual steps
python -m src.cli -c config/projects/mahkota_dewa_dn.yaml collect
python -m src.cli -c config/projects/mahkota_dewa_dn.yaml predict-targets
python -m src.cli -c config/projects/mahkota_dewa_dn.yaml get-disease-genes
python -m src.cli -c config/projects/mahkota_dewa_dn.yaml build-network
python -m src.cli -c config/projects/mahkota_dewa_dn.yaml enrich
python -m src.cli -c config/projects/mahkota_dewa_dn.yaml admet
python -m src.cli -c config/projects/mahkota_dewa_dn.yaml visualize

# Or run all at once
python -m src.cli -c config/projects/mahkota_dewa_dn.yaml run-all
```

**Using Python:**

```python
from src.config_loader import load_config
from src.compounds.collector import CompoundCollector
from src.network.builder import NetworkBuilder

# Load config
config = load_config("config/projects/mahkota_dewa_dn.yaml")

# Collect compounds
collector = CompoundCollector(config)
collector.collect_from_literature()
collector.enrich_from_pubchem()
collector.save_compounds()

# Build network
builder = NetworkBuilder(config)
builder.find_common_targets()
network = builder.build_ppi_network()
```

## Project Structure

```
herbal-network-pharmacology/
├── config/
│   ├── default.yaml              # Default config template
│   └── projects/                 # Project-specific configs
├── src/
│   ├── cli.py                    # Command-line interface
│   ├── config_loader.py          # Config management
│   ├── compounds/                # Compound data collection
│   ├── targets/                  # Target prediction
│   ├── disease/                  # Disease gene retrieval
│   ├── network/                  # Network analysis
│   ├── enrichment/               # Pathway enrichment
│   ├── admet/                    # ADMET prediction
│   └── visualization/            # Figure generation
├── notebooks/
│   └── colab/                    # Google Colab notebooks
│       ├── 03_molecular_docking.ipynb
│       └── 04_md_simulation.ipynb
├── data/                         # Data files
├── outputs/                      # Generated outputs
└── requirements.txt
```

## Workflow Overview

```
1. Compound Collection
   └── KNApSAcK, PubChem, Literature
   
2. Target Prediction
   └── SwissTargetPrediction
   
3. Disease Genes
   └── GeneCards, DisGeNET, OMIM
   
4. Network Analysis
   └── STRING PPI, Hub Genes
   
5. Enrichment
   └── KEGG, GO, Reactome
   
6. ADMET Prediction
   └── Lipinski, Veber Rules
   
7. Molecular Docking (Colab)
   └── AutoDock Vina
   
8. MD Simulation (Colab)
   └── GROMACS 50ns
```

## Heavy Computation (Google Colab)

For molecular docking and MD simulation, use the Colab notebooks in `notebooks/colab/`:

1. Upload notebook to Google Colab
2. Upload required input files (compounds, protein structures)
3. Run the notebook
4. Download results

## Configuration

Key configuration options in YAML:

```yaml
plant:
  name: "Plant Name"
  latin_name: "Genus species"
  known_compounds:
    - name: "Compound1"
      pubchem_cid: 12345
      smiles: "CCO"

disease:
  name: "Disease Name"
  disgenet_id: "C0011849"

analysis:
  target_probability_threshold: 0.1
  hub_gene_top_n: 10
  enrichment_pvalue_threshold: 0.05
```

## Output Files

After running the workflow, you'll find:

```
data/<project>/
├── raw/
│   ├── compounds.csv
│   └── disease_genes.csv
├── processed/
│   ├── predicted_targets.csv
│   └── common_targets.txt
└── results/
    ├── network.graphml
    ├── hub_genes.csv
    ├── kegg_enrichment.csv
    └── admet_predictions.csv

outputs/figures/<project>/
├── compound_target_network.png
├── ppi_network.png
├── venn_diagram.png
├── kegg_barplot.png
└── enrichment_bubble.png
```

## License

MIT License

## Citation

If you use this toolkit in your research, please cite:

```
[Your citation here]
```
