"""
Plotter Module
==============

Creates publication-ready figures for network pharmacology analysis.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from pathlib import Path
from typing import List, Dict, Optional, Any, Tuple
import json

from ..config_loader import Config


class Plotter:
    """Creates publication-ready visualizations."""
    
    def __init__(self, config: Config):
        """
        Initialize the plotter.
        
        Args:
            config: Configuration object
        """
        self.config = config
        self._setup_style()
    
    def _setup_style(self):
        """Set up matplotlib style for publication-quality figures."""
        plt.style.use('seaborn-v0_8-whitegrid')
        
        # Get colors from config
        colors = self.config.get("output.colors", {})
        
        self.colors = {
            "primary": colors.get("primary", "#2E86AB"),
            "secondary": colors.get("secondary", "#A23B72"),
            "accent": colors.get("accent", "#F18F01"),
            "background": colors.get("background", "#FFFFFF"),
            "network_node": colors.get("network_node", "#4ECDC4"),
            "network_edge": colors.get("network_edge", "#95A5A6"),
        }
        
        # Publication settings
        plt.rcParams.update({
            'figure.dpi': self.config.get("output.figure_dpi", 300),
            'savefig.dpi': self.config.get("output.figure_dpi", 300),
            'font.size': 10,
            'axes.titlesize': 12,
            'axes.labelsize': 10,
            'xtick.labelsize': 9,
            'ytick.labelsize': 9,
            'legend.fontsize': 9,
            'figure.figsize': (8, 6),
        })
    
    def _save_figure(self, fig: plt.Figure, name: str) -> Path:
        """Save figure to output directory."""
        fmt = self.config.get("output.figure_format", "png")
        output_path = self.config.figures_dir / f"{name}.{fmt}"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        fig.savefig(output_path, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        
        return output_path
    
    def plot_compound_target_network(self) -> Optional[Path]:
        """
        Plot bipartite compound-target network.
        
        Returns:
            Path to saved figure
        """
        targets_path = self.config.data_dir / "processed" / "predicted_targets.csv"
        
        if not targets_path.exists():
            print("Predicted targets file not found.")
            return None
        
        df = pd.read_csv(targets_path)
        
        # Build network
        G = nx.Graph()
        
        compounds = df["compound_name"].dropna().unique()
        targets = df["gene_symbol"].dropna().unique()
        
        for compound in compounds:
            G.add_node(compound, node_type="compound")
        
        for target in targets:
            G.add_node(target, node_type="target")
        
        for _, row in df.iterrows():
            if pd.notna(row.get("compound_name")) and pd.notna(row.get("gene_symbol")):
                G.add_edge(row["compound_name"], row["gene_symbol"])
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Layout
        pos = nx.spring_layout(G, k=2, iterations=50, seed=42)
        
        # Draw compounds (circles)
        compound_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "compound"]
        nx.draw_networkx_nodes(G, pos, nodelist=compound_nodes, 
                               node_color=self.colors["primary"],
                               node_size=500, alpha=0.8, ax=ax)
        
        # Draw targets (squares using scatter)
        target_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "target"]
        target_pos = {n: pos[n] for n in target_nodes}
        for node, (x, y) in target_pos.items():
            ax.scatter(x, y, s=300, c=self.colors["secondary"], marker='s', alpha=0.7, zorder=2)
        
        # Draw edges
        nx.draw_networkx_edges(G, pos, alpha=0.3, 
                               edge_color=self.colors["network_edge"], ax=ax)
        
        # Labels
        nx.draw_networkx_labels(G, pos, font_size=7, ax=ax)
        
        ax.set_title(f"Compound-Target Network\n{self.config.plant_name}", fontsize=14, fontweight='bold')
        ax.axis('off')
        
        # Legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=self.colors["primary"], label=f'Compounds ({len(compound_nodes)})'),
            Patch(facecolor=self.colors["secondary"], label=f'Targets ({len(target_nodes)})'),
        ]
        ax.legend(handles=legend_elements, loc='upper left')
        
        return self._save_figure(fig, "compound_target_network")
    
    def plot_ppi_network(self) -> Optional[Path]:
        """
        Plot PPI network with hub genes highlighted.
        
        Returns:
            Path to saved figure
        """
        network_path = self.config.data_dir / "results" / "network.graphml"
        hub_path = self.config.data_dir / "results" / "hub_genes.csv"
        
        if not network_path.exists():
            print("Network file not found.")
            return None
        
        G = nx.read_graphml(network_path)
        
        # Load hub genes
        hub_genes = []
        if hub_path.exists():
            hub_df = pd.read_csv(hub_path)
            hub_genes = hub_df["gene"].tolist()[:10]  # Top 10
        
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Layout
        pos = nx.spring_layout(G, k=1.5, iterations=50, seed=42)
        
        # Node sizes based on degree
        degrees = dict(G.degree())
        node_sizes = [degrees.get(n, 1) * 50 + 100 for n in G.nodes()]
        
        # Node colors - highlight hub genes
        node_colors = []
        for node in G.nodes():
            if node in hub_genes:
                node_colors.append(self.colors["accent"])
            else:
                node_colors.append(self.colors["network_node"])
        
        # Draw
        nx.draw_networkx_edges(G, pos, alpha=0.2, 
                               edge_color=self.colors["network_edge"], ax=ax)
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, 
                               node_size=node_sizes, alpha=0.8, ax=ax)
        
        # Label hub genes only
        hub_labels = {n: n for n in hub_genes if n in G.nodes()}
        nx.draw_networkx_labels(G, pos, labels=hub_labels, font_size=8, 
                                font_weight='bold', ax=ax)
        
        ax.set_title(f"PPI Network - {self.config.disease_name}\n({G.number_of_nodes()} nodes, {G.number_of_edges()} edges)", 
                     fontsize=14, fontweight='bold')
        ax.axis('off')
        
        # Legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=self.colors["accent"], label='Hub Genes'),
            Patch(facecolor=self.colors["network_node"], label='Other Genes'),
        ]
        ax.legend(handles=legend_elements, loc='upper left')
        
        return self._save_figure(fig, "ppi_network")
    
    def plot_venn_diagram(self) -> Optional[Path]:
        """
        Plot Venn diagram of drug targets vs disease genes.
        
        Returns:
            Path to saved figure
        """
        try:
            from matplotlib_venn import venn2
        except ImportError:
            print("matplotlib-venn not installed. Install with: pip install matplotlib-venn")
            return None
        
        venn_path = self.config.data_dir / "processed" / "venn_data.json"
        
        if not venn_path.exists():
            print("Venn data not found.")
            return None
        
        with open(venn_path, "r") as f:
            venn_data = json.load(f)
        
        drug_only = len(venn_data.get("drug_targets_only", []))
        disease_only = len(venn_data.get("disease_genes_only", []))
        common = len(venn_data.get("common", []))
        
        fig, ax = plt.subplots(figsize=(8, 8))
        
        venn = venn2(subsets=(drug_only, disease_only, common),
                     set_labels=(f'Drug Targets\n({drug_only + common})', 
                                f'Disease Genes\n({disease_only + common})'),
                     ax=ax)
        
        # Customize colors
        if venn.get_patch_by_id('10'):
            venn.get_patch_by_id('10').set_color(self.colors["primary"])
            venn.get_patch_by_id('10').set_alpha(0.6)
        if venn.get_patch_by_id('01'):
            venn.get_patch_by_id('01').set_color(self.colors["secondary"])
            venn.get_patch_by_id('01').set_alpha(0.6)
        if venn.get_patch_by_id('11'):
            venn.get_patch_by_id('11').set_color(self.colors["accent"])
            venn.get_patch_by_id('11').set_alpha(0.8)
        
        ax.set_title(f"Target-Disease Gene Intersection\n{self.config.plant_name} vs {self.config.disease_name}",
                     fontsize=14, fontweight='bold')
        
        return self._save_figure(fig, "venn_diagram")
    
    def plot_enrichment_barplot(self) -> Optional[Path]:
        """
        Plot bar chart of top enriched pathways.
        
        Returns:
            Path to saved figure
        """
        kegg_path = self.config.data_dir / "results" / "kegg_enrichment.csv"
        
        if not kegg_path.exists():
            print("KEGG enrichment results not found.")
            return None
        
        df = pd.read_csv(kegg_path)
        
        # Get top pathways
        top_n = min(15, len(df))
        df_top = df.head(top_n).copy()
        
        # Clean term names
        if "Term" in df_top.columns:
            df_top["Term"] = df_top["Term"].apply(lambda x: x.split("(")[0].strip() if pd.notna(x) else x)
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Bar plot
        y_pos = range(len(df_top))
        
        # Use -log10(p-value) or Adjusted P-value
        pval_col = "Adjusted P-value" if "Adjusted P-value" in df_top.columns else "P-value"
        if pval_col in df_top.columns:
            values = -np.log10(df_top[pval_col].clip(lower=1e-50))
            xlabel = "-log10(Adjusted P-value)"
        else:
            values = df_top.get("Combined Score", range(len(df_top)))
            xlabel = "Combined Score"
        
        bars = ax.barh(y_pos, values, color=self.colors["primary"], alpha=0.8)
        
        ax.set_yticks(y_pos)
        ax.set_yticklabels(df_top.get("Term", df_top.index))
        ax.invert_yaxis()
        
        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_title(f"Top KEGG Pathways\n{self.config.plant_name} - {self.config.disease_name}",
                     fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        
        return self._save_figure(fig, "kegg_barplot")
    
    def plot_enrichment_bubble(self) -> Optional[Path]:
        """
        Plot bubble chart for enrichment results.
        
        Returns:
            Path to saved figure
        """
        combined_path = self.config.data_dir / "results" / "top_pathways_combined.csv"
        
        if not combined_path.exists():
            print("Combined pathways file not found.")
            return None
        
        df = pd.read_csv(combined_path)
        
        if len(df) == 0:
            return None
        
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Prepare data
        if "Adjusted P-value" in df.columns:
            df["neg_log_pval"] = -np.log10(df["Adjusted P-value"].clip(lower=1e-50))
        else:
            df["neg_log_pval"] = 1
        
        if "Overlap" in df.columns:
            # Parse overlap (e.g., "5/100")
            try:
                df["gene_count"] = df["Overlap"].apply(lambda x: int(str(x).split("/")[0]) if pd.notna(x) else 1)
            except:
                df["gene_count"] = 1
        else:
            df["gene_count"] = 10
        
        # Color by source
        sources = df["source"].unique() if "source" in df.columns else ["KEGG"]
        color_map = {s: plt.cm.Set2(i) for i, s in enumerate(sources)}
        df["color"] = df.get("source", "KEGG").map(color_map)
        
        # Plot
        scatter = ax.scatter(df["neg_log_pval"], range(len(df)),
                            s=df["gene_count"] * 20,
                            c=df["color"],
                            alpha=0.7)
        
        ax.set_yticks(range(len(df)))
        term_col = "Term" if "Term" in df.columns else df.columns[0]
        labels = df[term_col].apply(lambda x: str(x)[:50] + "..." if len(str(x)) > 50 else str(x))
        ax.set_yticklabels(labels)
        ax.invert_yaxis()
        
        ax.set_xlabel("-log10(Adjusted P-value)", fontsize=11)
        ax.set_title(f"Pathway Enrichment Overview\n{self.config.plant_name} - {self.config.disease_name}",
                     fontsize=14, fontweight='bold')
        
        # Legend for sources
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor=color_map[s], label=s) for s in sources]
        ax.legend(handles=legend_elements, loc='lower right', title='Source')
        
        plt.tight_layout()
        
        return self._save_figure(fig, "enrichment_bubble")
    
    def plot_admet_summary(self) -> Optional[Path]:
        """
        Plot ADMET summary.
        
        Returns:
            Path to saved figure
        """
        admet_path = self.config.data_dir / "results" / "admet_predictions.csv"
        
        if not admet_path.exists():
            print("ADMET predictions not found.")
            return None
        
        df = pd.read_csv(admet_path)
        
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # Pie chart for drug-likeness
        if "drug_like" in df.columns:
            drug_like = df["drug_like"].sum()
            not_drug_like = len(df) - drug_like
            
            axes[0].pie([drug_like, not_drug_like],
                       labels=[f'Drug-like ({drug_like})', f'Non drug-like ({not_drug_like})'],
                       colors=[self.colors["primary"], self.colors["secondary"]],
                       autopct='%1.1f%%',
                       startangle=90)
            axes[0].set_title("Drug-likeness Assessment", fontsize=12, fontweight='bold')
        
        # Violations distribution
        if "lipinski_violations" in df.columns:
            viol_counts = df["lipinski_violations"].value_counts().sort_index()
            axes[1].bar(viol_counts.index, viol_counts.values, 
                       color=self.colors["primary"], alpha=0.8)
            axes[1].set_xlabel("Number of Lipinski Violations")
            axes[1].set_ylabel("Number of Compounds")
            axes[1].set_title("Lipinski Rule Violations", fontsize=12, fontweight='bold')
            axes[1].set_xticks(range(5))
        
        plt.suptitle(f"ADMET Summary - {self.config.plant_name}", fontsize=14, fontweight='bold', y=1.02)
        plt.tight_layout()
        
        return self._save_figure(fig, "admet_summary")
    
    def plot_all(self) -> List[Path]:
        """
        Generate all visualizations.
        
        Returns:
            List of saved figure paths
        """
        paths = []
        
        # Compound-target network
        path = self.plot_compound_target_network()
        if path:
            paths.append(path)
        
        # PPI network
        path = self.plot_ppi_network()
        if path:
            paths.append(path)
        
        # Venn diagram
        path = self.plot_venn_diagram()
        if path:
            paths.append(path)
        
        # Enrichment
        path = self.plot_enrichment_barplot()
        if path:
            paths.append(path)
        
        path = self.plot_enrichment_bubble()
        if path:
            paths.append(path)
        
        # ADMET
        path = self.plot_admet_summary()
        if path:
            paths.append(path)
        
        return paths
