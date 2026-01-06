"""
Network Builder Module
======================

Builds protein-protein interaction networks using STRING database.
"""

import pandas as pd
import requests
import networkx as nx
from pathlib import Path
from typing import List, Dict, Optional, Set, Any
from tqdm import tqdm
import json

from ..config_loader import Config


class NetworkBuilder:
    """Builds PPI networks from drug targets and disease genes."""
    
    STRING_API_URL = "https://string-db.org/api"
    
    def __init__(self, config: Config):
        """
        Initialize the builder.
        
        Args:
            config: Configuration object
        """
        self.config = config
        self.drug_targets: List[str] = []
        self.disease_genes: List[str] = []
        self.common_targets: List[str] = []
        self.network: Optional[nx.Graph] = None
        self._species = config.get("api.string.species", 9606)  # Human
        self._score_threshold = config.get("api.string.score_threshold", 400)
    
    def load_drug_targets(self, filename: str = "unique_targets.txt") -> List[str]:
        """
        Load predicted drug targets.
        
        Args:
            filename: Filename with target list
            
        Returns:
            List of target gene symbols
        """
        path = self.config.data_dir / "processed" / filename
        
        if path.exists():
            with open(path, "r") as f:
                self.drug_targets = [line.strip().upper() for line in f if line.strip()]
        else:
            # Try loading from CSV
            csv_path = self.config.data_dir / "processed" / "predicted_targets.csv"
            if csv_path.exists():
                df = pd.read_csv(csv_path)
                self.drug_targets = list(df["gene_symbol"].dropna().str.upper().unique())
        
        return self.drug_targets
    
    def load_disease_genes(self, filename: str = "disease_gene_symbols.txt") -> List[str]:
        """
        Load disease-associated genes.
        
        Args:
            filename: Filename with gene list
            
        Returns:
            List of gene symbols
        """
        path = self.config.data_dir / "raw" / filename
        
        if path.exists():
            with open(path, "r") as f:
                self.disease_genes = [line.strip().upper() for line in f if line.strip()]
        else:
            # Try loading from CSV
            csv_path = self.config.data_dir / "raw" / "disease_genes.csv"
            if csv_path.exists():
                df = pd.read_csv(csv_path)
                self.disease_genes = list(df["gene_symbol"].dropna().str.upper().unique())
        
        return self.disease_genes
    
    def find_common_targets(self) -> List[str]:
        """
        Find intersection between drug targets and disease genes.
        
        Returns:
            List of common gene symbols
        """
        if not self.drug_targets:
            self.load_drug_targets()
        if not self.disease_genes:
            self.load_disease_genes()
        
        drug_set = set(self.drug_targets)
        disease_set = set(self.disease_genes)
        
        self.common_targets = list(drug_set & disease_set)
        
        print(f"Drug targets: {len(self.drug_targets)}")
        print(f"Disease genes: {len(self.disease_genes)}")
        print(f"Common targets: {len(self.common_targets)}")
        
        # Save Venn diagram data
        venn_data = {
            "drug_targets_only": list(drug_set - disease_set),
            "disease_genes_only": list(disease_set - drug_set),
            "common": self.common_targets
        }
        
        venn_path = self.config.data_dir / "processed" / "venn_data.json"
        with open(venn_path, "w") as f:
            json.dump(venn_data, f, indent=2)
        
        # Save common targets
        common_path = self.config.data_dir / "processed" / "common_targets.txt"
        with open(common_path, "w") as f:
            f.write("\n".join(self.common_targets))
        
        return self.common_targets
    
    def _query_string_network(self, genes: List[str]) -> List[Dict]:
        """
        Query STRING database for PPI network.
        
        Args:
            genes: List of gene symbols
            
        Returns:
            List of interaction records
        """
        url = f"{self.STRING_API_URL}/json/network"
        
        params = {
            "identifiers": "%0d".join(genes),
            "species": self._species,
            "required_score": self._score_threshold,
            "network_type": "functional"
        }
        
        try:
            response = requests.post(url, data=params, timeout=120)
            
            if response.status_code == 200:
                return response.json()
            else:
                print(f"STRING API error: {response.status_code}")
                return []
                
        except Exception as e:
            print(f"STRING query failed: {e}")
            return []
    
    def _query_string_enrichment(self, genes: List[str]) -> Dict:
        """
        Get functional enrichment from STRING.
        
        Args:
            genes: List of gene symbols
            
        Returns:
            Enrichment results
        """
        url = f"{self.STRING_API_URL}/json/enrichment"
        
        params = {
            "identifiers": "%0d".join(genes),
            "species": self._species,
        }
        
        try:
            response = requests.post(url, data=params, timeout=120)
            
            if response.status_code == 200:
                return response.json()
            else:
                return {}
                
        except Exception:
            return {}
    
    def build_ppi_network(self, genes: Optional[List[str]] = None) -> nx.Graph:
        """
        Build PPI network from STRING database.
        
        Args:
            genes: Optional gene list. If None, uses common_targets.
            
        Returns:
            NetworkX graph
        """
        if genes is None:
            if not self.common_targets:
                self.find_common_targets()
            genes = self.common_targets
        
        if not genes:
            print("No genes to build network with.")
            return nx.Graph()
        
        print(f"Querying STRING database for {len(genes)} genes...")
        
        interactions = self._query_string_network(genes)
        
        # Build NetworkX graph
        self.network = nx.Graph()
        
        for interaction in interactions:
            node1 = interaction.get("preferredName_A", interaction.get("stringId_A", ""))
            node2 = interaction.get("preferredName_B", interaction.get("stringId_B", ""))
            score = interaction.get("score", 0)
            
            if node1 and node2:
                self.network.add_edge(node1, node2, weight=score)
        
        # Add isolated nodes (genes with no interactions)
        for gene in genes:
            if gene not in self.network:
                self.network.add_node(gene)
        
        print(f"Network built: {self.network.number_of_nodes()} nodes, {self.network.number_of_edges()} edges")
        
        return self.network
    
    def build_compound_target_network(self) -> nx.Graph:
        """
        Build bipartite network of compounds and their targets.
        
        Returns:
            NetworkX graph
        """
        targets_path = self.config.data_dir / "processed" / "predicted_targets.csv"
        
        if not targets_path.exists():
            print("Predicted targets file not found.")
            return nx.Graph()
        
        df = pd.read_csv(targets_path)
        
        G = nx.Graph()
        
        for _, row in df.iterrows():
            compound = row.get("compound_name", "")
            target = row.get("gene_symbol", "")
            
            if compound and target:
                # Add nodes with type attribute
                G.add_node(compound, node_type="compound")
                G.add_node(target, node_type="target")
                G.add_edge(compound, target)
        
        return G
    
    def export_to_cytoscape(self, filename: str = "network.graphml") -> Path:
        """
        Export network in GraphML format for Cytoscape.
        
        Args:
            filename: Output filename
            
        Returns:
            Path to saved file
        """
        if self.network is None:
            print("No network to export. Build network first.")
            return Path()
        
        output_path = self.config.data_dir / "results" / filename
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        nx.write_graphml(self.network, output_path)
        
        return output_path
    
    def save_network(self) -> Path:
        """
        Save network data and export for Cytoscape.
        
        Returns:
            Path to saved file
        """
        if self.network is None:
            print("No network to save.")
            return Path()
        
        # Save as GraphML
        graphml_path = self.export_to_cytoscape()
        
        # Save edge list as CSV
        edges_data = []
        for u, v, d in self.network.edges(data=True):
            edges_data.append({
                "source": u,
                "target": v,
                "weight": d.get("weight", 1.0)
            })
        
        edges_df = pd.DataFrame(edges_data)
        edges_path = self.config.data_dir / "results" / "network_edges.csv"
        edges_df.to_csv(edges_path, index=False)
        
        # Save node list
        nodes_data = []
        for node in self.network.nodes():
            nodes_data.append({
                "node": node,
                "degree": self.network.degree(node)
            })
        
        nodes_df = pd.DataFrame(nodes_data)
        nodes_path = self.config.data_dir / "results" / "network_nodes.csv"
        nodes_df.to_csv(nodes_path, index=False)
        
        return graphml_path
