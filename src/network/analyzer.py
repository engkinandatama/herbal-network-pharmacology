"""
Network Analyzer Module
=======================

Analyzes network topology and identifies hub genes.
"""

import pandas as pd
import networkx as nx
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Any

from ..config_loader import Config


class NetworkAnalyzer:
    """Analyzes PPI network topology and identifies key nodes."""
    
    def __init__(self, config: Config, network: Optional[nx.Graph] = None):
        """
        Initialize the analyzer.
        
        Args:
            config: Configuration object
            network: NetworkX graph to analyze
        """
        self.config = config
        self.network = network
        self.hub_genes: List[Dict[str, Any]] = []
        self.centrality_metrics: Dict[str, Dict[str, float]] = {}
    
    def load_network(self, filename: str = "network.graphml") -> nx.Graph:
        """
        Load network from GraphML file.
        
        Args:
            filename: GraphML filename
            
        Returns:
            NetworkX graph
        """
        path = self.config.data_dir / "results" / filename
        
        if path.exists():
            self.network = nx.read_graphml(path)
            return self.network
        else:
            raise FileNotFoundError(f"Network file not found: {path}")
    
    def calculate_degree_centrality(self) -> Dict[str, float]:
        """
        Calculate degree centrality for all nodes.
        
        Returns:
            Dictionary of node -> degree centrality
        """
        if self.network is None:
            raise ValueError("No network loaded.")
        
        dc = nx.degree_centrality(self.network)
        self.centrality_metrics["degree"] = dc
        return dc
    
    def calculate_betweenness_centrality(self) -> Dict[str, float]:
        """
        Calculate betweenness centrality for all nodes.
        
        Returns:
            Dictionary of node -> betweenness centrality
        """
        if self.network is None:
            raise ValueError("No network loaded.")
        
        bc = nx.betweenness_centrality(self.network)
        self.centrality_metrics["betweenness"] = bc
        return bc
    
    def calculate_closeness_centrality(self) -> Dict[str, float]:
        """
        Calculate closeness centrality for all nodes.
        
        Returns:
            Dictionary of node -> closeness centrality
        """
        if self.network is None:
            raise ValueError("No network loaded.")
        
        cc = nx.closeness_centrality(self.network)
        self.centrality_metrics["closeness"] = cc
        return cc
    
    def calculate_eigenvector_centrality(self) -> Dict[str, float]:
        """
        Calculate eigenvector centrality for all nodes.
        
        Returns:
            Dictionary of node -> eigenvector centrality
        """
        if self.network is None:
            raise ValueError("No network loaded.")
        
        try:
            ec = nx.eigenvector_centrality(self.network, max_iter=1000)
            self.centrality_metrics["eigenvector"] = ec
            return ec
        except nx.PowerIterationFailedConvergence:
            print("Warning: Eigenvector centrality did not converge.")
            return {}
    
    def calculate_all_centralities(self) -> pd.DataFrame:
        """
        Calculate all centrality metrics.
        
        Returns:
            DataFrame with all centrality values
        """
        if self.network is None:
            raise ValueError("No network loaded.")
        
        self.calculate_degree_centrality()
        self.calculate_betweenness_centrality()
        self.calculate_closeness_centrality()
        self.calculate_eigenvector_centrality()
        
        # Combine into DataFrame
        data = []
        for node in self.network.nodes():
            node_data = {
                "gene": node,
                "degree": self.network.degree(node),
                "degree_centrality": self.centrality_metrics.get("degree", {}).get(node, 0),
                "betweenness_centrality": self.centrality_metrics.get("betweenness", {}).get(node, 0),
                "closeness_centrality": self.centrality_metrics.get("closeness", {}).get(node, 0),
                "eigenvector_centrality": self.centrality_metrics.get("eigenvector", {}).get(node, 0),
            }
            data.append(node_data)
        
        return pd.DataFrame(data)
    
    def identify_hub_genes(self, top_n: Optional[int] = None, method: str = "degree") -> List[Dict[str, Any]]:
        """
        Identify hub genes based on centrality metrics.
        
        Args:
            top_n: Number of top genes to return. If None, uses config value.
            method: Centrality method to use (degree, betweenness, closeness, eigenvector, combined)
            
        Returns:
            List of hub gene dictionaries
        """
        if self.network is None:
            raise ValueError("No network loaded.")
        
        if top_n is None:
            top_n = self.config.get("analysis.hub_gene_top_n", 10)
        
        # Calculate centralities if not done
        centrality_df = self.calculate_all_centralities()
        
        # Sort by method
        if method == "combined":
            # Normalize and combine metrics
            for col in ["degree_centrality", "betweenness_centrality", "closeness_centrality", "eigenvector_centrality"]:
                if col in centrality_df.columns:
                    max_val = centrality_df[col].max()
                    if max_val > 0:
                        centrality_df[f"{col}_norm"] = centrality_df[col] / max_val
                    else:
                        centrality_df[f"{col}_norm"] = 0
            
            # Combined score (average of normalized metrics)
            norm_cols = [c for c in centrality_df.columns if c.endswith("_norm")]
            centrality_df["combined_score"] = centrality_df[norm_cols].mean(axis=1)
            
            sorted_df = centrality_df.sort_values("combined_score", ascending=False)
        else:
            sort_col = f"{method}_centrality" if method != "degree" else "degree"
            sorted_df = centrality_df.sort_values(sort_col, ascending=False)
        
        # Get top N
        hub_df = sorted_df.head(top_n)
        self.hub_genes = hub_df.to_dict("records")
        
        return self.hub_genes
    
    def get_network_statistics(self) -> Dict[str, Any]:
        """
        Calculate network statistics.
        
        Returns:
            Dictionary of network metrics
        """
        if self.network is None:
            raise ValueError("No network loaded.")
        
        stats = {
            "nodes": self.network.number_of_nodes(),
            "edges": self.network.number_of_edges(),
            "density": nx.density(self.network),
            "average_degree": sum(dict(self.network.degree()).values()) / self.network.number_of_nodes() if self.network.number_of_nodes() > 0 else 0,
        }
        
        # Check if connected
        if nx.is_connected(self.network):
            stats["diameter"] = nx.diameter(self.network)
            stats["average_path_length"] = nx.average_shortest_path_length(self.network)
        else:
            # Get largest component stats
            largest_cc = max(nx.connected_components(self.network), key=len)
            subgraph = self.network.subgraph(largest_cc)
            stats["largest_component_nodes"] = len(largest_cc)
            stats["number_of_components"] = nx.number_connected_components(self.network)
            if len(largest_cc) > 1:
                stats["diameter_largest"] = nx.diameter(subgraph)
                stats["average_path_length_largest"] = nx.average_shortest_path_length(subgraph)
        
        # Clustering coefficient
        stats["average_clustering"] = nx.average_clustering(self.network)
        
        return stats
    
    def perform_clustering(self) -> Dict[str, List[str]]:
        """
        Perform community detection using Louvain algorithm.
        
        Returns:
            Dictionary of cluster_id -> list of genes
        """
        if self.network is None:
            raise ValueError("No network loaded.")
        
        try:
            from networkx.algorithms import community
            communities = community.louvain_communities(self.network)
            
            clusters = {}
            for i, comm in enumerate(communities):
                clusters[f"cluster_{i+1}"] = list(comm)
            
            return clusters
            
        except Exception as e:
            print(f"Clustering failed: {e}")
            return {}
    
    def save_analysis(self) -> Path:
        """
        Save analysis results to files.
        
        Returns:
            Path to main output file
        """
        output_dir = self.config.data_dir / "results"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save centrality metrics
        centrality_df = self.calculate_all_centralities()
        centrality_path = output_dir / "network_centralities.csv"
        centrality_df.to_csv(centrality_path, index=False)
        
        # Save hub genes
        if self.hub_genes:
            hub_df = pd.DataFrame(self.hub_genes)
            hub_path = output_dir / "hub_genes.csv"
            hub_df.to_csv(hub_path, index=False)
        
        # Save network statistics
        stats = self.get_network_statistics()
        stats_path = output_dir / "network_statistics.json"
        import json
        with open(stats_path, "w") as f:
            json.dump(stats, f, indent=2)
        
        # Save clusters
        clusters = self.perform_clustering()
        if clusters:
            clusters_path = output_dir / "network_clusters.json"
            with open(clusters_path, "w") as f:
                json.dump(clusters, f, indent=2)
        
        return centrality_path
