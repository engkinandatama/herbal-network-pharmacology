"""
Enrichment Analysis Module
==========================

Performs pathway enrichment analysis using:
- KEGG pathways
- GO terms (BP, MF, CC)
- Reactome (optional)
"""

import pandas as pd
from pathlib import Path
from typing import List, Dict, Optional, Any
import json

try:
    import gseapy as gp
except ImportError:
    gp = None

from ..config_loader import Config


class EnrichmentAnalyzer:
    """Performs pathway and GO enrichment analysis."""
    
    def __init__(self, config: Config):
        """
        Initialize the analyzer.
        
        Args:
            config: Configuration object
        """
        self.config = config
        self.gene_list: List[str] = []
        self.kegg_results: Optional[pd.DataFrame] = None
        self.go_bp_results: Optional[pd.DataFrame] = None
        self.go_mf_results: Optional[pd.DataFrame] = None
        self.go_cc_results: Optional[pd.DataFrame] = None
        self.reactome_results: Optional[pd.DataFrame] = None
        
        self._pvalue_threshold = config.get("analysis.enrichment_pvalue_threshold", 0.05)
        self._top_n = config.get("analysis.enrichment_top_n", 20)
    
    def load_genes(self, filename: str = "common_targets.txt") -> List[str]:
        """
        Load gene list for enrichment.
        
        Args:
            filename: File with gene list
            
        Returns:
            List of gene symbols
        """
        path = self.config.data_dir / "processed" / filename
        
        if path.exists():
            with open(path, "r") as f:
                self.gene_list = [line.strip().upper() for line in f if line.strip()]
        else:
            # Try hub genes
            hub_path = self.config.data_dir / "results" / "hub_genes.csv"
            if hub_path.exists():
                df = pd.read_csv(hub_path)
                self.gene_list = df["gene"].tolist()
        
        return self.gene_list
    
    def kegg_enrichment(self, gene_list: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Perform KEGG pathway enrichment.
        
        Args:
            gene_list: Optional gene list. If None, uses loaded genes.
            
        Returns:
            DataFrame with enrichment results
        """
        if gp is None:
            print("gseapy not installed. Install with: pip install gseapy")
            return pd.DataFrame()
        
        if gene_list is None:
            if not self.gene_list:
                self.load_genes()
            gene_list = self.gene_list
        
        if not gene_list:
            print("No genes for enrichment analysis.")
            return pd.DataFrame()
        
        try:
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=["KEGG_2021_Human"],
                organism="human",
                outdir=None,
                cutoff=self._pvalue_threshold
            )
            
            self.kegg_results = enr.results
            return self.kegg_results
            
        except Exception as e:
            print(f"KEGG enrichment failed: {e}")
            return pd.DataFrame()
    
    def go_enrichment(self, gene_list: Optional[List[str]] = None) -> Dict[str, pd.DataFrame]:
        """
        Perform GO enrichment (BP, MF, CC).
        
        Args:
            gene_list: Optional gene list. If None, uses loaded genes.
            
        Returns:
            Dictionary with BP, MF, CC DataFrames
        """
        if gp is None:
            print("gseapy not installed. Install with: pip install gseapy")
            return {}
        
        if gene_list is None:
            if not self.gene_list:
                self.load_genes()
            gene_list = self.gene_list
        
        if not gene_list:
            print("No genes for enrichment analysis.")
            return {}
        
        results = {}
        
        # GO Biological Process
        try:
            enr_bp = gp.enrichr(
                gene_list=gene_list,
                gene_sets=["GO_Biological_Process_2021"],
                organism="human",
                outdir=None,
                cutoff=self._pvalue_threshold
            )
            self.go_bp_results = enr_bp.results
            results["BP"] = self.go_bp_results
        except Exception as e:
            print(f"GO BP enrichment failed: {e}")
        
        # GO Molecular Function
        try:
            enr_mf = gp.enrichr(
                gene_list=gene_list,
                gene_sets=["GO_Molecular_Function_2021"],
                organism="human",
                outdir=None,
                cutoff=self._pvalue_threshold
            )
            self.go_mf_results = enr_mf.results
            results["MF"] = self.go_mf_results
        except Exception as e:
            print(f"GO MF enrichment failed: {e}")
        
        # GO Cellular Component
        try:
            enr_cc = gp.enrichr(
                gene_list=gene_list,
                gene_sets=["GO_Cellular_Component_2021"],
                organism="human",
                outdir=None,
                cutoff=self._pvalue_threshold
            )
            self.go_cc_results = enr_cc.results
            results["CC"] = self.go_cc_results
        except Exception as e:
            print(f"GO CC enrichment failed: {e}")
        
        return results
    
    def reactome_enrichment(self, gene_list: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Perform Reactome pathway enrichment.
        
        Args:
            gene_list: Optional gene list. If None, uses loaded genes.
            
        Returns:
            DataFrame with enrichment results
        """
        if gp is None:
            print("gseapy not installed.")
            return pd.DataFrame()
        
        if gene_list is None:
            if not self.gene_list:
                self.load_genes()
            gene_list = self.gene_list
        
        if not gene_list:
            return pd.DataFrame()
        
        try:
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=["Reactome_2022"],
                organism="human",
                outdir=None,
                cutoff=self._pvalue_threshold
            )
            
            self.reactome_results = enr.results
            return self.reactome_results
            
        except Exception as e:
            print(f"Reactome enrichment failed: {e}")
            return pd.DataFrame()
    
    def run_all_enrichments(self, gene_list: Optional[List[str]] = None) -> Dict[str, pd.DataFrame]:
        """
        Run all enrichment analyses.
        
        Args:
            gene_list: Optional gene list
            
        Returns:
            Dictionary with all results
        """
        results = {}
        
        # KEGG
        kegg = self.kegg_enrichment(gene_list)
        if not kegg.empty:
            results["KEGG"] = kegg
        
        # GO
        go = self.go_enrichment(gene_list)
        results.update({f"GO_{k}": v for k, v in go.items() if not v.empty})
        
        # Reactome
        reactome = self.reactome_enrichment(gene_list)
        if not reactome.empty:
            results["Reactome"] = reactome
        
        return results
    
    def get_top_pathways(self, n: Optional[int] = None) -> pd.DataFrame:
        """
        Get top enriched pathways across all analyses.
        
        Args:
            n: Number of top pathways per analysis
            
        Returns:
            Combined DataFrame
        """
        if n is None:
            n = self._top_n
        
        dfs = []
        
        if self.kegg_results is not None and not self.kegg_results.empty:
            top_kegg = self.kegg_results.head(n).copy()
            top_kegg["source"] = "KEGG"
            dfs.append(top_kegg)
        
        if self.go_bp_results is not None and not self.go_bp_results.empty:
            top_bp = self.go_bp_results.head(n).copy()
            top_bp["source"] = "GO_BP"
            dfs.append(top_bp)
        
        if self.go_mf_results is not None and not self.go_mf_results.empty:
            top_mf = self.go_mf_results.head(n).copy()
            top_mf["source"] = "GO_MF"
            dfs.append(top_mf)
        
        if self.go_cc_results is not None and not self.go_cc_results.empty:
            top_cc = self.go_cc_results.head(n).copy()
            top_cc["source"] = "GO_CC"
            dfs.append(top_cc)
        
        if self.reactome_results is not None and not self.reactome_results.empty:
            top_react = self.reactome_results.head(n).copy()
            top_react["source"] = "Reactome"
            dfs.append(top_react)
        
        if dfs:
            return pd.concat(dfs, ignore_index=True)
        
        return pd.DataFrame()
    
    def save_results(self) -> Path:
        """
        Save all enrichment results.
        
        Returns:
            Path to output directory
        """
        output_dir = self.config.data_dir / "results"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save individual results
        if self.kegg_results is not None and not self.kegg_results.empty:
            self.kegg_results.to_csv(output_dir / "kegg_enrichment.csv", index=False)
        
        if self.go_bp_results is not None and not self.go_bp_results.empty:
            self.go_bp_results.to_csv(output_dir / "go_bp_enrichment.csv", index=False)
        
        if self.go_mf_results is not None and not self.go_mf_results.empty:
            self.go_mf_results.to_csv(output_dir / "go_mf_enrichment.csv", index=False)
        
        if self.go_cc_results is not None and not self.go_cc_results.empty:
            self.go_cc_results.to_csv(output_dir / "go_cc_enrichment.csv", index=False)
        
        if self.reactome_results is not None and not self.reactome_results.empty:
            self.reactome_results.to_csv(output_dir / "reactome_enrichment.csv", index=False)
        
        # Save combined top pathways
        top_pathways = self.get_top_pathways()
        if not top_pathways.empty:
            top_pathways.to_csv(output_dir / "top_pathways_combined.csv", index=False)
        
        return output_dir
