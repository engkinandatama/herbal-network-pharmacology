"""
Disease Gene Collector Module
=============================

Retrieves disease-associated genes from:
- GeneCards
- DisGeNET
- OMIM
"""

import os
import pandas as pd
import requests
import time
from pathlib import Path
from typing import List, Dict, Optional, Set, Any
from tqdm import tqdm

# Load environment variables from .env file
try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass  # python-dotenv not installed

from ..config_loader import Config


class DiseaseGeneCollector:
    """Collects disease-associated genes from multiple databases."""
    
    # OpenTargets GraphQL API - FREE, no auth required
    OPENTARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"
    
    def __init__(self, config: Config):
        """
        Initialize the collector.
        
        Args:
            config: Configuration object
        """
        self.config = config
        self.genes_opentargets: List[Dict[str, Any]] = []
        self.genes_genecards: List[Dict[str, Any]] = []
        self.genes_disgenet: List[Dict[str, Any]] = []
        self.genes_omim: List[Dict[str, Any]] = []
        self.all_genes: List[Dict[str, Any]] = []
        self._request_delay = 1.0
    
    def fetch_opentargets(self, score_threshold: float = 0.1) -> List[Dict[str, Any]]:
        """
        Fetch genes from OpenTargets Platform (FREE, no auth required).
        
        Args:
            score_threshold: Minimum association score (0-1)
            
        Returns:
            List of gene dictionaries
        """
        efo_id = self.config.get("disease.opentargets_id", "")
        disease_name = self.config.disease_name
        
        if not efo_id:
            # Try to construct from common patterns
            print(f"OpenTargets EFO ID not specified for {disease_name}")
            return self.genes_opentargets
        
        # GraphQL query for disease-gene associations
        query = """
        query diseaseAssociations($diseaseId: String!) {
          disease(efoId: $diseaseId) {
            id
            name
            associatedTargets(page: {index: 0, size: 1000}) {
              count
              rows {
                target {
                  id
                  approvedSymbol
                  approvedName
                  proteinIds {
                    id
                    source
                  }
                }
                score
              }
            }
          }
        }
        """
        
        try:
            response = requests.post(
                self.OPENTARGETS_URL,
                json={"query": query, "variables": {"diseaseId": efo_id}},
                timeout=120
            )
            
            if response.status_code == 200:
                data = response.json()
                
                if 'data' in data and data['data'].get('disease'):
                    disease = data['data']['disease']
                    targets = disease['associatedTargets']['rows']
                    
                    for row in targets:
                        if row['score'] >= score_threshold:
                            target = row['target']
                            
                            # Get UniProt ID if available
                            uniprot_id = ""
                            for pid in target.get('proteinIds', []):
                                if pid.get('source') == 'uniprot_swissprot':
                                    uniprot_id = pid.get('id', '')
                                    break
                            
                            gene_data = {
                                "gene_symbol": target.get('approvedSymbol', ''),
                                "gene_name": target.get('approvedName', ''),
                                "ensembl_id": target.get('id', ''),
                                "uniprot_id": uniprot_id,
                                "score": row['score'],
                                "disease_name": disease_name,
                                "source": "OpenTargets"
                            }
                            self.genes_opentargets.append(gene_data)
                    
                    print(f"OpenTargets: Found {len(self.genes_opentargets)} genes (score >= {score_threshold})")
                else:
                    print("Disease not found in OpenTargets")
                    
            else:
                print(f"OpenTargets API error: {response.status_code}")
                
        except Exception as e:
            print(f"OpenTargets fetch failed: {e}")
        
        return self.genes_opentargets
    
    def fetch_genecards(self) -> List[Dict[str, Any]]:
        """
        Fetch genes from GeneCards.
        
        Note: GeneCards requires API access for bulk queries.
        This method provides a template for manual search.
        
        Returns:
            List of gene dictionaries
        """
        disease_name = self.config.disease_name
        keywords = self.config.get("disease.search_keywords", [disease_name])
        
        print(f"Searching GeneCards for: {disease_name}")
        print("Note: GeneCards requires manual search or API subscription.")
        print(f"Visit: https://www.genecards.org/Search/Keyword?queryString={disease_name.replace(' ', '+')}")
        
        # For now, return empty - user needs to manually download
        # or use the GeneCards API if they have access
        
        # Check if manual data exists
        manual_path = self.config.data_dir / "raw" / "genecards_genes.csv"
        if manual_path.exists():
            df = pd.read_csv(manual_path)
            self.genes_genecards = df.to_dict("records")
            return self.genes_genecards
        
        return self.genes_genecards
    
    def fetch_disgenet(self, api_key: Optional[str] = None) -> List[Dict[str, Any]]:
        """
        Fetch genes from DisGeNET.
        
        Args:
            api_key: Optional DisGeNET API key. If not provided, 
                     loads from DISGENET_API_KEY environment variable.
            
        Returns:
            List of gene dictionaries
        """
        disease_id = self.config.get("disease.disgenet_id", "")
        disease_name = self.config.disease_name
        
        if not disease_id:
            print(f"Warning: DisGeNET ID not specified for {disease_name}")
            print("You can find the ID at: https://www.disgenet.org/search")
            return self.genes_disgenet
        
        # Load API key from environment if not provided
        if not api_key:
            api_key = os.getenv("DISGENET_API_KEY")
        
        if not api_key:
            print("Warning: No DisGeNET API key found.")
            print("Set DISGENET_API_KEY in .env file or pass directly.")
            return self.genes_disgenet
        
        # DisGeNET API endpoint
        base_url = "https://www.disgenet.org/api"
        
        headers = {
            "Accept": "application/json",
            "Authorization": f"Bearer {api_key}"
        }
        
        try:
            # Get gene-disease associations
            url = f"{base_url}/gda/disease/{disease_id}"
            
            response = requests.get(url, headers=headers, timeout=60)
            
            if response.status_code == 200:
                data = response.json()
                
                for item in data:
                    gene_data = {
                        "gene_symbol": item.get("gene_symbol", ""),
                        "gene_id": item.get("geneid", ""),
                        "uniprot_id": item.get("uniprotid", ""),
                        "score": item.get("score", 0),
                        "disease_name": disease_name,
                        "source": "DisGeNET"
                    }
                    self.genes_disgenet.append(gene_data)
            
            elif response.status_code == 401:
                print("DisGeNET API requires authentication.")
                print("Get free API key at: https://www.disgenet.org/api")
                print("Or manually download from: https://www.disgenet.org/search")
            
            else:
                print(f"DisGeNET API error: {response.status_code}")
                
        except Exception as e:
            print(f"DisGeNET fetch failed: {e}")
        
        # Check for manual data
        manual_path = self.config.data_dir / "raw" / "disgenet_genes.csv"
        if manual_path.exists() and not self.genes_disgenet:
            df = pd.read_csv(manual_path)
            self.genes_disgenet = df.to_dict("records")
        
        return self.genes_disgenet
    
    def fetch_omim(self) -> List[Dict[str, Any]]:
        """
        Fetch genes from OMIM.
        
        Note: OMIM requires API key for access.
        
        Returns:
            List of gene dictionaries
        """
        omim_id = self.config.get("disease.omim_id", "")
        
        if not omim_id:
            print("OMIM ID not specified.")
            return self.genes_omim
        
        print("OMIM requires API access.")
        print(f"Visit: https://www.omim.org/entry/{omim_id}")
        
        # Check for manual data
        manual_path = self.config.data_dir / "raw" / "omim_genes.csv"
        if manual_path.exists():
            df = pd.read_csv(manual_path)
            self.genes_omim = df.to_dict("records")
        
        return self.genes_omim
    
    def merge_genes(self) -> List[Dict[str, Any]]:
        """
        Merge genes from all sources and remove duplicates.
        
        Returns:
            Merged unique gene list
        """
        all_data = self.genes_opentargets + self.genes_genecards + self.genes_disgenet + self.genes_omim
        
        # Track unique genes
        seen_genes: Set[str] = set()
        self.all_genes = []
        
        for gene in all_data:
            symbol = gene.get("gene_symbol", "").upper()
            if symbol and symbol not in seen_genes:
                seen_genes.add(symbol)
                self.all_genes.append(gene)
        
        return self.all_genes
    
    def get_genes_df(self) -> pd.DataFrame:
        """
        Get merged genes as DataFrame.
        
        Returns:
            DataFrame with gene data
        """
        return pd.DataFrame(self.all_genes)
    
    def get_gene_symbols(self) -> List[str]:
        """
        Get list of unique gene symbols.
        
        Returns:
            List of gene symbols
        """
        return [g.get("gene_symbol", "").upper() for g in self.all_genes if g.get("gene_symbol")]
    
    def save_genes(self, filename: Optional[str] = None) -> Path:
        """
        Save genes to CSV.
        
        Args:
            filename: Optional custom filename
            
        Returns:
            Path to saved file
        """
        df = self.get_genes_df()
        
        if filename is None:
            filename = "disease_genes.csv"
        
        output_path = self.config.data_dir / "raw" / filename
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        df.to_csv(output_path, index=False)
        
        # Also save just the symbols
        symbols_path = self.config.data_dir / "raw" / "disease_gene_symbols.txt"
        symbols = self.get_gene_symbols()
        with open(symbols_path, "w") as f:
            f.write("\n".join(symbols))
        
        return output_path
    
    def load_genes(self, filename: str = "disease_genes.csv") -> pd.DataFrame:
        """
        Load previously saved genes.
        
        Args:
            filename: CSV filename
            
        Returns:
            DataFrame with gene data
        """
        input_path = self.config.data_dir / "raw" / filename
        
        if input_path.exists():
            df = pd.read_csv(input_path)
            self.all_genes = df.to_dict("records")
            return df
        else:
            raise FileNotFoundError(f"Disease genes file not found: {input_path}")


class ManualGeneDownload:
    """Helper class for manual gene data download."""
    
    def __init__(self, config: Config):
        self.config = config
    
    def print_instructions(self):
        """Print instructions for manual data download."""
        disease_name = self.config.disease_name
        disgenet_id = self.config.get("disease.disgenet_id", "")
        omim_id = self.config.get("disease.omim_id", "")
        
        print("\n" + "="*60)
        print("MANUAL DATA DOWNLOAD INSTRUCTIONS")
        print("="*60)
        
        print(f"\nDisease: {disease_name}\n")
        
        print("1. GeneCards")
        print(f"   URL: https://www.genecards.org/Search/Keyword?queryString={disease_name.replace(' ', '+')}")
        print(f"   Save as: {self.config.data_dir}/raw/genecards_genes.csv")
        print("   Required columns: gene_symbol, relevance_score")
        
        print("\n2. DisGeNET")
        if disgenet_id:
            print(f"   URL: https://www.disgenet.org/browser/0/1/0/{disgenet_id}/")
        print("   Or search at: https://www.disgenet.org/search")
        print(f"   Save as: {self.config.data_dir}/raw/disgenet_genes.csv")
        print("   Required columns: gene_symbol, score")
        
        print("\n3. OMIM (Optional)")
        if omim_id:
            print(f"   URL: https://www.omim.org/entry/{omim_id}")
        print(f"   Save as: {self.config.data_dir}/raw/omim_genes.csv")
        print("   Required columns: gene_symbol")
        
        print("\n" + "="*60)
