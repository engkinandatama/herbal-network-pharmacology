"""
Compound Collector Module
=========================

Collects compound data from various databases including:
- Literature (from config)
- PubChem
- KNApSAcK
"""

import pandas as pd
import requests
import time
from pathlib import Path
from typing import List, Dict, Optional, Any
from tqdm import tqdm

try:
    import pubchempy as pcp
except ImportError:
    pcp = None

from ..config_loader import Config


class CompoundCollector:
    """Collects and enriches compound data from multiple sources."""
    
    def __init__(self, config: Config):
        """
        Initialize the collector.
        
        Args:
            config: Configuration object
        """
        self.config = config
        self.compounds: List[Dict[str, Any]] = []
        self._pubchem_delay = 0.5  # Rate limiting
    
    def collect_from_literature(self) -> List[Dict[str, Any]]:
        """
        Collect compounds from literature (defined in config).
        
        Returns:
            List of compound dictionaries
        """
        known_compounds = self.config.known_compounds
        
        for compound in known_compounds:
            compound_data = {
                "name": compound.get("name", ""),
                "pubchem_cid": compound.get("pubchem_cid"),
                "smiles": compound.get("smiles"),
                "source": "literature",
                "plant": self.config.plant_latin_name,
            }
            self.compounds.append(compound_data)
        
        return self.compounds
    
    def enrich_from_pubchem(self) -> List[Dict[str, Any]]:
        """
        Enrich compound data with PubChem information.
        Fetches SMILES, molecular weight, formula, etc.
        
        Returns:
            Enriched compound list
        """
        if pcp is None:
            print("Warning: pubchempy not installed. Skipping PubChem enrichment.")
            return self.compounds
        
        for compound in tqdm(self.compounds, desc="Enriching from PubChem"):
            try:
                # Try by CID first
                if compound.get("pubchem_cid"):
                    cid = compound["pubchem_cid"]
                    pcp_compound = pcp.Compound.from_cid(cid)
                # Try by name
                elif compound.get("name"):
                    results = pcp.get_compounds(compound["name"], "name")
                    if results:
                        pcp_compound = results[0]
                        compound["pubchem_cid"] = pcp_compound.cid
                    else:
                        continue
                else:
                    continue
                
                # Extract data
                if not compound.get("smiles"):
                    compound["smiles"] = pcp_compound.isomeric_smiles or pcp_compound.canonical_smiles
                
                compound["molecular_formula"] = pcp_compound.molecular_formula
                compound["molecular_weight"] = pcp_compound.molecular_weight
                compound["xlogp"] = pcp_compound.xlogp
                compound["hbd"] = pcp_compound.h_bond_donor_count
                compound["hba"] = pcp_compound.h_bond_acceptor_count
                compound["tpsa"] = pcp_compound.tpsa
                compound["rotatable_bonds"] = pcp_compound.rotatable_bond_count
                
                time.sleep(self._pubchem_delay)  # Rate limiting
                
            except Exception as e:
                print(f"Warning: Could not enrich {compound.get('name')}: {e}")
        
        return self.compounds
    
    def search_knapsack(self) -> List[Dict[str, Any]]:
        """
        Search KNApSAcK database for additional compounds.
        
        Note: KNApSAcK doesn't have a public API, so this uses web scraping
        or falls back to searching by plant keywords.
        
        Returns:
            List of new compounds found
        """
        # KNApSAcK base URL
        base_url = "http://www.knapsackfamily.com/knapsack_core/result.php"
        
        search_keywords = self.config.get("plant.search_keywords", [])
        latin_name = self.config.plant_latin_name
        
        new_compounds = []
        
        # Search by Latin name
        try:
            params = {
                "sname": "organism",
                "word": latin_name
            }
            
            response = requests.get(base_url, params=params, timeout=30)
            
            if response.status_code == 200:
                # Parse the HTML response
                from bs4 import BeautifulSoup
                soup = BeautifulSoup(response.text, "html.parser")
                
                # Find compound entries in the table
                tables = soup.find_all("table")
                for table in tables:
                    rows = table.find_all("tr")
                    for row in rows[1:]:  # Skip header
                        cols = row.find_all("td")
                        if len(cols) >= 3:
                            try:
                                knapsack_id = cols[0].get_text(strip=True)
                                name = cols[1].get_text(strip=True)
                                
                                # Check if already exists
                                existing_names = [c.get("name", "").lower() for c in self.compounds]
                                if name.lower() not in existing_names:
                                    new_compound = {
                                        "name": name,
                                        "knapsack_id": knapsack_id,
                                        "source": "knapsack",
                                        "plant": latin_name,
                                    }
                                    new_compounds.append(new_compound)
                                    self.compounds.append(new_compound)
                            except Exception:
                                continue
                
        except Exception as e:
            print(f"Warning: KNApSAcK search failed: {e}")
            print("You may need to manually search at: http://www.knapsackfamily.com/knapsack_core/")
        
        return new_compounds
    
    def fetch_smiles(self, compound_name: str) -> Optional[str]:
        """
        Fetch SMILES for a compound by name.
        
        Args:
            compound_name: Name of the compound
            
        Returns:
            SMILES string or None
        """
        if pcp is None:
            return None
        
        try:
            results = pcp.get_compounds(compound_name, "name")
            if results:
                return results[0].isomeric_smiles or results[0].canonical_smiles
        except Exception:
            pass
        
        return None
    
    def get_compounds_df(self) -> pd.DataFrame:
        """
        Get compounds as a pandas DataFrame.
        
        Returns:
            DataFrame with compound data
        """
        return pd.DataFrame(self.compounds)
    
    def filter_by_smiles(self) -> List[Dict[str, Any]]:
        """
        Filter to only compounds with valid SMILES.
        
        Returns:
            Filtered compound list
        """
        return [c for c in self.compounds if c.get("smiles")]
    
    def save_compounds(self, filename: Optional[str] = None) -> Path:
        """
        Save compounds to CSV file.
        
        Args:
            filename: Optional custom filename
            
        Returns:
            Path to saved file
        """
        df = self.get_compounds_df()
        
        if filename is None:
            filename = "compounds.csv"
        
        output_path = self.config.data_dir / "raw" / filename
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        df.to_csv(output_path, index=False)
        
        return output_path
    
    def load_compounds(self, filename: str = "compounds.csv") -> pd.DataFrame:
        """
        Load previously saved compounds.
        
        Args:
            filename: CSV filename
            
        Returns:
            DataFrame with compound data
        """
        input_path = self.config.data_dir / "raw" / filename
        
        if input_path.exists():
            df = pd.read_csv(input_path)
            self.compounds = df.to_dict("records")
            return df
        else:
            raise FileNotFoundError(f"Compounds file not found: {input_path}")
    
    def __len__(self) -> int:
        return len(self.compounds)
    
    def __repr__(self) -> str:
        return f"CompoundCollector(n_compounds={len(self.compounds)})"
