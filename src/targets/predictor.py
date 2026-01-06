"""
Target Predictor Module
=======================

Predicts drug targets for compounds using:
- SwissTargetPrediction
- PharmMapper (optional)
"""

import pandas as pd
import requests
import time
from pathlib import Path
from typing import List, Dict, Optional, Any
from tqdm import tqdm
import json

from ..config_loader import Config


class TargetPredictor:
    """Predicts drug targets for compounds using various services."""
    
    def __init__(self, config: Config):
        """
        Initialize the predictor.
        
        Args:
            config: Configuration object
        """
        self.config = config
        self.compounds: List[Dict[str, Any]] = []
        self.predictions: List[Dict[str, Any]] = []
        self._request_delay = 2.0  # Rate limiting for web services
    
    def load_compounds(self, filename: str = "compounds.csv") -> pd.DataFrame:
        """
        Load compounds from saved CSV.
        
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
    
    def predict_swiss_target(self, smiles: str, compound_name: str = "") -> List[Dict[str, Any]]:
        """
        Predict targets using SwissTargetPrediction.
        
        Note: SwissTargetPrediction doesn't have a public API, so this
        simulates the web form submission. For actual use, you may need
        to use Selenium or submit manually.
        
        Args:
            smiles: SMILES string of the compound
            compound_name: Name of the compound
            
        Returns:
            List of predicted targets
        """
        # SwissTargetPrediction URL
        url = "http://www.swisstargetprediction.ch/predict.php"
        
        targets = []
        
        try:
            # Submit the form
            data = {
                "smiles": smiles,
                "organism": "Homo sapiens"
            }
            
            headers = {
                "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36",
                "Content-Type": "application/x-www-form-urlencoded"
            }
            
            response = requests.post(url, data=data, headers=headers, timeout=60)
            
            if response.status_code == 200:
                # Parse the response
                # Note: The actual parsing depends on the response format
                # This is a simplified version
                from bs4 import BeautifulSoup
                soup = BeautifulSoup(response.text, "html.parser")
                
                # Find the results table
                table = soup.find("table", {"id": "resultTable"})
                if table:
                    rows = table.find_all("tr")
                    for row in rows[1:]:  # Skip header
                        cols = row.find_all("td")
                        if len(cols) >= 4:
                            target_data = {
                                "compound_name": compound_name,
                                "compound_smiles": smiles,
                                "target_name": cols[0].get_text(strip=True),
                                "target_uniprot": cols[1].get_text(strip=True),
                                "gene_symbol": cols[2].get_text(strip=True),
                                "probability": float(cols[3].get_text(strip=True)) if cols[3].get_text(strip=True) else 0,
                                "source": "SwissTargetPrediction"
                            }
                            targets.append(target_data)
            
        except Exception as e:
            print(f"Warning: SwissTargetPrediction failed for {compound_name}: {e}")
            print("Consider using the web interface manually at: http://www.swisstargetprediction.ch/")
        
        return targets
    
    def predict_all(self) -> List[Dict[str, Any]]:
        """
        Predict targets for all compounds with SMILES.
        
        Returns:
            All predicted targets
        """
        probability_threshold = self.config.get("analysis.target_probability_threshold", 0.1)
        
        compounds_with_smiles = [c for c in self.compounds if c.get("smiles")]
        
        print(f"Predicting targets for {len(compounds_with_smiles)} compounds...")
        
        for compound in tqdm(compounds_with_smiles, desc="Predicting targets"):
            smiles = compound.get("smiles")
            name = compound.get("name", "Unknown")
            
            if smiles:
                targets = self.predict_swiss_target(smiles, name)
                
                # Filter by probability threshold
                filtered = [t for t in targets if t.get("probability", 0) >= probability_threshold]
                self.predictions.extend(filtered)
                
                time.sleep(self._request_delay)  # Rate limiting
        
        return self.predictions
    
    def get_predictions_df(self) -> pd.DataFrame:
        """
        Get predictions as DataFrame.
        
        Returns:
            DataFrame with target predictions
        """
        return pd.DataFrame(self.predictions)
    
    def get_unique_targets(self) -> List[str]:
        """
        Get list of unique target gene symbols.
        
        Returns:
            List of gene symbols
        """
        return list(set(p.get("gene_symbol") for p in self.predictions if p.get("gene_symbol")))
    
    def save_targets(self, filename: Optional[str] = None) -> Path:
        """
        Save predictions to CSV.
        
        Args:
            filename: Optional custom filename
            
        Returns:
            Path to saved file
        """
        df = self.get_predictions_df()
        
        if filename is None:
            filename = "predicted_targets.csv"
        
        output_path = self.config.data_dir / "processed" / filename
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        df.to_csv(output_path, index=False)
        
        # Also save unique targets
        unique_targets = self.get_unique_targets()
        targets_path = self.config.data_dir / "processed" / "unique_targets.txt"
        with open(targets_path, "w") as f:
            f.write("\n".join(unique_targets))
        
        return output_path


class ManualTargetPrediction:
    """
    Helper class for manual target prediction workflow.
    
    Since SwissTargetPrediction doesn't have a reliable API,
    this class helps manage batch submission via web interface.
    """
    
    def __init__(self, config: Config):
        self.config = config
    
    def generate_smiles_list(self, output_file: str = "smiles_for_prediction.txt") -> Path:
        """
        Generate a list of SMILES for manual submission.
        
        Args:
            output_file: Output filename
            
        Returns:
            Path to saved file
        """
        compounds_path = self.config.data_dir / "raw" / "compounds.csv"
        
        if not compounds_path.exists():
            raise FileNotFoundError("Compounds file not found. Run 'collect' first.")
        
        df = pd.read_csv(compounds_path)
        smiles_df = df[df["smiles"].notna()][["name", "smiles"]]
        
        output_path = self.config.data_dir / "processed" / output_file
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, "w") as f:
            for _, row in smiles_df.iterrows():
                f.write(f"{row['name']}\t{row['smiles']}\n")
        
        return output_path


class SwissTargetParser:
    """
    Parser for SwissTargetPrediction results.
    
    Handles both:
    - Copy-pasted table data (txt format)
    - Downloaded CSV files
    """
    
    def __init__(self, config: Config):
        self.config = config
        self.all_targets: List[Dict[str, Any]] = []
        self._probability_threshold = config.get("analysis.target_probability_threshold", 0.1)
    
    def parse_results(self, filepath: str) -> List[Dict[str, Any]]:
        """
        Parse SwissTargetPrediction results from file.
        
        Args:
            filepath: Path to results file (txt or csv)
            
        Returns:
            List of target dictionaries
        """
        filepath = Path(filepath)
        
        if filepath.suffix.lower() == ".csv":
            return self._parse_csv(filepath)
        else:
            return self._parse_txt(filepath)
    
    def _parse_txt(self, filepath: Path) -> List[Dict[str, Any]]:
        """Parse copy-pasted txt results."""
        current_compound = None
        
        with open(filepath, "r", encoding="utf-8") as f:
            lines = f.readlines()
        
        for line in lines:
            line = line.strip()
            
            # Skip empty lines
            if not line:
                continue
            
            # Check for compound header (e.g., "SwissTargetPrediction Quercetin")
            if "SwissTargetPrediction" in line:
                parts = line.replace("SwissTargetPrediction", "").strip()
                current_compound = parts if parts else None
                continue
            
            # Also check for compound name pattern at end of header
            if line.endswith("SwissTargetPrediction"):
                parts = line.replace("SwissTargetPrediction", "").strip()
                current_compound = parts if parts else None
                continue
            
            # Skip header row
            if line.startswith("Target\t") or line.startswith("Target	"):
                continue
            
            # Parse data row (tab-separated)
            parts = line.split("\t")
            if len(parts) >= 6:
                try:
                    target_full = parts[0]
                    gene_symbol = parts[1]
                    uniprot_id = parts[2]
                    chembl_id = parts[3] if len(parts) > 3 else ""
                    target_class = parts[4] if len(parts) > 4 else ""
                    probability = float(parts[5]) if len(parts) > 5 else 0.0
                    
                    # Apply threshold filter
                    if probability >= self._probability_threshold:
                        self.all_targets.append({
                            "compound_name": current_compound or "Unknown",
                            "target_name": target_full,
                            "gene_symbol": gene_symbol,
                            "uniprot_id": uniprot_id,
                            "chembl_id": chembl_id,
                            "target_class": target_class,
                            "probability": probability,
                            "source": "SwissTargetPrediction"
                        })
                except (ValueError, IndexError):
                    continue
        
        return self.all_targets
    
    def _parse_csv(self, filepath: Path) -> List[Dict[str, Any]]:
        """Parse CSV downloaded from SwissTargetPrediction."""
        df = pd.read_csv(filepath)
        
        # Get compound name from filename
        compound_name = filepath.stem.replace("_targets", "").replace("-", " ").title()
        
        for _, row in df.iterrows():
            probability = row.get("Probability*", row.get("Probability", 0))
            
            if probability >= self._probability_threshold:
                self.all_targets.append({
                    "compound_name": compound_name,
                    "target_name": row.get("Target", ""),
                    "gene_symbol": row.get("Common name", ""),
                    "uniprot_id": row.get("Uniprot ID", ""),
                    "chembl_id": row.get("ChEMBL ID", ""),
                    "target_class": row.get("Target Class", ""),
                    "probability": probability,
                    "source": "SwissTargetPrediction"
                })
        
        return self.all_targets
    
    def save_targets(self, filename: str = "predicted_targets.csv") -> Path:
        """Save parsed targets to CSV."""
        output_path = self.config.data_dir / "processed" / filename
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        df = pd.DataFrame(self.all_targets)
        df.to_csv(output_path, index=False)
        
        # Also save unique targets list
        unique_targets = df["gene_symbol"].dropna().unique().tolist()
        targets_file = output_path.parent / "unique_targets.txt"
        with open(targets_file, "w") as f:
            f.write("\n".join(unique_targets))
        
        return output_path


    """
    Automated target prediction using STITCH database.
    
    STITCH (Search Tool for Interactions of Chemicals) is a sister database
    of STRING and provides chemical-protein interactions with confidence scores.
    
    API Documentation: http://stitch.embl.de/cgi/help.pl?UserId=abcdef&sessionId=123
    """
    
    STITCH_API_URL = "http://stitch.embl.de/api"
    
    def __init__(self, config: Config):
        """
        Initialize STITCH predictor.
        
        Args:
            config: Configuration object
        """
        self.config = config
        self.compounds: List[Dict[str, Any]] = []
        self.predictions: List[Dict[str, Any]] = []
        self._species = 9606  # Homo sapiens
        self._score_threshold = config.get("analysis.target_probability_threshold", 0.4) * 1000  # STITCH uses 0-1000
        self._request_delay = 1.0
    
    def load_compounds(self, filename: str = "compounds.csv") -> pd.DataFrame:
        """Load compounds from saved CSV."""
        input_path = self.config.data_dir / "raw" / filename
        
        if input_path.exists():
            df = pd.read_csv(input_path)
            self.compounds = df.to_dict("records")
            return df
        else:
            raise FileNotFoundError(f"Compounds file not found: {input_path}")
    
    def _query_stitch_by_name(self, compound_name: str) -> List[Dict[str, Any]]:
        """
        Query STITCH by compound name.
        
        Args:
            compound_name: Name of the compound
            
        Returns:
            List of interaction data
        """
        url = f"{self.STITCH_API_URL}/json/interactors"
        
        params = {
            "identifiers": compound_name,
            "species": self._species,
            "limit": 100,
        }
        
        try:
            response = requests.get(url, params=params, timeout=60)
            
            if response.status_code == 200:
                return response.json()
            else:
                return []
                
        except Exception as e:
            print(f"STITCH query failed for {compound_name}: {e}")
            return []
    
    def _query_stitch_interactions(self, compound_name: str) -> List[Dict[str, Any]]:
        """
        Get protein interactions for a compound from STITCH.
        
        Args:
            compound_name: Name of the compound
            
        Returns:
            List of target dictionaries
        """
        url = f"{self.STITCH_API_URL}/json/interaction_partners"
        
        params = {
            "identifiers": compound_name,
            "species": self._species,
            "required_score": int(self._score_threshold),
            "limit": 50,
        }
        
        targets = []
        
        try:
            response = requests.get(url, params=params, timeout=60)
            
            if response.status_code == 200:
                data = response.json()
                
                for item in data:
                    # STITCH returns both chemical-protein and protein-protein
                    # We only want chemical-protein interactions
                    string_id_a = item.get("stringId_A", "")
                    string_id_b = item.get("stringId_B", "")
                    
                    # Chemical IDs start with "CID" or "s" prefix in STITCH
                    if string_id_a.startswith("CID") or string_id_a.startswith("s"):
                        protein_id = string_id_b
                        protein_name = item.get("preferredName_B", "")
                    elif string_id_b.startswith("CID") or string_id_b.startswith("s"):
                        protein_id = string_id_a
                        protein_name = item.get("preferredName_A", "")
                    else:
                        continue  # Skip protein-protein interactions
                    
                    target_data = {
                        "compound_name": compound_name,
                        "target_name": protein_name,
                        "gene_symbol": protein_name.upper(),
                        "string_id": protein_id,
                        "score": item.get("score", 0) / 1000,  # Normalize to 0-1
                        "source": "STITCH"
                    }
                    targets.append(target_data)
                    
        except Exception as e:
            print(f"STITCH interaction query failed for {compound_name}: {e}")
        
        return targets
    
    def predict_all(self) -> List[Dict[str, Any]]:
        """
        Predict targets for all compounds using STITCH.
        
        Returns:
            All predicted targets
        """
        if not self.compounds:
            self.load_compounds()
        
        print(f"Querying STITCH for {len(self.compounds)} compounds...")
        
        for compound in tqdm(self.compounds, desc="STITCH prediction"):
            name = compound.get("name", "")
            
            if name:
                targets = self._query_stitch_interactions(name)
                self.predictions.extend(targets)
                time.sleep(self._request_delay)  # Rate limiting
        
        print(f"Found {len(self.predictions)} compound-target interactions")
        return self.predictions
    
    def get_predictions_df(self) -> pd.DataFrame:
        """Get predictions as DataFrame."""
        return pd.DataFrame(self.predictions)
    
    def get_unique_targets(self) -> List[str]:
        """Get list of unique target gene symbols."""
        return list(set(p.get("gene_symbol") for p in self.predictions if p.get("gene_symbol")))
    
    def save_targets(self, filename: Optional[str] = None) -> Path:
        """Save predictions to CSV."""
        df = self.get_predictions_df()
        
        if filename is None:
            filename = "predicted_targets_stitch.csv"
        
        output_path = self.config.data_dir / "processed" / filename
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        df.to_csv(output_path, index=False)
        
        # Also save unique targets
        unique_targets = self.get_unique_targets()
        targets_path = self.config.data_dir / "processed" / "unique_targets.txt"
        with open(targets_path, "w") as f:
            f.write("\n".join(unique_targets))
        
        return output_path


class CombinedTargetPredictor:
    """
    Combines multiple target prediction sources.
    
    Uses STITCH (automated) as primary source, with option to
    merge with manually obtained SwissTargetPrediction results.
    """
    
    def __init__(self, config: Config):
        self.config = config
        self.stitch_predictor = STITCHPredictor(config)
        self.all_predictions: List[Dict[str, Any]] = []
    
    def predict_with_stitch(self) -> List[Dict[str, Any]]:
        """Run automated prediction with STITCH."""
        self.stitch_predictor.load_compounds()
        predictions = self.stitch_predictor.predict_all()
        self.all_predictions.extend(predictions)
        return predictions
    
    def merge_swiss_results(self, swiss_csv_path: str) -> pd.DataFrame:
        """
        Merge manually obtained SwissTargetPrediction results.
        
        Args:
            swiss_csv_path: Path to SwissTargetPrediction results CSV
            
        Returns:
            Combined DataFrame
        """
        swiss_df = pd.read_csv(swiss_csv_path)
        swiss_df["source"] = "SwissTargetPrediction"
        
        stitch_df = pd.DataFrame(self.all_predictions)
        
        combined = pd.concat([stitch_df, swiss_df], ignore_index=True)
        
        # Remove duplicates (same compound-target pair)
        combined = combined.drop_duplicates(
            subset=["compound_name", "gene_symbol"],
            keep="first"
        )
        
        return combined
    
    def get_unique_targets(self) -> List[str]:
        """Get all unique targets from all sources."""
        return list(set(p.get("gene_symbol") for p in self.all_predictions if p.get("gene_symbol")))
    
    def save_all(self, filename: str = "predicted_targets.csv") -> Path:
        """Save all predictions to CSV."""
        df = pd.DataFrame(self.all_predictions)
        
        output_path = self.config.data_dir / "processed" / filename
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        df.to_csv(output_path, index=False)
        
        # Save unique targets
        unique_targets = self.get_unique_targets()
        targets_path = self.config.data_dir / "processed" / "unique_targets.txt"
        with open(targets_path, "w") as f:
            f.write("\n".join(unique_targets))
        
        return output_path
