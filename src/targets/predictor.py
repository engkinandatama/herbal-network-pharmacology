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
        
        print(f"SMILES list saved to: {output_path}")
        print(f"\nInstructions:")
        print("1. Go to http://www.swisstargetprediction.ch/")
        print("2. Submit each SMILES one by one")
        print("3. Download results and place in data/processed/")
        print("4. Run the import function to parse results")
        
        return output_path
    
    def import_results(self, results_dir: str) -> pd.DataFrame:
        """
        Import manually downloaded results.
        
        Args:
            results_dir: Directory containing result files
            
        Returns:
            Combined DataFrame
        """
        results_path = Path(results_dir)
        all_results = []
        
        for file in results_path.glob("*.csv"):
            df = pd.read_csv(file)
            df["source_file"] = file.name
            all_results.append(df)
        
        if all_results:
            combined = pd.concat(all_results, ignore_index=True)
            
            # Save combined results
            output_path = self.config.data_dir / "processed" / "predicted_targets.csv"
            combined.to_csv(output_path, index=False)
            
            return combined
        
        return pd.DataFrame()
