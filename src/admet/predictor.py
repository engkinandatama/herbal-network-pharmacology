"""
ADMET Predictor Module
======================

Predicts ADMET properties and evaluates drug-likeness.
"""

import pandas as pd
import requests
from pathlib import Path
from typing import List, Dict, Optional, Any
from tqdm import tqdm
import time

from ..config_loader import Config


class ADMETPredictor:
    """Predicts ADMET properties for compounds."""
    
    def __init__(self, config: Config):
        """
        Initialize the predictor.
        
        Args:
            config: Configuration object
        """
        self.config = config
        self.compounds: List[Dict[str, Any]] = []
        self.predictions: List[Dict[str, Any]] = []
        self._lipinski_max_violations = config.get("analysis.lipinski_violations_max", 1)
    
    def load_compounds(self, filename: str = "compounds.csv") -> pd.DataFrame:
        """
        Load compounds from file.
        
        Args:
            filename: CSV filename
            
        Returns:
            DataFrame with compound data
        """
        path = self.config.data_dir / "raw" / filename
        
        if path.exists():
            df = pd.read_csv(path)
            self.compounds = df.to_dict("records")
            return df
        else:
            raise FileNotFoundError(f"Compounds file not found: {path}")
    
    def calculate_lipinski(self, compound: Dict[str, Any]) -> Dict[str, Any]:
        """
        Calculate Lipinski's Rule of Five.
        
        Args:
            compound: Compound dictionary with properties
            
        Returns:
            Dictionary with Lipinski evaluation
        """
        violations = 0
        details = {}
        
        # Molecular weight <= 500
        mw = compound.get("molecular_weight", 0)
        if mw and mw > 500:
            violations += 1
            details["MW"] = f"{mw:.2f} (>500)"
        else:
            details["MW"] = f"{mw:.2f}" if mw else "N/A"
        
        # LogP <= 5
        logp = compound.get("xlogp", 0)
        if logp and logp > 5:
            violations += 1
            details["LogP"] = f"{logp:.2f} (>5)"
        else:
            details["LogP"] = f"{logp:.2f}" if logp else "N/A"
        
        # H-bond donors <= 5
        hbd = compound.get("hbd", 0)
        if hbd and hbd > 5:
            violations += 1
            details["HBD"] = f"{hbd} (>5)"
        else:
            details["HBD"] = str(hbd) if hbd else "N/A"
        
        # H-bond acceptors <= 10
        hba = compound.get("hba", 0)
        if hba and hba > 10:
            violations += 1
            details["HBA"] = f"{hba} (>10)"
        else:
            details["HBA"] = str(hba) if hba else "N/A"
        
        return {
            "lipinski_violations": violations,
            "lipinski_pass": violations <= self._lipinski_max_violations,
            "lipinski_details": details
        }
    
    def calculate_veber(self, compound: Dict[str, Any]) -> Dict[str, Any]:
        """
        Calculate Veber's rules for oral bioavailability.
        
        Args:
            compound: Compound dictionary with properties
            
        Returns:
            Dictionary with Veber evaluation
        """
        violations = 0
        details = {}
        
        # Rotatable bonds <= 10
        rotb = compound.get("rotatable_bonds", 0)
        if rotb and rotb > 10:
            violations += 1
            details["RotB"] = f"{rotb} (>10)"
        else:
            details["RotB"] = str(rotb) if rotb else "N/A"
        
        # TPSA <= 140
        tpsa = compound.get("tpsa", 0)
        if tpsa and tpsa > 140:
            violations += 1
            details["TPSA"] = f"{tpsa:.2f} (>140)"
        else:
            details["TPSA"] = f"{tpsa:.2f}" if tpsa else "N/A"
        
        return {
            "veber_violations": violations,
            "veber_pass": violations == 0,
            "veber_details": details
        }
    
    def predict_all(self) -> List[Dict[str, Any]]:
        """
        Predict ADMET properties for all compounds.
        
        Returns:
            List of predictions
        """
        if not self.compounds:
            self.load_compounds()
        
        self.predictions = []
        
        for compound in tqdm(self.compounds, desc="Calculating ADMET"):
            prediction = {
                "name": compound.get("name", ""),
                "smiles": compound.get("smiles", ""),
                "molecular_weight": compound.get("molecular_weight"),
                "xlogp": compound.get("xlogp"),
                "hbd": compound.get("hbd"),
                "hba": compound.get("hba"),
                "tpsa": compound.get("tpsa"),
                "rotatable_bonds": compound.get("rotatable_bonds"),
            }
            
            # Calculate rules
            lipinski = self.calculate_lipinski(compound)
            veber = self.calculate_veber(compound)
            
            prediction.update(lipinski)
            prediction.update(veber)
            
            # Overall drug-likeness
            prediction["drug_like"] = lipinski["lipinski_pass"] and veber["veber_pass"]
            
            self.predictions.append(prediction)
        
        return self.predictions
    
    def get_predictions_df(self) -> pd.DataFrame:
        """
        Get predictions as DataFrame.
        
        Returns:
            DataFrame with predictions
        """
        if not self.predictions:
            self.predict_all()
        
        # Flatten details for DataFrame
        flat_predictions = []
        for p in self.predictions:
            flat = {k: v for k, v in p.items() if not k.endswith("_details")}
            flat_predictions.append(flat)
        
        return pd.DataFrame(flat_predictions)
    
    def get_drug_like_compounds(self) -> List[Dict[str, Any]]:
        """
        Get compounds that pass drug-likeness filters.
        
        Returns:
            List of drug-like compounds
        """
        return [p for p in self.predictions if p.get("drug_like")]
    
    def generate_summary(self) -> Dict[str, Any]:
        """
        Generate ADMET summary statistics.
        
        Returns:
            Dictionary with summary
        """
        if not self.predictions:
            self.predict_all()
        
        total = len(self.predictions)
        lipinski_pass = sum(1 for p in self.predictions if p.get("lipinski_pass"))
        veber_pass = sum(1 for p in self.predictions if p.get("veber_pass"))
        drug_like = sum(1 for p in self.predictions if p.get("drug_like"))
        
        return {
            "total_compounds": total,
            "lipinski_pass": lipinski_pass,
            "lipinski_pass_pct": f"{lipinski_pass/total*100:.1f}%" if total else "0%",
            "veber_pass": veber_pass,
            "veber_pass_pct": f"{veber_pass/total*100:.1f}%" if total else "0%",
            "drug_like": drug_like,
            "drug_like_pct": f"{drug_like/total*100:.1f}%" if total else "0%",
        }
    
    def save_predictions(self, filename: Optional[str] = None) -> Path:
        """
        Save predictions to file.
        
        Args:
            filename: Optional custom filename
            
        Returns:
            Path to saved file
        """
        df = self.get_predictions_df()
        
        if filename is None:
            filename = "admet_predictions.csv"
        
        output_path = self.config.data_dir / "results" / filename
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        df.to_csv(output_path, index=False)
        
        # Save summary
        summary = self.generate_summary()
        summary_path = output_path.parent / "admet_summary.json"
        import json
        with open(summary_path, "w") as f:
            json.dump(summary, f, indent=2)
        
        return output_path


class SwissADMEHelper:
    """
    Helper for SwissADME web service.
    
    SwissADME doesn't have a public API, so this helps with manual workflow.
    """
    
    def __init__(self, config: Config):
        self.config = config
    
    def generate_smiles_batch(self, output_file: str = "smiles_for_swissadme.txt") -> Path:
        """
        Generate SMILES list for SwissADME batch submission.
        
        Args:
            output_file: Output filename
            
        Returns:
            Path to saved file
        """
        compounds_path = self.config.data_dir / "raw" / "compounds.csv"
        
        if not compounds_path.exists():
            raise FileNotFoundError("Compounds file not found.")
        
        df = pd.read_csv(compounds_path)
        smiles_list = df[df["smiles"].notna()]["smiles"].tolist()
        
        output_path = self.config.data_dir / "processed" / output_file
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, "w") as f:
            f.write("\n".join(smiles_list))
        
        print(f"SMILES list saved to: {output_path}")
        print(f"\nInstructions:")
        print("1. Go to http://www.swissadme.ch/")
        print("2. Paste the SMILES list")
        print("3. Submit and download results")
        print(f"4. Save to: {self.config.data_dir}/results/swissadme_results.csv")
        
        return output_path
