"""
Configuration Loader Module
===========================

Handles loading and merging of YAML configuration files.
"""

import os
import yaml
from pathlib import Path
from typing import Any, Dict, Optional
from copy import deepcopy


class Config:
    """Configuration manager for network pharmacology analysis."""
    
    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize configuration.
        
        Args:
            config_path: Path to project-specific YAML config file.
                        If None, only default config is loaded.
        """
        self._config: Dict[str, Any] = {}
        self._project_root = self._find_project_root()
        
        # Load default config first
        self._load_default_config()
        
        # Merge with project config if provided
        if config_path:
            self._load_project_config(config_path)
        
        # Ensure output directories exist
        self._setup_directories()
    
    def _find_project_root(self) -> Path:
        """Find the project root directory."""
        current = Path(__file__).resolve().parent
        while current != current.parent:
            if (current / "config").exists() or (current / "requirements.txt").exists():
                return current
            current = current.parent
        return Path.cwd()
    
    def _load_default_config(self) -> None:
        """Load the default configuration file."""
        default_path = self._project_root / "config" / "default.yaml"
        if default_path.exists():
            with open(default_path, "r", encoding="utf-8") as f:
                self._config = yaml.safe_load(f) or {}
    
    def _load_project_config(self, config_path: str) -> None:
        """
        Load and merge project-specific configuration.
        
        Args:
            config_path: Path to project config file.
        """
        path = Path(config_path)
        if not path.is_absolute():
            path = self._project_root / path
        
        if not path.exists():
            raise FileNotFoundError(f"Config file not found: {path}")
        
        with open(path, "r", encoding="utf-8") as f:
            project_config = yaml.safe_load(f) or {}
        
        # Deep merge project config into default
        self._config = self._deep_merge(self._config, project_config)
    
    def _deep_merge(self, base: Dict, override: Dict) -> Dict:
        """
        Deep merge two dictionaries.
        
        Args:
            base: Base dictionary (default config)
            override: Override dictionary (project config)
            
        Returns:
            Merged dictionary
        """
        result = deepcopy(base)
        
        for key, value in override.items():
            if key in result and isinstance(result[key], dict) and isinstance(value, dict):
                result[key] = self._deep_merge(result[key], value)
            else:
                result[key] = deepcopy(value)
        
        return result
    
    def _setup_directories(self) -> None:
        """Create output directories if they don't exist."""
        output = self.get("output", {})
        
        dirs_to_create = [
            output.get("data_dir", "data"),
            output.get("figures_dir", "outputs/figures"),
            output.get("tables_dir", "outputs/tables"),
            "logs",
        ]
        
        for dir_path in dirs_to_create:
            full_path = self._project_root / dir_path
            full_path.mkdir(parents=True, exist_ok=True)
            
            # Create subdirectories for data
            if "data" in dir_path:
                (full_path / "raw").mkdir(exist_ok=True)
                (full_path / "processed").mkdir(exist_ok=True)
                (full_path / "results").mkdir(exist_ok=True)
    
    def get(self, key: str, default: Any = None) -> Any:
        """
        Get a configuration value.
        
        Args:
            key: Dot-notation key (e.g., "plant.name" or "analysis.hub_gene_top_n")
            default: Default value if key not found
            
        Returns:
            Configuration value
        """
        keys = key.split(".")
        value = self._config
        
        try:
            for k in keys:
                value = value[k]
            return value
        except (KeyError, TypeError):
            return default
    
    def get_path(self, key: str) -> Path:
        """
        Get a path configuration value as Path object.
        
        Args:
            key: Configuration key for a path value
            
        Returns:
            Path object resolved relative to project root
        """
        path_str = self.get(key, "")
        return self._project_root / path_str
    
    @property
    def project_root(self) -> Path:
        """Get the project root directory."""
        return self._project_root
    
    @property
    def plant_name(self) -> str:
        """Get the plant name."""
        return self.get("plant.name", "Unknown Plant")
    
    @property
    def plant_latin_name(self) -> str:
        """Get the plant Latin name."""
        return self.get("plant.latin_name", "")
    
    @property
    def disease_name(self) -> str:
        """Get the disease name."""
        return self.get("disease.name", "Unknown Disease")
    
    @property
    def known_compounds(self) -> list:
        """Get list of known compounds."""
        return self.get("plant.known_compounds", [])
    
    @property
    def data_dir(self) -> Path:
        """Get the data directory path."""
        return self.get_path("output.data_dir")
    
    @property
    def figures_dir(self) -> Path:
        """Get the figures directory path."""
        return self.get_path("output.figures_dir")
    
    @property
    def tables_dir(self) -> Path:
        """Get the tables directory path."""
        return self.get_path("output.tables_dir")
    
    def to_dict(self) -> Dict[str, Any]:
        """Get the full configuration as dictionary."""
        return deepcopy(self._config)
    
    def __repr__(self) -> str:
        project = self.get("project.name", "unnamed")
        return f"Config(project='{project}')"


def load_config(config_path: Optional[str] = None) -> Config:
    """
    Convenience function to load configuration.
    
    Args:
        config_path: Path to project config file
        
    Returns:
        Config object
    """
    return Config(config_path)
