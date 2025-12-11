"""
Configuration management for ONT target sequencing pipeline
"""

import os
import json
from typing import Dict, Any, Optional


class Config:
    """Configuration manager for the pipeline"""
    
    DEFAULT_CONFIG = {
        'preprocessing': {
            'min_quality': 7.0,
            'min_length': 100,
            'max_length': 10000,
            'trim_adapters': True,
            'adapter_sequences': []
        },
        'alignment': {
            'preset': 'map-ont',
            'threads': 4,
            'min_mapping_quality': 20
        },
        'variant_calling': {
            'min_coverage': 10,
            'min_variant_frequency': 0.2,
            'min_base_quality': 7
        },
        'coverage_analysis': {
            'min_coverage': 10,
            'report_format': 'bed'
        },
        'output': {
            'output_dir': './ont_targseq_output',
            'keep_intermediates': True
        }
    }
    
    def __init__(self, config_file: Optional[str] = None):
        """
        Initialize configuration
        
        Args:
            config_file: Path to JSON configuration file (optional)
        """
        self.config = self.DEFAULT_CONFIG.copy()
        
        if config_file and os.path.exists(config_file):
            self.load_config(config_file)
    
    def load_config(self, config_file: str):
        """Load configuration from JSON file"""
        with open(config_file, 'r') as f:
            user_config = json.load(f)
            self._update_config(self.config, user_config)
    
    def _update_config(self, base: Dict, update: Dict):
        """Recursively update configuration"""
        for key, value in update.items():
            if key in base and isinstance(base[key], dict) and isinstance(value, dict):
                self._update_config(base[key], value)
            else:
                base[key] = value
    
    def save_config(self, config_file: str):
        """Save configuration to JSON file"""
        with open(config_file, 'w') as f:
            json.dump(self.config, f, indent=2)
    
    def get(self, section: str, key: Optional[str] = None, default: Any = None) -> Any:
        """
        Get configuration value
        
        Args:
            section: Configuration section
            key: Configuration key (optional)
            default: Default value if not found
            
        Returns:
            Configuration value
        """
        if section not in self.config:
            return default
        
        if key is None:
            return self.config[section]
        
        return self.config[section].get(key, default)
    
    def set(self, section: str, key: str, value: Any):
        """
        Set configuration value
        
        Args:
            section: Configuration section
            key: Configuration key
            value: Value to set
        """
        if section not in self.config:
            self.config[section] = {}
        
        self.config[section][key] = value
    
    def create_default_config(self, output_file: str):
        """Create a default configuration file"""
        with open(output_file, 'w') as f:
            json.dump(self.DEFAULT_CONFIG, f, indent=2)
        
        return output_file
