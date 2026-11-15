import pandas as pd
from chembl_webresource_client.new_client import new_client
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
logger = logging.getLogger(__name__)

class ChEMBLDataFetcher:
    """Module for automated extraction of bioactivity data from ChEMBL REST API."""
    
    def __init__(self, target_chembl_id: str, standard_type: str = "IC50"):
        self.target_chembl_id = target_chembl_id
        self.standard_type = standard_type

    def fetch_activities(self) -> pd.DataFrame:
        """Fetch valid bioactivity records (Binding assays) for the target."""
        logger.info(f"Connecting to ChEMBL for target {self.target_chembl_id} ({self.standard_type})...")
        
        activity = new_client.activity
        res = activity.filter(
            target_chembl_id=self.target_chembl_id,
            standard_type=self.standard_type,
            standard_relation="=",
            assay_type="B"
        ).only([
            "molecule_chembl_id", 
            "canonical_smiles", 
            "standard_type", 
            "standard_value", 
            "standard_units"
        ])

        df = pd.DataFrame(list(res))
        if df.empty:
            logger.warning("No data found for the specified parameters.")
            return df
            
        logger.info(f"Successfully retrieved {len(df)} raw records from ChEMBL.")
        return df
