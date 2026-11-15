import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
logger = logging.getLogger(__name__)

class DataCurator:
    """Module for molecular curation, standardization, and binarization using RDKit."""
    
    def __init__(self, df: pd.DataFrame, activity_threshold_nM: float = 1000.0):
        self.df = df.copy()
        self.activity_threshold_nM = activity_threshold_nM

    def clean_and_standardize(self) -> pd.DataFrame:
        """Clean data, standardize SMILES, calculate properties, and binarize activity."""
        if self.df.empty:
            logger.warning("Input dataframe is empty. Skipping curation.")
            return self.df

        logger.info("Starting chemical curation pipeline...")
        
        initial_count = len(self.df)
        self.df = self.df.dropna(subset=["canonical_smiles", "standard_value"])
        logger.info(f"Dropped {initial_count - len(self.df)} rows with missing values.")

        self.df["standard_value"] = pd.to_numeric(self.df["standard_value"], errors="coerce")
        self.df = self.df.dropna(subset=["standard_value"])

        canonical_smiles_list = []
        mw_list = []
        valid_indices = []

        for idx, row in self.df.iterrows():
            mol = Chem.MolFromSmiles(row["canonical_smiles"])
            if mol is not None:
                canon_smiles = Chem.MolToSmiles(mol, isomericSmiles=False)
                canonical_smiles_list.append(canon_smiles)
                mw_list.append(Descriptors.MolWt(mol))
                valid_indices.append(idx)

        filtered_df = self.df.loc[valid_indices].copy()
        filtered_df["canonical_smiles"] = canonical_smiles_list
        filtered_df["molecular_weight"] = mw_list

        filtered_df = filtered_df.drop_duplicates(subset=["canonical_smiles"])
        filtered_df["active"] = (filtered_df["standard_value"] <= self.activity_threshold_nM).astype(int)

        logger.info(f"Curation complete. Final dataset contains {len(filtered_df)} compounds.")
        return filtered_df
