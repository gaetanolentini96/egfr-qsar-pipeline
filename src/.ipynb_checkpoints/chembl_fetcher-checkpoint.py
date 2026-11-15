{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "ebc4605f-4ec8-4398-b271-d41243dac5f8",
   "metadata": {},
   "outputs": [],
   "source": [
    "import pandas as pd\n",
    "from chembl_webresource_client.new_client import new_client\n",
    "import logging\n",
    "\n",
    "logging.basicConfig(level=logging.INFO, format=\"%(asctime)s [%(levelname)s] %(message)s\")\n",
    "logger = logging.getLogger(__name__)\n",
    "\n",
    "class ChEMBLDataFetcher:\n",
    "    \"\"\"Module for automated extraction of bioactivity data from ChEMBL REST API.\"\"\"\n",
    "    \n",
    "    def __init__(self, target_chembl_id: str, standard_type: str = \"IC50\"):\n",
    "        self.target_chembl_id = target_chembl_id\n",
    "        self.standard_type = standard_type\n",
    "\n",
    "    def fetch_activities(self) -> pd.DataFrame:\n",
    "        \"\"\"Fetch valid bioactivity records (Binding assays) for the target.\"\"\"\n",
    "        logger.info(f\"Connecting to ChEMBL for target {self.target_chembl_id} ({self.standard_type})...\")\n",
    "        \n",
    "        activity = new_client.activity\n",
    "        res = activity.filter(\n",
    "            target_chembl_id=self.target_chembl_id,\n",
    "            standard_type=self.standard_type,\n",
    "            standard_relation=\"=\",\n",
    "            assay_type=\"B\"  # Binding assays only\n",
    "        ).only([\n",
    "            \"molecule_chembl_id\", \n",
    "            \"canonical_smiles\", \n",
    "            \"standard_type\", \n",
    "            \"standard_value\", \n",
    "            \"standard_units\"\n",
    "        ])\n",
    "\n",
    "        df = pd.DataFrame(list(res))\n",
    "        if df.empty:\n",
    "            logger.warning(\"No data found for the specified parameters.\")\n",
    "            return df\n",
    "            \n",
    "        logger.info(f\"Successfully retrieved {len(df)} raw records from ChEMBL.\")\n",
    "        return df"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.11.7"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
