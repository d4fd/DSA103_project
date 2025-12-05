import numpy as np
import pandas as pd
import random
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors
from IPython.display import display


#import custom functions
from helper_functions import take_most_common, populate_chemistry
from derive_chemistry import derive_chemistry



df = pd.read_csv ("data/raw_data/mtbs_tropical_annotations.tsv", sep = '\t')

#compounds, rename some of the cols
compounds = df.rename(columns={'structure_smiles':'smiles', "structure_taxonomy_npclassifier_01pathway":"class","structure_taxonomy_npclassifier_03class": "my_class"}).groupby('smiles',                    
  sort=False).agg(lambda s: s.dropna().mode().iat[0] if not s.dropna().mode().empty else             
  np.nan).reset_index()

#generate fake presAbs data for testing purposes (would normally be generated in previous part of pipeline)
samples = ["sample1","sample2","sample3","sample4","sample5"]
presAbs = pd.DataFrame(0, index=samples, columns=compounds["structure_inchikey"])
for s in samples:
    for cmpd in compounds["structure_inchikey"]:
            presAbs.loc[s, cmpd] = random.randint(0,1) #0 symbolizes absence, 1 presence
mtbs = {"compounds":compounds, "presAbs":presAbs}

results = derive_chemistry(mtbs)

output_dir = Path("data/processed")
output_dir.mkdir(parents=True, exist_ok=True)

for name, result_df in results.items():
    result_df.to_csv(output_dir / f"{name}.csv")
