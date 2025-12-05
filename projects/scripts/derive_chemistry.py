import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from helper_functions import take_most_common, populate_chemistry


def derive_chemistry(mtbs):# Am schluss folgenden code noch einpacken in function
    isPresent = mtbs["compounds"]['smiles'].notna() & (mtbs["compounds"]['smiles'] != "")
    compoundsOK = mtbs["compounds"].loc[isPresent]  

    # we want all columns from 'class' to 'my_class' (inclusive)
    cols = list(compoundsOK.columns)
    start = cols.index('class')
    end = cols.index('my_class')
    cols_to_agg = cols[start:end + 1]

    compounds = (
        compoundsOK
        .groupby("smiles")[cols_to_agg]
        .agg(take_most_common)    # collapse each group to most common values
        .reset_index()            # `rownames_to_column("smiles")` analogue
    )

    #collapse metabolite data to smiles level
    # select columns
    met = mtbs["presAbs"].iloc[:, np.asarray(isPresent)]
    # transpose and convert to DataFrame
    met = met.T
    # add smiles information
    met["smiles"] = compoundsOK["smiles"].values
    # group by smiles and take mean across all columns
    metabolites = (met.groupby("smiles", as_index=False).mean())


    ## Prepare data tables
    # select smiles column only
    compReady = metabolites[["smiles"]]
    # left join to compounds
    compReady = compReady.merge(compounds, on="smiles", how="left")

    # generate SID (=compound IDs)
    compReady["SID"] = ["S" + str(i) for i in range(1, len(compReady) + 1)]
    # reorder columns
    cols = ["SID", "smiles"] + list(compounds.columns[compounds.columns.get_loc("class"):compounds.columns.get_loc("my_class") + 1])
    compReady = compReady[cols]

    mtbsReady = metabolites.drop(columns=["smiles"]).transpose()
    mtbsReady.columns = compReady["SID"].values
    mtbsReady[mtbsReady > 0] = 1


    ## Derive chemistry
    smilesOnly = compReady.set_index("SID")["smiles"].to_dict()
    smilesParsed = {sid: Chem.MolFromSmiles(sm) for sid, sm in smilesOnly.items()}

    redundant = [2, 7, 8, 11, 15, 17, 18, 20, 21, 24, 29] + list(range(33, 39)) + [41] + list(range(43, 46))
    all_desc = [d[0] for d in Descriptors.descList]
    descriptors = [d for i, d in enumerate(all_desc) if i not in redundant]

    rawChem = ({
        sid: {name: Descriptors.__dict__[name](mol) for name in descriptors}
        for sid, mol in smilesParsed.items()
    })

    rawChem_df = pd.DataFrame.from_dict(rawChem, orient="index").reset_index().rename(columns={"index": "SID"})



    compOut = (
        compReady
        .merge(rawChem_df, on = "SID", how="left")
    )
    compOut = compOut[["SID"] + list(compOut.loc[:, "smiles":"my_class"].columns) + list(compOut.loc[:, "MaxAbsEStateIndex":"fr_urea"].columns)]
    compOut.index = compOut["SID"]
    
    ## Populate chemistry across presence-absence matrix ----

    chemSpecies = pd.DataFrame({
        col: populate_chemistry(mtbsReady, compOut[col])
        for col in compOut.loc[:, "MaxAbsEStateIndex":"fr_urea"].columns
    })

    #delete species where all absent (shouldnt exist)
    chemSpecies = chemSpecies.loc[:, chemSpecies.sum(axis=0) > 0]

    out = {"presAbs" :mtbsReady, "chemistry": chemSpecies, "compounds":compOut}
    return out



