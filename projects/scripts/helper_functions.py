import numpy as np
import pandas as pd

def take_most_common(s: pd.Series):
    # same idea as R::which.max(table(x))
    counts = s.value_counts(dropna=False)
    return counts.idxmax()

def populate_chemistry(
    mtbs: pd.DataFrame,
    trait: pd.Series
) -> pd.Series:
    # Populate per-species chemistry by averaging traits over present compounds
    trait_series = pd.Series(trait, index=mtbs.columns)
    values = []
    for _, row in mtbs.iterrows():
        present = row == 1
        if present.any():
            vals = trait_series[present].dropna()
            if vals.empty:
                values.append(np.nan)
                continue
            values.append(vals.mean())
        else:
            values.append(np.nan)

    return pd.Series(values, index=mtbs.index)
