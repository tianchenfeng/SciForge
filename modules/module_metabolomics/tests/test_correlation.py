import numpy as np
import pandas as pd

from modules.module_metabolomics.correlation import correlate_16s_metabolomics, correlate_features


def test_phenotype_and_lightweight_multiomics_correlations():
    samples = [f"S{i}" for i in range(8)]; trend = np.arange(8, dtype=float)
    abundance = pd.DataFrame([trend, trend[::-1]], index=["M1", "M2"], columns=samples)
    phenotype = pd.DataFrame({"Sample ID": samples, "pH": trend})
    long, corr, q = correlate_features(abundance, phenotype, "spearman")
    assert long.iloc[0]["FDR"] < .05 and corr.shape == (2, 1) and q.shape == (2, 1)
    taxa = pd.DataFrame([trend], index=["Taxon1"], columns=samples)
    multi = correlate_16s_metabolomics(taxa, abundance)
    assert set(["Taxon", "Feature ID", "Spearman", "FDR"]).issubset(multi.columns)
