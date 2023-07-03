# PD_Hallucinations
Code repository for the paper - 'Abnormal higher-order network interactions in Parkinson's disease visual hallucinations' (Tan et al. 2023)

Code - contains the gradients and t-SNE analysis code. Gradients code was adapted from the Brainspace toolbox (https://github.com/MICA-MNI/BrainSpace).

Templates - contains the regional base templates used in the analysis. Templates include the Margulies gradients in Schafer-400 space from Margulies et al. 2016 (https://www.pnas.org/doi/10.1073/pnas.1608282113) and the Thomas Yeo's 7-Network Atlas (2011) in the Schaefer-400 space (DOI: 10.1152/jn.00338.2011; DOI: 10.1093/cercor/bhx179). We also have included the Anterior-Posterior axis and the Network hierarchy described in Zarkali et al. 2021 (https://doi.org/10.1038/s42003-020-01622-9)

Data - subject functional connectivity data in Schaefer-400 space. After preprocessing and denoising, subject data was bandpass filtered, z-scored and extracted into Schaefer-400 space.
