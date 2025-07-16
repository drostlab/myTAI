# Global variables for R CMD check
# This file declares global variables used in NSE contexts to avoid 
# R CMD check warnings.

# Global variables used in dplyr/ggplot2 NSE contexts
utils::globalVariables(c(
    # Data frame column names used in dplyr operations
    "Stratum", "GeneID", "Expression", "Condition", "TXI", "Group", "Stage",
    "Id", "Sample", "Gene", "Angle", "Run", "Threshold", "Set", "Size", "Type",
    "PartialTAI", "Ratio", "ColourVar", "PC1", "PC2", "angle", "strata",
    "highlight", "label", "p_value", "stage_label", "reference", "stage_quantiles",
    "p_sel", "p_all", "log_obs_exp", "PS_num", "Rank", "AvgPS", "V1", "V2", 
    "group", "se", "y", "age", "sample_size", "pval",
    
    # Functions commonly used without explicit namespace
    "var", "sd", "quantile", "prcomp", "colorRampPalette", "setNames", "p.adjust"
))
