####################################################
######### Merge MT and Expression z-scores #########
####################################################

library(dplyr)
library(tibble)

setwd("./stunning-fiesta_processed/")

mut_df <- read.csv("luad_tcga_pancancer_mut_coded_for_merge.csv")

setwd("..")
setwd("..")
setwd("./Heatmap/Th17 RNA-Seq Zscore Gene Expression Data/")
expression_z <- read.csv("th17_genes_mrna_seq_v2_rsem_zscores_ref_all_samples.csv", row.names = 1)

# Transpose z-score data #

express_z_transposed <- as.data.frame(t(expression_z))
express_z_transposed <- rownames_to_column(express_z_transposed, var = "Patient.ID")  

# Merge #

merged <- merge(express_z_transposed, mut_df) 

setwd("../Merged Expression and MT Data Th17/")
write.csv(merged, "merged_th17_zscore_and_mut_data.csv")
