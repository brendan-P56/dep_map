# reading in data proteomics
protein_quant_current_normalized.csv <- read.csv("./data/in/protein_quant_current_normalized.csv.gz", header=T)
Table_S1_Sample_Information <- read.csv("./data/in/Table_S1_Sample_Information.csv")

# reading in dep data
Model <- read.csv("./data/in/Model.csv")
CRISPRGeneDependency <- read.csv("./data/in/CRISPRGeneDependency.csv.gz")

# subtracting the differences and chekcing assumption
Model <- Model[-which(!Model$ModelID %in% CRISPRGeneDependency$ModelID), ]
identical(CRISPRGeneDependency$ModelID, Model$ModelID)

