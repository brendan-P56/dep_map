# reading in data proteomics
protein_quant_current_normalized.csv <- read.csv("/media/chenglab-bp/aquila/dep/data/in/protein_quant_current_normalized.csv.gz", header=T)
Table_S1_Sample_Information <- read.csv("/media/chenglab-bp/aquila/dep/data/in/Table_S1_Sample_Information.csv")

# reading in dep data
Model <- read.csv("/media/chenglab-bp/aquila/dep/data/in/Model.csv")
CRISPRGeneDependency <- read.csv("/media/chenglab-bp/aquila/dep/data/in/CRISPRGeneDependency.csv")



colnames(protein_quant_current_normalized.csv)
CRISPRGeneDependency$ModelID
Table_S1_Sample_Information

