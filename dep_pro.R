library(stringr)


### ingest data 
# reading in data proteomics
protein_quant_current_normalized.csv <- read.csv("./data/in/protein_quant_current_normalized.csv.gz", header=T)
Table_S1_Sample_Information <- read.csv("./data/in/Table_S1_Sample_Information.csv")

# reading in dep data
Model <- read.csv("./data/in/Model.csv")
CRISPRGeneDependency <- read.csv("./data/in/CRISPRGeneDependency.csv.gz")




### purging unequiv data
## purging unequiv cell lines
# subtracting the differences and chekcing assumption
Model <- Model[-which(!Model$ModelID %in% CRISPRGeneDependency$ModelID), ]
identical(CRISPRGeneDependency$ModelID, Model$ModelID)

# replacing the model id with cell line names
CRISPRGeneDependency$ModelID <- Model$CCLEName



## handling duplicates
# note that duplicates exist in proteomics and in crispr screen
# there are dups in crispr but we need not handle them since all have a unique tag
sum(duplicated(protein_quant_current_normalized.csv$Gene_Symbol))
sum(duplicated(sapply(str_split(colnames(CRISPRGeneDependency), '\\.'), function(x) x[1])))

# simple solution is to paste protein ID to the gene ID 
hold_dups_pro <- protein_quant_current_normalized.csv[which(duplicated(protein_quant_current_normalized.csv$Gene_Symbol)), ]
hold_dups_pro$Gene_Symbol <- paste(hold_dups_pro$Gene_Symbol, hold_dups_pro$Uniprot_Acc, sep='_')

# subtracting dups then mashing hold dups and proteomic data frame together
protein_quant_current_normalized.csv <- protein_quant_current_normalized.csv[-which(duplicated(protein_quant_current_normalized.csv$Gene_Symbol)), ]
protein_quant_current_normalized.csv <- rbind(protein_quant_current_normalized.csv,hold_dups_pro)

# checking duplicated data
sum(duplicated(protein_quant_current_normalized.csv$Gene_Symbol))

# making rownames in protein gene symbol to split frame
rownames(protein_quant_current_normalized.csv) <- protein_quant_current_normalized.csv$Gene_Symbol

# splitting proteomics data into meta and expression values
protein_meta <- protein_quant_current_normalized.csv[grep('^TenPx|Protein|Gene_Sym|Group_I|Uniprot|Description', colnames(protein_quant_current_normalized.csv))]
protein_quant_current_normalized.csv <- protein_quant_current_normalized.csv[-grep('^TenPx|Protein|Gene_Sym|Group_I|Uniprot|Description', colnames(protein_quant_current_normalized.csv))]





### purging unequiv cell lines
# making samples equiv
CRISPRGeneDependency <- CRISPRGeneDependency[which(CRISPRGeneDependency$ModelID %in% sapply(str_split(colnames(protein_quant_current_normalized.csv), '_'), function(x) paste(x[1:length(x)-1], collapse = '_'))),]
protein_quant_current_normalized.csv <- protein_quant_current_normalized.csv[which(sapply(str_split(colnames(protein_quant_current_normalized.csv), '_'), function(x) paste(x[1:length(x)-1], collapse = '_')) %in% CRISPRGeneDependency$ModelID )]

# note there are duplicates of cell lines
hold_dups <- protein_quant_current_normalized.csv[which(duplicated(sapply(str_split(colnames(protein_quant_current_normalized.csv), '_'), function(x) paste(x[1:length(x)-1], collapse = '_'))))]
hold_dups_t <- protein_quant_current_normalized.csv[which(sapply(str_split(colnames(protein_quant_current_normalized.csv), '_'), function(x) paste(x[1:length(x)-1], collapse = '_'))  %in% sapply(str_split(colnames(hold_dups), '_'), function(x) paste(x[1:length(x)-1], collapse = '_')) )]

# checking correlation between two cell lines: its pretty low
cor(hold_dups_t$CAL120_BREAST_TenPx02,hold_dups_t$CAL120_BREAST_TenPx28, use ='pairwise.complete.obs')
plot(hold_dups_t$CAL120_BREAST_TenPx02,hold_dups_t$CAL120_BREAST_TenPx28)

cor(hold_dups_t$HCT15_LARGE_INTESTINE_TenPx18,hold_dups_t$HCT15_LARGE_INTESTINE_TenPx30, use ='pairwise.complete.obs')
plot(hold_dups_t$HCT15_LARGE_INTESTINE_TenPx18,hold_dups_t$HCT15_LARGE_INTESTINE_TenPx30)

# going to just dedup since these are so different
protein_quant_current_normalized.csv <- protein_quant_current_normalized.csv[-which(duplicated(sapply(str_split(colnames(protein_quant_current_normalized.csv), '_'), function(x) paste(x[1:length(x)-1], collapse = '_'))))]

# ordering data for ident check 
CRISPRGeneDependency <- CRISPRGeneDependency[order(CRISPRGeneDependency$ModelID), ]
protein_quant_current_normalized.csv <- protein_quant_current_normalized.csv[order(colnames(protein_quant_current_normalized.csv))]

identical(sapply(str_split(colnames(protein_quant_current_normalized.csv), '_'), function(x) paste(x[1:length(x)-1], collapse = '_')) , CRISPRGeneDependency$ModelID)

# cleaning
rm(hold_dups, hold_dups_crispr, hold_dups_pro, hold_dups_t)
