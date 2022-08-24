#### Import data

DGE_folder <- "../../data/scRNAseq/"
mat <- readMM(paste0(DGE_folder, "DGE.mtx"))
cell_meta <- read.delim(paste0(DGE_folder, "cell_metadata.csv"),
                        stringsAsFactor = FALSE, sep = ",")
genes <- read.delim(paste0(DGE_folder, "all_genes.csv"),
                    stringsAsFactor = FALSE, sep = ",")

cell_meta$bc_wells <- make.unique(cell_meta$bc_wells, sep = "_dup")

rownames(cell_meta) <- cell_meta$bc_index
genes$gene_name <- make.unique(genes$gene_name, sep = "_dup")

# Set column and rownames to expression matrix
colnames(mat) <- genes$gene_name
rownames(mat) <- rownames(cell_meta)
mat_t <- t(mat)

# Remove empty rownames, if they exist
mat_t <- mat_t[(rownames(mat_t) != ""), ]
df_raw <- CreateSeuratObject(mat_t, min_cells = 2, meta.data = cell_meta)

# Setting initial cell class to a single type (instead of plate well numbers), this will change after clustering. 
df_raw@meta.data$orig.ident <- factor(rep("df", nrow(df_raw@meta.data)))
Idents(df_raw) <- df_raw@meta.data$orig.ident

# Add % mitochondrial reads, condition and genotype metadata
df_raw[["condition"]] <- gsub("(w*)_\\d*_E\\d", "\\1", df_raw@meta.data$sample)
df_raw[["genotype"]] <- gsub(".*_(E\\d)", "\\1", df_raw@meta.data$sample)
df_raw[["percent.mt"]] <- PercentageFeatureSet(df_raw, pattern = "^mt-")
