## =======================================================================================
## Reorder and rename clusters in averaged Seurat object
## =======================================================================================

rename_clusters <- c(`1` = "Airway epithelial A", `2` = "Airway epithelial B", 
                     `3` = "Granulocytes", `4` = "Myofibroblasts", 
                     `5` = "Lipofibroblasts", `6` = "AT1", 
                     `7` = "misc_01", `8` = "Ciliated cells", 
                     `9` = "Mesothelial cells", `10` = "Pericytes", 
                     `11` = "misc_02", `12` = "Col14a1pos fibroblasts", 
                     `13` = "misc_06" , `14` = "Vcam1pos ECs A", 
                     `15` = "NK cells", `16` = "B cells", 
                     `17` = "T cells A", `18` = "T cells B", 
                     `19` = "T cells C", `20` = "Neuronal", 
                     `21` = "Capillary ECs", `22` = "Vascular ECs A", 
                     `23` = "Other ECs", `24` = "Vcam1pos ECs B", 
                     `25` = "Vascular ECs B", `26` = "misc_03", 
                     `27` = "AT2", `28` = "Alveolar mø A", 
                     `29` = "Alveolar mø B", `30` = "Alveolar mø prolif", 
                     `31` = "DCs", `32` = "misc_04", `33` = "Monocytes A", 
                     `34` = "Monocytes B", `35` = "Interstitial mø")

df <- RenameIdents(object = df, rename_clusters)
df$celltype <- Idents(object = df)