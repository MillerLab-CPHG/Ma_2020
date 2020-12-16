# download the source package of scCATCH-2.1.tar.gz and install it
# ensure the right directory for scCATCH-2.1.tar.gz

library(scCATCH) # devtools::install_github('ZJUFanLab/scCATCH')


stanford <- readRDS(file = "final_stanford_labeled.rds")
clu_markers <- findmarkergenes(
  stanford,
  species = "Human",
  cluster = 'All',
  match_CellMatch = FALSE, # set T for large dataset
  cancer = NULL,
  tissue = NULL,
  cell_min_pct = 0.025,
  logfc = 0.025,
  pvalue = 0.05
)

clu_ann <- scCATCH(clu_markers$clu_markers,
                              species = "Human",
                              cancer = NULL,
                              tissue = "Blood vessel")
write.csv(clu_ann, file = "scCATCH_vs_singleR_bloodvessel.csv")


clu_ann <- scCATCH(clu_markers$clu_markers,
                   species = "Human",
                   cancer = NULL,
                   tissue = "Blood")
write.csv(clu_ann, file = "scCATCH_vs_singleR_blood.csv")

clu_ann <- scCATCH(clu_markers$clu_markers,
                   species = "Human",
                   cancer = NULL,
                   tissue = "Serum")
write.csv(clu_ann, file = "scCATCH_vs_singleR_serum.csv")

clu_ann <- scCATCH(clu_markers$clu_markers,
                   species = "Human",
                   cancer = NULL,
                   tissue = "Heart")
write.csv(clu_ann, file = "scCATCH_vs_singleR_heart.csv")

clu_ann <- scCATCH(clu_markers$clu_markers,
                   species = "Human",
                   cancer = NULL,
                   tissue = "Myocardium")
write.csv(clu_ann, file = "scCATCH_vs_singleR_myocardium.csv")