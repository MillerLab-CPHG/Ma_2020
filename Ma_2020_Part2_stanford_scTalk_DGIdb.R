# devtools::install_github("VCCRI/scTalk", build = TRUE, build_vignettes = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
#### libraries and db building ----
library(scTalk)
library(STRINGdb)
library(Seurat)

# only keep stanford seurat object
rm(list=setdiff(ls(), "stanford"))
string_db <- STRINGdb$new( version="11", species=9606, score_threshold=200, input_directory="")
## to call up the vignette
# browseVignettes("scTalk")

#### set the active identity to singleR assignments ----
Idents(stanford) <- stanford@meta.data[["SingleR.calls"]]

#### Overall Approach ----
# edge weights >> network paths >> permutation testing >> visualisation
# scTALK works well on seurat obj with cell identities

#### STEP 1: edge weights ----
## Define clusters to use in the analysis
# here we use the stanford TCF21 dataset
# i changed it here so it calls for the singleR annotated instead
populations.use <- names(stanford@meta.data[["SingleR.calls"]])

## Define a label for output files
file.lab <- "TIP"

GenerateEdgeWeights(seurat.object = stanford,
                    file.label = file.lab, # output file
                    species = "human",
                    populations.use = populations.use,
                    string.ver = "11.0",
                    verbose = TRUE) # this is a bug, need to run this line

edge.table <- read.csv(paste0(file.lab, "_all_ligand_receptor_network_edges.csv"))
edge.table <- edge.table[order(edge.table$weight, decreasing = TRUE), ]
head(edge.table)

#### STEP 2: network paths ----
GenerateNetworkPaths(file.label = file.lab,
                     min.weight = 1.5,
                     ncores = 4)
paths.table <- read.csv(paste0(file.lab, "_network_paths_weight1.5.csv"), row.names=1, stringsAsFactors = FALSE)
paths.table <- paths.table[order(paths.table$Weight, decreasing = TRUE), ]
head(paths.table)

#### STEP 3: select for significant cell-cell interactions ----
# permutation test on the number of cell-cell connections that pass the minimum weight threshold
results.table <- EvaluateConnections(file.label = file.lab, 
                                     ncores = 4,
                                     return.results = TRUE)
head(results.table)

#### STEP 4: plots ----

# set color scheme
manual_color_list <-
  c("rosybrown2", # this is the same code for 'manual_color_list'
    "palevioletred1",
    "cadetblue1",
    "lemonchiffon3",
    "darkseagreen",
    "skyblue3",
    "cadetblue3",
    "lemonchiffon4",
    "darkseagreen1",
    "darkseagreen2",
    "rosybrown3",
    "thistle2",
    "salmon1",
    "palevioletred3",
    "palevioletred4",
    "lightsteelblue3",
    "cadetblue2",
    "thistle3"
  )

p.labels <- sort(levels(stanford@active.ident))

names(manual_color_list) <- p.labels

## Generate a circle plot

cell.network.file = paste0("Permutation_tests_", file.lab, "_network.csv")

pdf("Figure_images/scTalk_circoplot.pdf", width=8, height=6)
circle <- CellCirclePlot(input.file = cell.network.file, 
               adj_pval_thresh = 0.01,
               col.set = manual_color_list)
dev.off()

## Generate an inbound vs outbound plot
pdf("Figure_images/scTalk_inboundoutboundplots.pdf", width=8, height=6)
InboundOutboundPlot(input.file = cell.network.file, 
                    p.adj.thresh = 0.01, 
                    lab.weight.thresh = 0.1,
                    col.set = manual_color_list) +
  ylim(c(0,2000)) +
  xlim(c(0,2500))
dev.off()

## Plot top ligands for SMC population
path.file <- paste0(file.lab, "_network_paths_weight1.5.csv")
cell.type <- "SMC"
PlotTopLigands(input.file = path.file, 
               col.use = "#0099cc",
               cell.identity = cell.type)
ggsave("Figure_images/scTALK_SMCoutgoingligand.pdf", width = 6, height = 5, dpi = 300)


#### STEP 5: tree plots ----
## First, select ligands upregulated in source population
source.marker.table <- Seurat::FindMarkers(stanford, ident.1 = "SMC", 
                                           logfc.threshold = 1.5, 
                                           only.pos = TRUE)
source.marker.genes <- rownames(source.marker.table)
path.file <- paste0(file.lab, "_network_paths_weight1.5.csv")
edge.score.file <- paste0(file.lab, "_all_ligand_receptor_network_edges.csv")

NetworkTreePlot(path.file = path.file, 
                edge.score.file = edge.score.file,
                source.population = "SMC", # change start source as needed
                target.populations = c("FB", "Osteoblasts", "CH"), # change cell target as needed
                source.marker.genes = source.marker.genes,
                population.cols = c("cadetblue2", "cadetblue3", "salmon1", "cadetblue1")) # add more as needed
ggsave("Figure_images/scTALK_treeplot.pdf", width = 4, height = 3, dpi = 300)


#### DGIdb druggability ----
# STEP 1: generate Gene list ----
# BiocManager::install("rDGIdb")
genes <- c(
  "FBLN1", "APOD", "DCN", "ITGB1", "LEPR", "EGFR", "TCF21")
result <- queryDGIdb(genes,
                     sourceDatabases = c("DrugBank", "FDA", "ChemblInteractions"),
                     geneCategories = "CLINICALLY ACTIONABLE")

# to look at the results
resultSummary(result)
detailedResults(result)

# STEP 2: basic visualization
# this shows which db has interactions
plotInteractionsBySource(result, main = "Number of interactions by source")
resourceVersions() # this command shows the db version used
