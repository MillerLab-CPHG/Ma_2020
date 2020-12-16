#### Library and data loading ----
library(Seurat)
library(patchwork)
library(readr)
library(scCATCH)
library(SingleR)
library(tidyverse)
library(monocle3)
library(SeuratData)
library(magrittr)
library(colourpicker)
library(rDGIdb)
library(dyno)

#### step1: loading the data ####
# in Seurat the columns represent cells, whereas in dyno the rows represent cells. 
# you need to transpose your input matrix:
# object_counts <- Matrix::t(as(as.matrix(object@assays$RNA@counts), 'sparseMatrix'))
# object_expression <- Matrix::t(as(as.matrix(object@assays$RNA@data), 'sparseMatrix'))
# object_dyn<- wrap_expression(
#   counts = object_counts, 
#   expression = object_expression
# )

stanford <- readRDS(file = "final_stanford_labeled.rds")


#### step2: selecting interested cells for a subset ####
## uncomment line31 to enable subsetting
stanford.subset <- stanford
stanford.subset <- subset(stanford, idents = c("SMC", "FB", "CH"))

#### step3: construct dyno obj ####
object_counts <- Matrix::t(as(as.matrix(stanford.subset@assays$RNA@counts), 'sparseMatrix'))
object_expression <- Matrix::t(as(as.matrix(stanford.subset@assays$RNA@data), 'sparseMatrix'))
object_cellinfo <- stanford.subset@meta.data[["SingleR.labels"]]


stanford.dyno <- wrap_expression(
  counts = object_counts,
  expression = object_expression)

#### step4: choosing the right methods ####
# this will call a shiny app #
# dynguidelines::guidelines_shiny()
#### slingshot: construct the model ####
model <- infer_trajectory(stanford.dyno, "slingshot")

#### slingshot: project the model ###
# add dim reduction
model <- model %>% 
  add_dimred(dimred = stanford@reductions[["umap"]],
             expression_source = stanford.dyno$expression)

pdf("Figure_images/dyno_slingshot_full.pdf", width=7, height=6)
plot_dimred(
  model, 
  expression_source = stanford.dyno$expression, 
  grouping = object_cellinfo # basically stanford@meta.data[["SingleR.labels"]]
)
dev.off()

#### slingshot: show a gene expression
plot_dimred(
  model, 
  expression_source = stanford.dyno$expression, 
  feature_oi = "FN1"
)

## Coloring by grouping



#### scorpius: construct the model ####
model <- infer_trajectory(stanford.dyno, "scorpius")

#### scorpius: project the model ###
# add dim reduction
model <- model %>% 
  add_dimred(dyndimred::dimred_umap,
             expression_source = stanford.dyno$expression)
plot_dimred(
  model, 
  expression_source = stanford.dyno$expression, 
  grouping = object_cellinfo # basically stanford@meta.data[["SingleR.labels"]]
)

#### scorpius: show a gene expression
plot_dimred(
  model, 
  expression_source = stanford.dyno$expression, 
  feature_oi = "FN1",
  alpha = 0.6
)

plot_dimred(
  model, 
  expression_source = stanford.dyno$expression, 
  feature_oi = "MYH11",
  alpha = 0.6
)

#### PAGA: construct the model ####
model <- infer_trajectory(stanford.dyno, "paga")

#### PAGA: project the model ###
# add dim reduction
model %>% add_dimred(dimred = stanford@reductions[["umap"]])

model <- model %>% 
  add_dimred(dyndimred::dimred_umap,
             expression_source = stanford.dyno$expression)

plot_dimred(
  model, 
  expression_source = stanford.dyno$expression, 
  grouping = object_cellinfo # basically stanford@meta.data[["SingleR.labels"]]
)

#### scorpius: show a gene expression
plot_dimred(
  model, 
  expression_source = stanford.dyno$expression, 
  feature_oi = "FN1",
  alpha = 0.6
)

plot_dimred(
  model, 
  expression_source = stanford.dyno$expression, 
  feature_oi = "MYH11",
  alpha = 0.6
)



