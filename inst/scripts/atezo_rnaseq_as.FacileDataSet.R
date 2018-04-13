### For testing
### Delete before OSS

library(FacileData)
library(S4Vectors)
library(tidyverse)

## load the list of ExpressionSets previously used
atezo.se.list <- readRDS("~/Documents/Analysis/FE/facileatezo-eSets-list_withIMvigor211_irAEcm_20180407_priorConOnly.rds")

## needs to satisfy:
# req.cols <- c(
#   feature_type="character",
#   feature_id="character",
#   name="character",
#   meta="character",
#   effective_length="numeric",
#   source="character"
# )
## currently have:
names(fData(atezo.se.list$es[[1]]))

## update some pData and fData columns for consistency
for (i in seq_along(atezo.se.list$es)) {
  pData(atezo.se.list$es[[i]])$group <- as.numeric(as.character(pData(atezo.se.list$es[[i]])$group))
  fData(atezo.se.list$es[[i]])$feature_type <- "entrez"
  fData(atezo.se.list$es[[i]])$name <- fData(atezo.se.list$es[[i]])$symbol ## ?
  fData(atezo.se.list$es[[i]])$meta <- "" ## ?
  fData(atezo.se.list$es[[i]])$effective_length <- 0L ## ?
}

## convert to SummarizedExperiments
atezo.se.list <- atezo.se.list %>% 
  mutate(se = map(es, ~as(., "SummarizedExperiment")))

## extract SummarizedExperiments into a named list
ds <- atezo.se.list$se
names(ds) <- atezo.se.list$dataset

## add metadata
metadata(ds[[1]]) = list(url = "http://gene.com", description = "This is birch")
metadata(ds[[2]]) = list(url = "http://gene.com", description = "This is fir")
metadata(ds[[3]]) = list(url = "http://gene.com", description = "This is imvigor210")
metadata(ds[[4]]) = list(url = "http://gene.com", description = "This is pcd")
metadata(ds[[5]]) = list(url = "http://gene.com", description = "This is poplar")
metadata(ds[[6]]) = list(url = "http://gene.com", description = "This is imvigor211")

## create new FacileAtezoDataSet .h5, .sqlite, .yaml
DIR = "~/Documents/Analysis/FE/testAsFacileDataSet/"
unlink(DIR, recursive = TRUE)
newFADS = as.FacileDataSet(ds,
                           path = DIR,
                           assay_name = "rnaseq",
                           assay_type = "rnaseq",
                           source_assay = "exprs",
                           dataset_name = "atezo",
                           organism = "Homo sapiens"
)

#### TESTING #### 

library(reshape2)
samples <- sample_info_tbl(newFADS) %>%
  filter(dataset == 'imvigor210') %>%
  collect
genes <- c(TIGIT='201633', CD274='29126')
exprs <- newFADS %>%
  fetch_assay_data(genes, samples, 'rnaseq', normalized=TRUE)
ew <- exprs %>%
  dcast(dataset + sample_id ~ feature_name, value.var='value') %>%
  with_sample_covariates(c('subtype_tcga_bladder'), .fds=newFADS)

## FacileExplorer pre refactor has this @ 0.64
cor(ew$CD274, ew$TIGIT, method='spearman') ## now 0.6705419 (!?)

## Correlation by subtype
ocor <- tibble(
  subtype_tcga_bladder=c("I", "II", 'III', "IV", "na"),
  n=c(99, 97, 60, 62, 36),
  spearman=c(0.54, 0.71, 0.54, 0.53, 0.57))
ncor <- ew %>%
  group_by(subtype_tcga_bladder) %>%
  summarize(n=n(),
            spearman=round(cor(CD274, TIGIT, method='spearman'), 2))
all.equal(ocor$spearman, ncor$spearman) ## YES (almost)!

library(ggplot2)
theme_set(theme_bw())
ggplot(ew, aes(CD274, TIGIT)) +
  geom_point() ## match with original FacileExplorer
ggplot(ew, aes(CD274, TIGIT)) +
  geom_point() +
  facet_wrap(~ subtype_tcga_bladder) ## match with original facile explorer
