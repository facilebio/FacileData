### For testing
### Delete before OSS

library(FacileData)
library(ExpressionPlot)

ngs114 = ep.project.from.id("PRJ0013166")
ngs171 = ep.project.from.id("PRJ0013155")

ds = list(
    NGS114 = ep.SummarizedExperiment(ngs114, attach.annot = TRUE),
    NGS171 = ep.SummarizedExperiment(ngs171, attach.annot = TRUE)
)

ds = lapply(ds,
            function(x) {
                rownames(x) = paste0("GeneID:",rownames(x))
                f = colnames(mcols(x)) = c("aliases", "effective_length", "feature_type", "name", "meta")
                mcols(x)$source = "IGIS"
                mcols(x)$feature_type = "entrez"
                mcols(x)$feature_id = rownames(x)
                colData(x) = colData(x)[,c("SAMPLE_ID","CLID","AGE","TISSUE_METACLASS_ONCOLOGY")]
                colnames(colData(x)) = c("samid","clid","age","tissue")
                metadata(colData(x)) = list(
                    samid = list(label = "Sample Hub ID",
                                 description = "This is a sample ID",
                                 type = "general"
                                 ),
                    clid = list(label = "gCell Cell Line ID",
                                description = "This is another ID",
                                type = "general"
                                ),
                    age = list(label = "age in years",
                               description = "The age of the patient",
                               type = "general"
                               ),
                    tissue = list(label = "Tissue Group",
                                  description = "Rollup of tissue type to defined vocab",
                                  type = "clinical"
                                  )
                )
                x
            })

metadata(ds[[1]]) = list(url = "http://gene.com", description = "This is NGS114")
metadata(ds[[2]]) = list(url = "http://gene.com", description = "This is NGS171")

DIR = "/gne/research/workspace/phaverty/FDS"
unlink(DIR, recursive = TRUE)
gcell_fds = as.FacileDataSet(ds,
                       path = DIR,
                       assay_name = "rnaseq",
                       assay_type = "rnaseq",
                       source_assay = "counts",
                       dataset_name = "gCell",
                       organism = "Homo sapiens"
                       )


## Try it out
library(dplyr)
library(reshape2)

samples <- sample_info_tbl(gcell_fds) %>% collect
samples <- fetch_sample_covariates(gcell_fds, samples, "tissue") %>%
    filter(variable == "tissue" & value == "Breast") %>%
    collect
genes <- c(ERBB2 = 'GeneID:2064', GRB7 = 'GeneID:2886')
exprs <- gcell_fds %>% fetch_assay_data(genes, foo, 'rnaseq', normalized = TRUE)
ew <- exprs %>% dcast(dataset + sample_id ~ feature_name, value.var = 'value')
stopifnot(cor(ew[,"ERBB2"], ew[,"GRB7"]) > 0.6)
