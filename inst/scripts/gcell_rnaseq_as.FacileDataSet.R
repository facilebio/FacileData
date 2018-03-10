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
                f = colnames(mcols(x)) = c("aliases","effective_length", "feature_type","name","meta")
                mcols(x)$source = "IGIS"
                mcols(x)$feature_type = "entrez"
                mcols(x)$feature_id = rownames(x)
                colData(x) = colData(x)[,c("SAMPLE_ID","CLID","AGE","TISSUE_METACLASS_ONCOLOGY")]
                colnames(colData(x)) = c("samid","clid","age","tissue")
                metadata(colData(x)) = list(samid = "Sample Hub ID", clid = "gCell Cell Line ID", age = "age in years", tissue = "Tissue Group")
                x
            })

metadata(ds[[1]]) = list(url = "http://gene.com", description = "This is NGS114")
metadata(ds[[2]]) = list(url = "http://gene.com", description = "This is NGS171")

DIR = "/gne/research/workspace/phaverty/FDS"
unlink(DIR, recursive = TRUE)
gcell_fds = as.FacileDataSet(ds,
                       path = DIR,
                       assay_name = "Expression",
                       assay_type = "rnaseq",
                       source_assay = "rpkms",
                       dataset_name = "gCell",
                       organism = "Homo sapiens"
                       )
