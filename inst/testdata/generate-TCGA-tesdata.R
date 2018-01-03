# Using FacileTCGADataSet to create some testdata
library(FacileTCGADataSet)
library(magrittr)
tcga <- FacileTCGADataSet()

## let's get some 20 samples
set.seed(0xBEEF)
bsamples.all <- tcga %>%
  filter_samples(indication %in% c("BLCA", "BRCA")) %>%
  with_sample_covariates %>%
  filter(sample_type != 'tumor_metastatic')

bsamples <- bsamples.all %>%
  group_by(dataset, sample_type) %>%
  sample_n(5) %>%
  ungroup %>%
  set_fds(tcga)

# pData for testing entity-attribute-value encodings
# Create a covariate pData object with non-default factor levels to test
scovs <- bsamples %>%
  mutate(stage = factor(stage),
         sex = factor(sex, c("male", "female"))) %>%
  select(dataset, sample_id, stage, sex, age, sample_type,
         subtype_molecular_bladder, subtype_receptor_breast,
         tte_OS, event_OS)

# Let's fill the categorical variables with all levels, even though our sampling
# can't possibly do that.
scovs %<>% mutate(stage = sub("[ab]$", "", stage))
scovs %<>% mutate(stage = factor(stage, paste("stage", c("i", "ii", "iii", "iv"))))
is.blca.tumor <- with(bsamples, dataset == "BLCA" & sample_type == "tumor")
is.brca.tumor <- with(bsamples, dataset == "BRCA" & sample_type == "tumor")
blca.sub.lvls <- c("luminal", "basal")
brca.sub.lvls <- c("ER+/PR+", "Her2+", "TNBC")
scovs %<>%
  mutate(subtype_molecular_bladder = ifelse(
    is.blca.tumor,
    sample(blca.sub.lvls, sum(is.blca.tumor), replace = TRUE), NA))
scovs %<>%
  mutate(subtype_receptor_breast = ifelse(
    is.brca.tumor,
    sample(brca.sub.lvls, sum(is.brca.tumor), replace = TRUE), NA))
scovs %<>%
  mutate(
    # keep subtype_molecular_bladder as just a character
    # subtype_molecular_bladder = factor(subtype_molecular_bladder, blca.sub.lvls),
    subtype_receptor_breast = factor(subtype_receptor_breast, brca.sub.lvls))
saveRDS(scovs, "test-sample-covariates.rds")

