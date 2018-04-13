-- This schema supports creating a multiassay FacileDataSet.
--
-- The `assay` value is the name of the hdf5 directory in the hdf5 datastore,
-- and the names of the matrices in each `assay` directory are the `dataset`
-- values.
--
-- I don't think that The different assays across the same dataset need to be of
-- the same width (ie. have the same exact samples measured across each assay).
-- Although the same samples on different assays are expected to have the same
-- sample_id.
--
-- Although different assays (fData of varying "heights") can be assigned to
-- the a single dataset. A given dataset only has ONE pData data.frame

-- The `feature_info` table provides meta information about the types of
-- features that one can measure. Each `assay` matrix will measure one
-- `feature_type`, and the details of the features in that `feature_type` are
-- outlined here.
CREATE TABLE feature_info (
  feature_type TEXT , -- "entrez", "ensgid", "enstid", ...
  feature_id TEXT ,   -- "29126", "ENSG00000120217", etc.
  name TEXT ,         -- "CD274" (if this is a gene, this will be 'symbol')
  meta TEXT ,         -- some unforseen information (maybe PTM (phosphorylated(?)))
--  seqnames TEXT,
--  start INTEGER,
--  end INTEGER,
--  strand INTEGER,    -- 1, -1, 0: positive, negative, unstranded
  effective_length INTEGER,
  source TEXT,
  PRIMARY KEY (feature_type, feature_id));
CREATE INDEX feature_info__feature_id ON feature_info (feature_id);
CREATE INDEX feature_info__name ON feature_info (name);
CREATE INDEX feature_info__feature_type ON feature_info (feature_type);

-- The `dataset_info` provides meta information for each dataset in a
-- FacileDataSet
-- CREATE TABLE dataset_info (
--   dataset TEXT,
--   organism TEXT,
--   PRIMARY KEY (dataset));

-- This `sample_info` table feels too slim. I wanted to put hd5_index into here
-- but we have occasions where we don't run the same assays over the same
-- samples, ie. fluidigm and rnaseq assays don't have 1:1 coverage over our
-- sample space, otherwise it would be great if we could put the hd5_index in
-- here.
CREATE TABLE sample_info (
  dataset TEXT,
  sample_id TEXT,
  parent_id TEXT, -- patient_id or something similar; imagine multiple samples from the same "animal"
  PRIMARY KEY (dataset, sample_id));
CREATE INDEX sample_info__dataset_parent_id ON sample_info (dataset, parent_id);

-- an element of the pData
CREATE TABLE sample_covariate (
  dataset TEXT,
  sample_id TEXT,
  variable TEXT,
  value TEXT,
  class TEXT,
  type TEXT,
  date_entered INTEGER,
  PRIMARY KEY (dataset, sample_id, variable));
CREATE INDEX sample_covariate__variable ON sample_covariate (variable);
CREATE INDEX sample_covariate__value ON sample_covariate (value);
CREATE INDEX sample_covariate__class ON sample_covariate (class);

-- The `assay_info` table holds meta information about the *types* of data
-- stored in each of the top level HDF5 directories. There will be as many
-- entries as there are assays run across the FacileDataSet.
-- Each dataset will have a matrix in one or more of the different
-- assay_info::assay directories.
--
-- All of the features in a given assay must be the same. Currently this isn't
-- enforced *within* the schema.
CREATE TABLE assay_info (
  assay TEXT,
  assay_type TEXT,   -- rnaseq, nanostring, luminex, etc
  feature_type TEXT, -- what's quantitated in assay ("entrez", "ensgid", "enstid", ..)
  description TEXT,
  nfeatures INTEGER,
  storage_mode TEXT,
  PRIMARY KEY (assay));
CREATE INDEX assay_info__assay_type ON assay_info (assay_type);
CREATE INDEX assay_info__feature_type ON assay_info (feature_type);

-- The `assay_sample_info` table holds information for each column (sample)
-- in an assay matrix. It tracks which column a sample belongs to in the hdf5
-- array(s), and currently holds meta information specfic to sample-level
-- variables that are required for RNA-seq.
CREATE TABLE assay_sample_info (
  assay TEXT,
  dataset TEXT,        -- (dataset, sample_id) should be foreign keys into sample_info
  sample_id TEXT,
  hdf5_index INTEGER,
  libsize DOUBLE PRECISION,
  normfactor DOUBLE PRECISION,
  PRIMARY KEY (assay, dataset, sample_id));
CREATE INDEX assay_sample_info__dataset_samid ON assay_sample_info (dataset, sample_id);

-- TODO: add assay_sample_covariate to provide assay specific information

-- The `assay_feature_info` table holds sample level information for each
-- feature in the assay. Most important is the mapping from a feature_id to its
-- cognate row in its HDF5 assay matrix.
CREATE TABLE assay_feature_info (
  assay TEXT,         -- foreign key to assay_info::assay column
  feature_id TEXT,    -- compound key into feature_info::feature_type,feature_id (feature_type is pulled from assay_info table)
  hdf5_index INTEGER, -- the corresponding row of the assay matrix
  PRIMARY KEY (assay, feature_id));
CREATE INDEX assay_feature_info__feature_id ON assay_feature_info (feature_id);


CREATE TABLE feature_set (
  collection TEXT,   -- collection,name are largely the primary key
  name TEXT,         --
  feature_type TEXT, -- "entrez" (including this in PK is largely for future proofing)
  feature_id TEXT,
  PRIMARY KEY (collection, name, feature_type, feature_id));

