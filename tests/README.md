Tests and Coverage
================
01 June, 2018 23:41:44

-   [Coverage](#coverage)
-   [Unit Tests](#unit-tests)

This output is created by [covrpage](https://github.com/yonicd/covrpage).

Coverage
--------

Coverage summary is created using the [covr](https://github.com/r-lib/covr) package.

| Object                                                      | Coverage (%) |
|:------------------------------------------------------------|:------------:|
| FacileData                                                  |     41.64    |
| [R/build-db.R](../R/build-db.R)                             |     0.00     |
| [R/features.R](../R/features.R)                             |     0.00     |
| [R/NSE-filter-features.R](../R/NSE-filter-features.R)       |     0.00     |
| [R/NSE-filter-samples.R](../R/NSE-filter-samples.R)         |     0.00     |
| [R/sql.R](../R/sql.R)                                       |     0.00     |
| [R/construction.R](../R/construction.R)                     |     4.44     |
| [R/as.FacileDataSet.R](../R/as.FacileDataSet.R)             |     5.32     |
| [R/samples.R](../R/samples.R)                               |     27.47    |
| [R/utilities.R](../R/utilities.R)                           |     44.44    |
| [R/FacileDataSet.R](../R/FacileDataSet.R)                   |     50.46    |
| [R/db-and-table-functions.R](../R/db-and-table-functions.R) |     56.91    |
| [R/as.BiocAssayContainers.R](../R/as.BiocAssayContainers.R) |     57.55    |
| [R/entity-attribute-value.R](../R/entity-attribute-value.R) |     60.17    |
| [R/validation.R](../R/validation.R)                         |     60.38    |
| [R/assay-data.R](../R/assay-data.R)                         |     65.30    |
| [R/sample-covariates.R](../R/sample-covariates.R)           |     75.83    |
| [R/test-helpers.R](../R/test-helpers.R)                     |     92.86    |
| [R/sample-info.R](../R/sample-info.R)                       |     93.33    |

<br>

Unit Tests
----------

Unit Test summary is created using the [testthat](https://github.com/r-lib/testthat) package.

| file                                                                    |    n|   time|  error|  failed|  skipped|  warning|
|:------------------------------------------------------------------------|----:|------:|------:|-------:|--------:|--------:|
| [test-as.FacileDataSet.R](testthat/test-as.FacileDataSet.R)             |    0|  0.084|      1|       0|        0|        0|
| [test-assay-data.R](testthat/test-assay-data.R)                         |    4|  1.910|      0|       0|        0|        0|
| [test-bioc-assay-containers.R](testthat/test-bioc-assay-containers.R)   |   14|  0.706|      0|       0|        0|        0|
| [test-entity-attribute-value.R](testthat/test-entity-attribute-value.R) |   44|  0.077|      1|       0|        0|        0|
| [test-FacileDataSet.R](testthat/test-FacileDataSet.R)                   |    3|  0.013|      0|       0|        0|        0|
| [test-sample-covariates.R](testthat/test-sample-covariates.R)           |   24|  1.043|      0|       0|        0|        0|
| [test-samples.R](testthat/test-samples.R)                               |    1|  0.034|      0|       0|        0|        0|
| [test-utilities.R](testthat/test-utilities.R)                           |    2|  0.016|      0|       0|        0|        0|

| file                                                                    | test                                                                | context                                          | status |    n|   time|
|:------------------------------------------------------------------------|:--------------------------------------------------------------------|:-------------------------------------------------|:-------|----:|------:|
| [test-as.FacileDataSet.R](testthat/test-as.FacileDataSet.R)             | We can get pdata metadata                                           | as.FacileDataSet                                 | ERROR  |    0|  0.084|
| [test-assay-data.R](testthat/test-assay-data.R)                         | fetch\_assay\_data limits samples correctly                         | Fetching assay level data                        | PASS   |    2|  0.834|
| [test-assay-data.R](testthat/test-assay-data.R)                         | spreading data works with\_assay\_data                              | Fetching assay level data                        | PASS   |    1|  0.597|
| [test-assay-data.R](testthat/test-assay-data.R)                         | fetch\_assay\_data(..., aggregate.by='ewm') provides scores         | Fetching assay level data                        | PASS   |    1|  0.479|
| [test-bioc-assay-containers.R](testthat/test-bioc-assay-containers.R)   | fetch\_assay\_data results converted to DGEList                     | Testing conversion to Bioc Expression Containers | PASS   |   14|  0.706|
| [test-entity-attribute-value.R](testthat/test-entity-attribute-value.R) | pData -&gt; meta.yaml covariate encoding works (simple & compound)  | Entity-Attribute-Value conversions               | PASS   |   44|  0.076|
| [test-entity-attribute-value.R](testthat/test-entity-attribute-value.R) | basic encoding and decoding of EAV columns works                    | Entity-Attribute-Value conversions               | ERROR  |    0|  0.001|
| [test-FacileDataSet.R](testthat/test-FacileDataSet.R)                   | Fetching various database tables from FacileDataSet                 | Basic FacileDataSet functions                    | PASS   |    3|  0.013|
| [test-sample-covariates.R](testthat/test-sample-covariates.R)           | fetch\_sample\_covariates retrieves all covariates if not specified | Sample Covariates                                | PASS   |    1|  0.069|
| [test-sample-covariates.R](testthat/test-sample-covariates.R)           | fetch\_sample\_covariates::samples arg limits samples correctly     | Sample Covariates                                | PASS   |    1|  0.046|
| [test-sample-covariates.R](testthat/test-sample-covariates.R)           | cast\_covariate converts simple variables to correct type           | Sample Covariates                                | PASS   |    6|  0.092|
| [test-sample-covariates.R](testthat/test-sample-covariates.R)           | cast\_covariate converts right\_censored data correctly             | Sample Covariates                                | PASS   |    8|  0.057|
| [test-sample-covariates.R](testthat/test-sample-covariates.R)           | spread\_covariates casts simple covariates to correct class         | Sample Covariates                                | PASS   |    4|  0.086|
| [test-sample-covariates.R](testthat/test-sample-covariates.R)           | spread\_covariates works with both simple and complex types         | Sample Covariates                                | PASS   |    2|  0.110|
| [test-sample-covariates.R](testthat/test-sample-covariates.R)           | with\_sample\_covariates returns long input with wide covariates    | Sample Covariates                                | PASS   |    1|  0.492|
| [test-sample-covariates.R](testthat/test-sample-covariates.R)           | successive with\_sample\_covariate calls build correct frame        | Sample Covariates                                | PASS   |    1|  0.091|
| [test-samples.R](testthat/test-samples.R)                               | fetch samples returns valid sample descriptor absent assay spec     | Retrieving arbitrary samples                     | PASS   |    1|  0.034|
| [test-utilities.R](testthat/test-utilities.R)                           | We can bind\_rows data.frames with Surv columns                     |                                                  | PASS   |    2|  0.016|
