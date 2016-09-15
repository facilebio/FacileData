library(FacileRepo)
library(RSQLite)
library(reshape2)

fn <- '/Users/lianogls/workspace/data/facile/test/TcgaDb-test.sqlite'
to <- sub('.sqlite', '-wide.sqlite', fn)

file.copy(fn, to)
db <- dbConnect(RSQLite::SQLite(), to)

exprs <- dbGetQuery(db, 'select * from expression')
E <- dcast(exprs, sample_id ~ feature_id, value.var='count')
colnames(E)[-1] <- paste0('fid_', colnames(E)[-1])

dbGetQuery(db, 'DROP INDEX expression_gene')
dbGetQuery(db, 'DROP INDEX expression_sample_gene')
dbGetQuery(db, 'DROP TABLE expression')

col.types <- c("TEXT", rep("INTEGER", ncol(E) - 1))
names(col.types) <- colnames(E)

dbWriteTable(db, 'expression', E,
             row.names=FALSE,
             field.types=col.types)
dbSendQuery(db, "CREATE UNIQUE INDEX expression_sample ON expression (sample_id)")
dbDisconnect(db)

src <- src_sqlite(to)
exprs <- tbl(src, 'expression')

x = copy_to(src, E, name='expression')


