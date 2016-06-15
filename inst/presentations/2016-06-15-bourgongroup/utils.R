plot.cor <- function(dat, ind, sub, cor.method=c('spearman', 'pearson'),
                     method='ellipse', type='upper', order='hclust', addrect=2,
                     ...) {
  cor.method <- match.arg(cor.method)
  cor.cols <- c(
    "#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
    "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061") %>%
    rev %>%
    colorRampPalette
  if (!missing(ind) && !missing(sub)) {
    X <- filter(dat, indication == ind & subtype == sub) %>%
      mcast(sample ~ name, value.var='score')
  } else {
    X <- mcast(dat, sample ~ name, value.var='score')
  }
  XC <- cor(X, method=cor.method)
  colnames(XC) <- NULL
  corrplot(XC, method=method, type=type, order=order, addrect=addrect,
           col=cor.cols(50), ...)
}
