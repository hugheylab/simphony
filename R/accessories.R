#' Default function for mapping expected counts to dispersion.
#'
#' The function was estimated from circadian RNA-seq data from mouse liver
#' (PRJNA297287), using local regression in DESeq2. In a negative binomial
#' distribution, variance = mean + mean^2 * dispersion.
#'
#' @format A vectorized function with some attributes.
#'
#' @examples
#' means = 2^(6:10)
#' dispersions = defaultDispFunc(means)
#'
#' @seealso `\link{simphony}`
'defaultDispFunc'


#' Calculate expected expression
#'
#' Calculate expected expression for multiple genes at multiple timepoints in
#' multiple conditions.
#'
#' @param geneMetadata `data.table` with columns `gene`, `base`, `rhyFunc`,
#'   `amp`, and `phase`, where every row corresponds to a gen. If `byCondGroup` ==
#'   `TRUE`, then must also have columns `cond` and `group`.
#' @param period Integer for the period of simulated rhythms.
#' @param times Numeric vector of the times (in the same units as `period`) at
#'   which to calculate expected expression for each row in `geneMetadata`.
#' @param sampleMetadata `data.table` with columns `sample`, `cond`, and
#'   `time`. Either `times` or `sampleMetadata` must be provided, and the former
#'   takes precedence.
#' @param byCondGroup Logical for whether to speed up the calculation by
#'   grouping by the columns `cond` and `group`. Primarily for internal use.
#'
#' @return `data.table` derived from `geneMetadata` (but with more rows),
#'   with additional columns `time` and `mu` and possibly others. If sampling
#'   will use the negative binomial family, `mu` corresponds to log2 counts.
#'
#' @examples
#' library('data.table')
#' geneMetadata = data.table(gene = c('gene_1', 'gene_2'), base = 0,
#'                           amp = c(0, 1), phase = 0, rhyFunc = sin)
#' exprDt = getExpectedExpr(geneMetadata, times = 6:17)
#'
#' @seealso `\link{simphony}`, `\link{getSampledExpr}`
#'
#' @export
getExpectedExpr = function(geneMetadata, period = 24,
                           times = NULL, sampleMetadata = NULL,
                           byCondGroup = is.null(times)) {
  if (!is.null(times)) {
    d = data.table(geneMetadata)[rep(1:.N, each = length(times))]
    d[, time := rep(times, times = nrow(geneMetadata))]
  } else if (!is.null(sampleMetadata)) {
    d = merge(data.table(geneMetadata), sampleMetadata, by = 'cond',
              allow.cartesian = TRUE)
  } else {
    stop('Either times or sampleMetadata must not be NULL.')}

  if (byCondGroup) {
    d[, mu := base + amp * rhyFunc[[1]]((time + phase) * 2 * pi / ..period),
      by = c('cond', 'group')]
  } else {
    d[, mu := base + amp * rhyFunc[[1]]((time + phase) * 2 * pi / ..period),
      by = 1:nrow(d)]}
  return(data.table::copy(d))}


#' Sample expression values
#'
#' Sample gene expression values from the given distributions. This function
#' is used internally by `simphony()`, and should not usually need to be
#' called directly.
#'
#' @param exprDt `data.table` of expected expression. If `family` == 'gaussian',
#'   required columns are `gene`, `sample`, `mu`, and `sd`. If `family` ==
#'   'negbinom', required columns are `gene`, `sample`, `mu`, `dispFunc`, `cond`,
#'   and `group`.
#' @param family Character string for the family of distributions from which
#'   to generate the expression values. Must be 'gaussian' or 'negbinom'.
#' @param inplace Logical for whether to modify in-place `exprDt`, adding a
#'   column `expr` containing the expression values.
#'
#' @return Matrix of expression values, where rows correspond to genes and
#'   columns correspond to samples.
#'
#' @examples
#' library('data.table')
#' set.seed(6022)
#' exprDt = data.table(gene = 'gene_1', sample = c('sample_1', 'sample_2'),
#'                     mu = c(0, 5), sd = 1)
#' exprMat = getSampledExpr(exprDt)
#'
#' @seealso `\link{simphony}`, `\link{getExpectedExpr}`
#'
#' @export
getSampledExpr = function(exprDt, family = 'gaussian', inplace = FALSE) {
  if (!inplace) {
    exprDt = data.table(exprDt)}

  if (family == 'gaussian') {
    exprDt[, expr := stats::rnorm(.N, mu, sd)]
  } else if (family == 'negbinom') {
    # dispFunc is identical for genes of the same group in the same condition
    # this is the way I've figured out how to call functions that are columns
    exprDt[, expr := stats::rnbinom(.N, mu = 2^mu, size = 1/dispFunc[[1]](2^mu)),
           by = c('cond', 'group')]
  } else {
    stop("family must be 'gaussian' or 'negbinom'.")}

  data.table::setorderv(exprDt, c('sample', 'gene'))
  genes = unique(exprDt$gene)
  samples = unique(exprDt$sample)
  exprMat = matrix(exprDt$expr, nrow = length(genes),
                   dimnames = list(genes, samples))
  return(exprMat)}


#' Merge expression data, gene metadata, and sample metadata
#'
#' Merge a simulation's expression data, gene metadata, and sample metadata
#' into one `data.table`. This function is useful for making plots using
#' ggplot2.
#'
#' @param simData List with the following elements, such as returned by
#' `simphony()`:
#' \describe{
#'   \item{exprData}{Matrix of expression values, with rownames for genes and
#'   colnames for samples.}
#'   \item{sampleMetadata}{`data.table` with columns `sample` and `cond`.}
#'   \item{geneMetadata}{`data.table` with columns `gene` and `cond`.}
#' }
#' @param genes Character vector of genes for which to get expression data. If
#'   NULL, then all genes.
#'
#' @return `data.table`.
#'
#' @examples
#' library('data.table')
#' exprGroups = data.table(amp = c(0, 1))
#' simData = simphony(exprGroups)
#' mergedSimData = mergeSimData(simData, simData$geneMetadata$gene[1:2])
#'
#' @seealso `\link{simphony}`
#'
#' @export
mergeSimData = function(simData, genes = NULL) {
  if (is.null(genes)) {
    genes = rownames(simData$exprData)}

  d = data.table(simData$exprData, keep.rownames = TRUE)
  setnames(d, 'rn', 'gene')
  d = melt(d, id.vars = 'gene', variable.name = 'sample', value.name = 'expr')

  d = merge(d, simData$sampleMetadata, by = 'sample')
  d = merge(d, simData$geneMetadata, by = c('gene', 'cond'))
  return(d)}


#' Split differential exprGroups
#'
#' Split a diffExprGroups data.frame into a list of two exprGroups data.frames,
#' which can then be passed to `simphony()`.
#'
#' @param diffExprGroups `data.frame` with optional columns `meanBase`,
#'   `dBase`, `meanSd`, `dSd`, `meanAmp`, `dAmp`, `meanPhase`, and `dPhase`
#'   describing the changes in expression between two conditions. Each row
#'   corresponds to a group of genes.
#' @param checkValid Logical for whether to only return rows for which both
#'   amplitudes are greater than or equal to zero and both standard deviations
#'   are greater than zero.
#'
#' @return List of two `data.table`s with possible columns `base`, `sd`, `amp`,
#'   and `phase`, depending on the columns in `diffExprGroups`.
#'
#' @examples
#' dGroups = data.frame(meanAmp = c(1, 1, 1, 1), dAmp = c(1, 1, 2, 2),
#'                      meanPhase = c(0, 0, 0, 0), dPhase = c(0, 3, 0, 3))
#' exprGroups = splitDiffExprGroups(dGroups)
#'
#' @seealso `\link{simphony}`
#'
#' @export
splitDiffExprGroups = function(diffExprGroups, checkValid = TRUE) {
  dGroups = data.table(diffExprGroups)

  capCols = c('Base', 'Amp', 'Phase', 'Sd')
  cols = tolower(capCols)
  meanCols = paste0('mean', capCols)
  dCols = paste0('d', capCols)

  d1 = data.table(.dummy = rep(1, nrow(dGroups)))
  d2 = data.table(.dummy = rep(1, nrow(dGroups)))
  for (ii in 1:length(cols)) {
    if (all(c(meanCols[ii], dCols[ii]) %in% colnames(dGroups))) {
      d1[[cols[ii]]] = dGroups[[meanCols[ii]]] - 0.5 * dGroups[[dCols[ii]]]
      d2[[cols[ii]]] = dGroups[[meanCols[ii]]] + 0.5 * dGroups[[dCols[ii]]]}}
  d1[, .dummy := NULL]
  d2[, .dummy := NULL]

  heldbackCols = setdiff(colnames(dGroups), c(meanCols, dCols))
  if (length(heldbackCols) > 0) {
    dHeldback = dGroups[, heldbackCols, with = FALSE]
    d1 = cbind(d1, dHeldback)
    d2 = cbind(d2, dHeldback)}

  if (checkValid) {
    idx = rep(TRUE, nrow(d1))
    if ('amp' %in% colnames(d1)) {
      idx = idx & (d1$amp >= 0) & (d2$amp >= 0)}
    if ('sd' %in% colnames(d1)) {
      idx = idx & (d1$sd > 0) & (d2$sd > 0)}
    d1 = d1[idx]
    d2 = d2[idx]}

  return(list(d1, d2))}
