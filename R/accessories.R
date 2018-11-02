#' Calculate dispersion given expected counts
#'
#' Calculate dispersion for a negative binomial distribution given expected
#' counts. Trends were estimated by using DESeq2 and circadian RNA-seq data from
#' mice and humans.
#'
#' @param x Numeric vector of expected expression counts.
#'
#' @return Numeric vector of dispersions.
#'
#' @examples
#'
#' @seealso `\link{simphony}`
#'
#' @export
defaultDispFunc = function(x) {
  return(3/x)}


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
#' @param genes Character vector of genes for which to get expression data.
#'
#' @return `data.table`.
#'
#' @examples
#'
#' @seealso `\link{simphony}`
#'
#' @export
mergeSimData = function(simData, genes) {
  d = data.table(t(simData$exprData[genes, ]), keep.rownames = TRUE)
  d = merge(data.table(simData$sampleMetadata), d, by.x = 'sample', by.y = 'rn')
  d = data.table::melt(d, measure.vars = genes, variable.name = 'gene',
                       value.name = 'expr')
  d = merge(d, simData$geneMetadata, by = c('gene', 'cond'))
  return(d)}


#' Split differential exprGroups
#'
#' Split a diffExprGroups data.frame into a list of two exprGroups data.frames,
#' which can then be passed to `simphony()`.
#'
#' @param diffExprGroups `data.frame` with optional columns `meanBase`,
#' `dBase`, `meanSd`, `dSd`, `meanAmp`, `dAmp`, `meanPhase`, and `dPhase`
#' describing the changes in expression between two conditions. Each row
#' corresponds to a group of genes.
#' @param checkValid Logical for whether to only return rows for which both
#' amplitudes are greater than or equal to zero and both standard deviations are
#' greater than zero.
#'
#' @return List of two `data.table`s with possible columns `base`, `sd`,
#' `amp`, and `phase`, depending on the columns in `diffExprGroups`.
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
