#' @importFrom data.table data.table ":="
#' @importFrom foreach foreach "%do%"
globalVariables(c('base', 'amp', 'phase', 'group', 'rhyFunc', 'sd', 'cond',
                  'dispFunc', 'exprGroups', 'time', 'expr', '..cond', '.N',
                  '.dummy', 'gene', 'mu', '..period'))


#' Simulate gene expression data
#'
#' Simulate experiments in which rhythmic and non-rhythmic gene expression is
#' measured at multiple timepoints in one or more conditions.
#'
#' @param exprGroupsList `data.frame` or `data.table` (for a single condition)
#' or list of `data.frame`s or `data.table`s (for multiple conditions), where
#' each row corresponds to a group of genes to simulate. The following
#' columns are all optional:
#' \describe{
#'   \item{fracGenes}{Fraction of simulated genes to allocate to each group.
#'   Defaults to 1/(number of groups).}
#'   \item{base}{Average expression. Defaults to 0 if `family` == 'gaussian'
#'   and to 7 (mean log2 counts) if `family` == 'negbinom'.}
#'   \item{sd}{Standard deviation of sampled expression values. Defaults to 1.
#'   Only used if `family` == 'gaussian'.}
#'   \item{dispFunc}{Function to calculate dispersion of sampled expression
#'   values, given expected expression in counts. Only used if `family` ==
#'   'negbinom'.}
#'   \item{amp}{Amplitude of rhythmic expression. Defaults to 0 (i.e.,
#'   non-rhythmic.}
#'   \item{phase}{Phase of rhythmic expression, in the same units as `period`.
#'   Defaults to 0.}
#'   \item{rhyFunc}{Function to generate rhythmic expression. Must have a period
#'   of 2*pi. Defaults to `sin`.}
#' }
#' @param fracGenes Fraction of simulated genes to allocate to each group.
#' Defaults to 1/(number of groups). Only used if the first `exprGroupsList`
#' `data.frame` lacks a `fracGenes` column.
#' @param nGenes Integer for the total number of genes to simulate.
#' @param period Integer for the period of simulated rhythms.
#' @param timepointsType Character string for how to set the timepoints
#' for the simulation. Must be 'auto' (default), 'specified', or 'random'.
#' @param interval Integer for the amount of time between consecutive
#' timepoints, in the same units as `period`. The first timepoint is 0. Only
#' used if `timepointsType` == 'auto'.
#' @param nReps Integer for the number of replicates per timepoint. Only used
#' if `timepointsType` == 'auto'.
#' @param timepoints Numeric vector of exact timepoints to simulate, including
#' any replicates. Only used if `timepointsType` == 'specified'.
#' @param nSamplesPerCond Integer for the number of samples per condition,
#' which will be randomly uniformly spaced between 0 and `period` and different
#' for each condition. Only used if timepointsType == 'random'.
#' @param dispFunc Function to calculate dispersion of sampled expression
#' values, given expected expression in counts. Defaults to `defaultDispFunc`.
#' Only used if `family` == 'negbinom' and a `data.frame` in `exprGroupsList`
#' lacks a `dispFunc` column.
#' @param rhyFunc Function to generate rhythmic expression. Must have a period
#' of 2*pi. Defaults to `sin`. Only used if a `data.frame` in `exprGroupsList`
#' lacks a `rhyFunc` column.
#' @param family Character string for the family of distributions from
#' which to generate the expression values. Must be 'gaussian' or 'negbinom'.
#'
#' @return List with the following elements:
#' \describe{
#'   \item{exprData}{Matrix of expression values (counts, if `family` ==
#'   'negbinom'), with rownames for genes and colnames for samples.}
#'   \item{sampleMetadata}{`data.table` with one row per sample.}
#'   \item{geneMetadata}{`data.table` with one row per gene per condition.}
#' }
#'
#' @examples
#'
#' @seealso `\link{defaultDispFunc}`
#'
#' @export
simphony = function(exprGroupsList, fracGenes = NULL, nGenes = 10, period = 24,
                    timepointsType = 'auto', interval = 2, nReps = 1,
                    timepoints = NULL, nSamplesPerCond = NULL,
                    dispFunc = defaultDispFunc, rhyFunc = sin,
                    family = 'gaussian') {
  if (!family %in% c('gaussian', 'negbinom')) {
    stop("family must be 'gaussian' or 'negbinom'.")}

  if (is.data.frame(exprGroupsList)) {
    exprGroupsList = list(exprGroupsList)}
  if (length(unique(sapply(exprGroupsList, nrow))) != 1) {
    stop('Each exprGroups data.frame must have the same number of rows.')}

  exprGroupsList = foreach(exprGroups = exprGroupsList) %do% {
    setDefaultExprGroups(exprGroups, nGenes, dispFunc, rhyFunc, family)}
  nCond = length(exprGroupsList)

  nGenesPerGroup = getNGenesPerGroup(exprGroupsList[[1]], fracGenes, nGenes)

  times = getTimes(timepointsType, interval, nReps, timepoints,
                   nSamplesPerCond, nCond, period)

  genes = sprintf(sprintf('gene_%%0%dd', floor(log10(nGenes)) + 1), 1:nGenes)

  gm = foreach(exprGroups = exprGroupsList, cond = 1:nCond, .combine = rbind) %do% {
    gmNow = exprGroups[rep(1:.N, times = nGenesPerGroup)]
    gmNow[, cond := ..cond]
    gmNow[, gene := genes]
    data.table::setcolorder(gmNow, c('cond', 'group', 'gene'))
    gmNow}

  nSamples = prod(dim(times))
  nSamplesPerCond = ncol(times)
  sm = foreach(cond = 1:nCond, .combine = rbind) %do% {
    sampleIds = ((cond - 1) * nSamplesPerCond + 1):(nSamplesPerCond * cond)
    sampleNames = sprintf(sprintf('sample_%%0%dd', floor(log10(nSamples)) + 1),
                          sampleIds)
    data.table(sample = sampleNames, cond = cond, time = times[cond, ])}

  exprDt = getExpectedExpr(gm, period, sampleMetadata = sm)
  exprMat = getObservedExpr(exprDt, family, inplace = TRUE)

  return(list(exprData = exprMat, sampleMetadata = sm, geneMetadata = gm))}


#' Calculate expected expression
#'
#' Calculate expected expression for multiple genes at multiple timepoints in
#' multiple conditions.
#'
#' @param geneMetadata `data.table` with columns `base`, `rhyFunc`, `amp`, and
#' `phase`.
#' @param period Integer for the period of simulated rhythms.
#' @param times Numeric vector of the times (in the same units as `period`) at
#' which to calculate expected expression for each row in `geneMetadata`.
#' @param sampleMetadata `data.table` with columns `sample`, `cond`, and
#' `time`. Either `times` or `sampleMetadata` must be provided, and the former
#' takes precedence.
#' @param byCondGroup Logical for whether to speed up the calculation by
#' grouping by the columns `cond` and `group`. Primarily for internal use.
#'
#' @return `data.table` derived from `geneMetadata` (but with more rows),
#' with additional columns `time` and `mu` and possibly others. If sampling
#' will use the negative binomial family, `mu` corresponds to log2 counts.
#'
#' @examples
#'
#' @seealso `\link{simphony}`, `\link{getObservedExpr}`
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
#' @param exprDt `data.table`. If `family` == 'gaussian', required columns are
#' `gene`, `sample`, `mu`, and `sd`. If `family` == 'negbinom', required columns
#' are `gene`, `sample`, `mu`, `dispFunc`, `cond`, and `group`.
#' @param family Character string for the family of distributions from which
#' to generate the expression values. Must be 'gaussian' or 'negbinom'.
#' @param inplace Logical for whether to modify in-place `exprDt`, adding a
#' column `expr` containing the expression values.
#'
#' @return Matrix of expression values, where rows correspond to genes and
#' columns correspond to samples.
#'
#' @examples
#'
#' @seealso `\link{simphony}`, `\link{getExpectedExpr}`
#'
#' @export
getObservedExpr = function(exprDt, family, inplace = FALSE) {
  if (!inplace) {
    exprDt = data.table(exprDt)}

  if (family == 'gaussian') {
    exprDt[, expr := stats::rnorm(.N, mu, sd)]
  } else {
    # dispFunc is identical for genes of the same group in the same condition
    # this is the way I've figured out how to call functions that are columns
    exprDt[, expr := stats::rnbinom(.N, mu = 2^mu, size = 1/dispFunc[[1]](2^mu)),
           by = c('cond', 'group')]}

  exprDtCast = data.table::dcast(exprDt, gene ~ sample, value.var = 'expr')
  exprMat = as.matrix(exprDtCast, rownames = 1)
  return(exprMat)}
