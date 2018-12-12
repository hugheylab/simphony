#' @importFrom data.table data.table ":="
#' @importFrom foreach foreach "%do%"
NULL


globalVariables(c('base', 'amp', 'phase', 'group', 'rhyFunc', 'sd', 'cond',
                  'dispFunc', 'exprGroups', 'time', 'expr', '..cond', '.N',
                  '.dummy', 'gene', 'mu', '..period'))


#' Simulate gene expression data
#'
#' Simulate experiments in which rhythmic and non-rhythmic gene expression is
#' measured at multiple timepoints in one or more conditions.
#'
#' @param exprGroupsList `data.frame` or `data.table` (for a single condition)
#'   or list of `data.frame`s or `data.table`s (for multiple conditions), where
#'   each row corresponds to a group of genes to simulate. The following
#'   columns are all optional:
#'   \describe{
#'     \item{fracGenes}{Fraction of simulated genes to allocate to each group.
#'       Defaults to 1/(number of groups).}
#'     \item{base}{Average expression. Defaults to 0 if `family` == 'gaussian'
#'       and to 7 (mean log2 counts) if `family` == 'negbinom'.}
#'     \item{sd}{Standard deviation of sampled expression values. Defaults to 1.
#'       Only used if `family` == 'gaussian'.}
#'     \item{dispFunc}{Function to calculate dispersion of sampled expression
#'       values, given expected expression in counts. Only used if `family` ==
#'       'negbinom'.}
#'     \item{amp}{Amplitude of rhythmic expression. Defaults to 0 (i.e.,
#'       non-rhythmic.}
#'     \item{phase}{Phase of rhythmic expression, in the same units as `period`.
#'       Defaults to 0.}
#'     \item{rhyFunc}{Function to generate rhythmic expression. Must have a
#'       period of 2*pi. Defaults to `sin`.}
#'   }
#' @param fracGenes Fraction of simulated genes to allocate to each group.
#'   Defaults to 1/(number of groups). Only used if the first `exprGroupsList`
#'   `data.frame` lacks a `fracGenes` column.
#' @param nGenes Integer for the total number of genes to simulate.
#' @param period Integer for the period of simulated rhythms.
#' @param timepointsType Character string for how to set the timepoints
#'   for the simulation. Must be 'auto' (default), 'specified', or 'random'.
#' @param interval Integer for the amount of time between consecutive
#'   timepoints, in the same units as `period`. The first timepoint is 0. Only
#'   used if `timepointsType` == 'auto'.
#' @param nReps Integer for the number of replicates per timepoint. Only used
#'   if `timepointsType` == 'auto'.
#' @param timepoints Numeric vector of exact timepoints to simulate, including
#'   any replicates. Only used if `timepointsType` == 'specified'.
#' @param nSamplesPerCond Integer for the number of samples per condition,
#'   which will be randomly uniformly spaced between 0 and `period` and different
#'   for each condition. Only used if timepointsType == 'random'.
#' @param dispFunc Function to calculate dispersion of sampled expression
#'   values, given expected expression in counts. Defaults to `defaultDispFunc`.
#'   Only used if `family` == 'negbinom' and a `data.frame` in `exprGroupsList`
#'   lacks a `dispFunc` column.
#' @param rhyFunc Function to generate rhythmic expression. Must have a period
#'   of 2*pi. Defaults to `sin`. Only used if a `data.frame` in `exprGroupsList`
#'   lacks a `rhyFunc` column.
#' @param family Character string for the family of distributions from
#'   which to generate the expression values. Must be 'gaussian' or 'negbinom'.
#'
#' @return List with the following elements:
#' \describe{
#'   \item{exprData}{Matrix of expression values (counts, if `family` ==
#'   'negbinom'), with genes as rownames and samples as colnames.}
#'   \item{sampleMetadata}{`data.table` with one row per sample.}
#'   \item{geneMetadata}{`data.table` with one row per gene per condition.}
#'   \item{simMetadata}{List of arguments that were passed to `simphony`.}
#' }
#'
#' @examples
#' library('data.table')
#'
#' # Simulate data for genes having one of three sets of rhythmic parameters.
#' exprGroups = data.table(amp = c(0, 1, 1), phase = c(0, 0, 6),
#'                         rhyFunc = c(cos, cos, sin))
#' simData = simphony(exprGroups)
#'
#' # Simulate data for an experiment with specified timepoints and replicates.
#' exprGroups = data.table(amp = c(0, 1))
#' simData = simphony(exprGroups, timepointsType = 'specified',
#'                    timepoints = rep(seq(0, 6, 2), each = 2))
#'
#' # Simulate data for genes whose rhythmicity varies between two conditions.
#' exprGroupsList = list(data.table(amp = c(1, 2), phase = c(0, -3)),
#'                       data.table(amp = c(3, 2), phase = c(0, 3)))
#' simData = simphony(exprGroupsList)
#'
#' # Simulate data for 100 genes, half non-rhythmic and half rhythmic, with
#' # amplitudes for rhythmic genes sampled from a distribution whose parameters
#' # were estimated, using limma-voom (q <= 0.01) and fitdistrplus, from
#' # circadian RNA-seq data from mouse liver (PRJNA297287).
#' nGenes = 100
#' rhyFrac = 0.5
#' nRhyGenes = round(rhyFrac * nGenes)
#' rhyAmps = 2^rnorm(nRhyGenes, mean = -0.278, sd = 0.563)
#' fracGenes = c(1 - rhyFrac, rep(rhyFrac / nRhyGenes, nRhyGenes))
#' exprGroups = data.table(amp = c(0, rhyAmps), fracGenes = fracGenes)
#' simData = simphony(exprGroups, nGenes = nGenes, family = 'negbinom')
#'
#' # Simulate data for 100 rhythmic genes, with baseline log2 expected counts
#' # and residual log dispersion sampled from distributions whose parameters
#' # were estimated, using DESeq2 and fitdistrplus, from circadian RNA-seq data
#' # from mouse liver (PRJNA297287).
#' nGenes = 100
#' baseLog2Counts = rnorm(nGenes, mean = 8.63, sd = 2.73)
#' dispFactors = exp(rnorm(nGenes, sd = 0.819))
#' dispFuncs = sapply(dispFactors, function(z) {function(x) defaultDispFunc(x) * z})
#' exprGroups = data.table(base = baseLog2Counts, dispFunc = dispFuncs, amp = 1)
#' simData = simphony(exprGroups, nGenes = nGenes, family = 'negbinom')
#'
#' @seealso `\link{defaultDispFunc}`, `\link{getExpectedExpr}`,
#' `\link{getSampledExpr}`, `\link{mergeSimData}`
#'
#' @export
simphony = function(exprGroupsList, fracGenes = NULL, nGenes = 10, period = 24,
                    timepointsType = 'auto', interval = 2, nReps = 1,
                    timepoints = NULL, nSamplesPerCond = NULL, rhyFunc = sin,
                    dispFunc = defaultDispFunc, family = 'gaussian') {
  if (!family %in% c('gaussian', 'negbinom')) {
    stop("family must be 'gaussian' or 'negbinom'.")}

  if (is.data.frame(exprGroupsList)) {
    exprGroupsList = list(exprGroupsList)
  } else if (!is.list(exprGroupsList) || !all(sapply(exprGroupsList, is.data.frame))) {
    stop('exprGroupsList must be a data.frame or a list of data.frames.')}

  if (length(unique(sapply(exprGroupsList, nrow))) != 1) {
    stop('Each exprGroups data.frame must have the same number of rows.')}

  exprGroupsList = foreach(exprGroups = exprGroupsList) %do% {
    setDefaultExprGroups(exprGroups, nGenes, dispFunc, rhyFunc, family)}

  times = getTimes(timepointsType, interval, nReps, timepoints,
                   nSamplesPerCond, length(exprGroupsList), period)
  sm = getSampleMetadata(times)
  gm = getGeneMetadata(exprGroupsList, fracGenes, nGenes)

  exprDt = getExpectedExpr(gm, period, sampleMetadata = sm)
  exprMat = getSampledExpr(exprDt, family, inplace = TRUE)

  # call = sys.call()
  simMetadata = list(exprGroupsList = exprGroupsList, fracGenes = fracGenes,
                     nGenes = nGenes, period = period,
                     timepointsType = timepointsType, interval = interval,
                     nReps = nReps, timepoints = timepoints,
                     nSamplesPerCond = nSamplesPerCond, rhyFunc = rhyFunc,
                     dispFunc = dispFunc, family = family)

  return(list(exprData = exprMat, sampleMetadata = sm, geneMetadata = gm,
              simMetadata = simMetadata))}
