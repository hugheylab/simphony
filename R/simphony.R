#' @importFrom data.table data.table ":="
#' @importFrom foreach foreach "%do%"
NULL


globalVariables(c('base', 'amp', 'phase', 'group', 'rhyFunc', 'sd', 'cond',
                  'dispFunc', 'featureGroups', 'time', 'abund', '..cond', '.N',
                  '.dummy', 'feature', 'mu', '..period'))


#' Simulate feature abundance data
#'
#' Simulate experiments in which rhythmic and non-rhythmic feature abundance is
#' measured at multiple timepoints in one or more conditions.
#'
#' @param featureGroupsList `data.frame` or `data.table` (for a single condition)
#'   or list of `data.frame`s or `data.table`s (for multiple conditions), where
#'   each row corresponds to a group of features to simulate. The following
#'   columns are all optional:
#'   \describe{
#'     \item{fracFeatures}{Fraction of simulated features to allocate to each group.
#'       Defaults to 1/(number of groups).}
#'     \item{base}{Average abundance. Defaults to 0 if `family` == 'gaussian'
#'       and to 8 (mean log2 counts) if `family` == 'negbinom'.}
#'     \item{sd}{Standard deviation of sampled abundance values. Defaults to 1.
#'       Only used if `family` == 'gaussian'.}
#'     \item{dispFunc}{Function to calculate dispersion of sampled abundance
#'       values, given expected abundance in counts. Only used if `family` ==
#'       'negbinom'.}
#'     \item{amp}{Amplitude of rhythmic abundance. Defaults to 0 (i.e.,
#'       non-rhythmic.}
#'     \item{phase}{Phase of rhythmic abundance, in the same units as `period`.
#'       Defaults to 0.}
#'     \item{rhyFunc}{Function to generate rhythmic abundance. Must have a
#'       period of 2*pi. Defaults to `sin`.}
#'   }
#' @param fracFeatures Fraction of simulated features to allocate to each group.
#'   Defaults to 1/(number of groups). Only used if the first `featureGroupsList`
#'   `data.frame` lacks a `fracFeatures` column.
#' @param nFeatures Integer for the total number of features to simulate.
#' @param period Integer for the period of simulated rhythms.
#' @param timepointsType Character string for how to set the timepoints
#'   for the simulation. Must be 'auto' (default), 'specified', or 'random'.
#' @param timeRange Optional 2-element vector controlling the range of times to
#'   sample in simulated data. Defaults to c(0, period).
#' @param interval Integer for the amount of time between consecutive
#'   timepoints, in the same units as `period`. The first timepoint is 0. Only
#'   used if `timepointsType` == 'auto'.
#' @param nReps Integer for the number of replicates per timepoint. Only used
#'   if `timepointsType` == 'auto'.
#' @param timepoints Numeric vector of exact timepoints to simulate, including
#'   any replicates. Only used if `timepointsType` == 'specified'.
#' @param nSamplesPerCond Integer for the number of samples per condition,
#'   which will be randomly uniformly spaced between 0 and `period` and
#'   different for each condition. Only used if timepointsType == 'random'.
#' @param dispFunc Function to calculate dispersion of sampled abundance
#'   values, given expected abundance in counts. Defaults to `defaultDispFunc`.
#'   Only used if `family` == 'negbinom' and a `data.frame` in `featureGroupsList`
#'   lacks a `dispFunc` column.
#' @param rhyFunc Function to generate rhythmic abundance. Must have a period
#'   of 2*pi. Defaults to `sin`. Only used if a `data.frame` in `featureGroupsList`
#'   lacks a `rhyFunc` column.
#' @param family Character string for the family of distributions from
#'   which to generate the abundance values. Must be 'gaussian' or 'negbinom'.
#'
#' @return List with the following elements:
#' \describe{
#'   \item{abundData}{Matrix of abundance values (counts, if `family` ==
#'   'negbinom'), with features as rownames and samples as colnames.}
#'   \item{sampleMetadata}{`data.table` with one row per sample.}
#'   \item{featureMetadata}{`data.table` with one row per feature per condition.}
#'   \item{experMetadata}{List of arguments that were passed to `simphony`.}
#' }
#'
#' @examples
#' library('data.table')
#'
#' # Simulate data for features having one of three sets of rhythmic parameters.
#' featureGroups = data.table(amp = c(0, 1, 1), phase = c(0, 0, 6),
#'                         rhyFunc = c(cos, cos, sin))
#' simData = simphony(featureGroups)
#'
#' # Simulate data for an experiment with specified timepoints and replicates.
#' featureGroups = data.table(amp = c(0, 1))
#' simData = simphony(featureGroups, timepointsType = 'specified',
#'                    timepoints = rep(seq(0, 6, 2), each = 2))
#'
#' # Simulate data for features whose rhythmicity varies between two conditions.
#' featureGroupsList = list(data.table(amp = c(1, 2), phase = c(0, -3)),
#'                       data.table(amp = c(3, 2), phase = c(0, 3)))
#' simData = simphony(featureGroupsList)
#'
#' # Simulate data for 100 features, half non-rhythmic and half rhythmic, with
#' # amplitudes for rhythmic features sampled from a log-normal distribution.
#' nFeatures = 100
#' rhyFrac = 0.5
#' nRhyFeatures = round(rhyFrac * nFeatures)
#' rhyAmps = exp(rnorm(nRhyFeatures, mean = 0, sd = 0.25))
#' fracFeatures = c(1 - rhyFrac, rep(rhyFrac / nRhyFeatures, nRhyFeatures))
#' featureGroups = data.table(amp = c(0, rhyAmps), fracFeatures = fracFeatures)
#' simData = simphony(featureGroups, nFeatures = nFeatures)
#'
#' # Simulate data for 100 rhythmic features, with baseline log2 expected counts
#' # and residual log dispersion sampled from distributions whose parameters
#' # were estimated, using DESeq2 and fitdistrplus, from circadian RNA-seq data
#' # from mouse liver (PRJNA297287).
#' nFeatures = 100
#' baseLog2Counts = rnorm(nFeatures, mean = 8.63, sd = 2.73)
#' dispFactors = exp(rnorm(nFeatures, sd = 0.819))
#' dispFuncs = sapply(dispFactors, function(z) {function(x) defaultDispFunc(x) * z})
#' featureGroups = data.table(base = baseLog2Counts, dispFunc = dispFuncs, amp = 1)
#' simData = simphony(featureGroups, nFeatures = nFeatures, family = 'negbinom')
#'
#' @seealso `\link{defaultDispFunc}`, `\link{getExpectedAbund}`,
#' `\link{getSampledAbund}`, `\link{mergeSimData}`
#'
#' @export
simphony = function(featureGroupsList, fracFeatures = NULL, nFeatures = 10, period = 24,
                    timepointsType = c('auto', 'specified', 'random'),
                    timeRange = c(0, period),
                    interval = 2, nReps = 1, timepoints = NULL,
                    nSamplesPerCond = NULL, rhyFunc = sin,
                    dispFunc = defaultDispFunc,
                    family = c('gaussian', 'negbinom')) {
  family = match.arg(family)
  timepointsType = match.arg(timepointsType)

  if (is.data.frame(featureGroupsList)) {
    featureGroupsList = list(featureGroupsList)
  } else if (!is.list(featureGroupsList) || !all(sapply(featureGroupsList, is.data.frame))) {
    stop('featureGroupsList must be a data.frame or a list of data.frames.')}

  if (length(unique(sapply(featureGroupsList, nrow))) != 1) {
    stop('Each featureGroups data.frame must have the same number of rows.')}

  featureGroupsList = foreach(featureGroups = featureGroupsList) %do% {
    setDefaultFeatureGroups(featureGroups, nFeatures, dispFunc, rhyFunc, family)}

  times = getTimes(timepointsType, interval, nReps, timepoints, timeRange,
                   nSamplesPerCond, length(featureGroupsList))
  sm = getSampleMetadata(times)
  fm = getFeatureMetadata(featureGroupsList, fracFeatures, nFeatures)

  abundDt = getExpectedAbund(fm, period, sampleMetadata = sm)
  abundMat = getSampledAbund(abundDt, family, inplace = TRUE)

  experMetadata = list(featureGroupsList = featureGroupsList, fracFeatures = fracFeatures,
                       nFeatures = nFeatures, period = period,
                       timepointsType = timepointsType, interval = interval,
                       nReps = nReps, timepoints = timepoints, rhyFunc = rhyFunc,
                       nSamplesPerCond = nSamplesPerCond,
                       dispFunc = dispFunc, family = family)

  return(list(abundData = abundMat, sampleMetadata = sm, featureMetadata = fm,
              experMetadata = experMetadata))}
