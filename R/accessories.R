#' Default function for mapping expected counts to dispersion.
#'
#' The function was estimated from circadian RNA-seq data from mouse liver
#' (PRJNA297287), using local regression in DESeq2. In a negative binomial
#' distribution, variance = mean + mean^2 * dispersion.
#'
#' @param x Numeric vector of mean counts.
#'
#' @return Numeric vector of dispersions.
#'
#' @examples
#' means = 2^(6:10)
#' dispersions = defaultDispFunc(means)
#'
#' @seealso `\link{simphony}`
'defaultDispFunc'


#' Calculate expected abundance
#'
#' Calculate expected abundance for multiple features at multiple timepoints in
#' multiple conditions.
#'
#' @param featureMetadata `data.table` with columns `feature`, `base`,
#'   `rhyFunc`, `amp`, `period`, and `phase`, where every row corresponds to a
#'   gene. If `byCondGroup` is `TRUE`, then must also have columns `cond` and
#'   `group`.
#' @param times Numeric vector of the times at which to calculate expected
#'   abundance for each row in `featureMetadata`.
#' @param sampleMetadata `data.table` with columns `sample`, `cond`, and
#'   `time`. Either `times` or `sampleMetadata` must be provided, and the former
#'   takes precedence.
#' @param byCondGroup Logical for whether to speed up the calculation by
#'   grouping by the columns `cond` and `group`. Primarily for internal use.
#'
#' @return `data.table` derived from `featureMetadata` (but with more rows),
#'   with additional columns `time` and `mu` and possibly others. If sampling
#'   will use the negative binomial family, `mu` corresponds to log2 counts.
#'
#' @examples
#' library('data.table')
#' featureMetadata = data.table(feature = c('feature_1', 'feature_2'),
#'                              base = function(x) 0,
#'                              amp = c(function(x) 0, function(x) 1),
#'                              period = 24,
#'                              phase = 0, rhyFunc = sin)
#' abundDt = getExpectedAbund(featureMetadata, times = 6:17)
#'
#' @seealso `\link{simphony}`, `\link{getSampledAbund}`
#'
#' @export
getExpectedAbund = function(featureMetadata, times = NULL,
                            sampleMetadata = NULL,
                            byCondGroup = is.null(times)) {
  if (!is.null(times)) {
    d = data.table(featureMetadata)[rep(1:.N, each = length(times))]
    d[, time := rep(times, times = nrow(featureMetadata))]
  } else if (!is.null(sampleMetadata)) {
    d = merge(data.table(featureMetadata), sampleMetadata, by = 'cond',
              allow.cartesian = TRUE)
  } else {
    stop('Either times or sampleMetadata must not be NULL.')}

  if (isTRUE(byCondGroup)) {
    d[, mu := base[[1]](time) +
              amp[[1]](time) * rhyFunc[[1]]((time + phase) * 2 * pi / period),
      by = c('cond', 'group')]
  } else {
    d[, mu := base[[1]](time) +
              amp[[1]](time) * rhyFunc[[1]]((time + phase) * 2 * pi / period),
      by = 1:nrow(d)]}
  return(data.table::copy(d))}


#' Sample abundance values
#'
#' Sample feature abundance values from the given distributions. This function
#' is used internally by `simphony()`, and should not usually need to be
#' called directly.
#'
#' @param abundDt `data.table` of expected abundance. If `family` is 'gaussian',
#'   required columns are `feature`, `sample`, `mu`, and `sd`. If `family` is
#'   'negbinom', required columns are `feature`, `sample`, `mu`, `dispFunc`,
#'   `cond`, and `group`. If `family` is 'bernoulli' or 'poisson', required
#'   columns are `feature`, `sample`, and `mu`.
#' @param family Character string for the family of distributions from which
#'   to sample the abundance values.
#' @param inplace Logical for whether to modify `abundDt` in-place, adding a
#'   column `abund` containing the abundance values.
#'
#' @return Matrix of abundance values, where rows correspond to features and
#'   columns correspond to samples.
#'
#' @examples
#' library('data.table')
#' set.seed(6022)
#' abundDt = data.table(feature = 'feature_1', sample = c('sample_1', 'sample_2'),
#'                     mu = c(0, 5), sd = 1)
#' abundMat = getSampledAbund(abundDt)
#'
#' @seealso `\link{simphony}`, `\link{getExpectedAbund}`
#'
#' @export
getSampledAbund = function(abundDt,
                           family = c('gaussian', 'negbinom', 'bernoulli', 'poisson'),
                           inplace = FALSE) {
  family = match.arg(family)
  if (isFALSE(inplace)) {
    abundDt = data.table(abundDt)}

  if (family == 'gaussian') {
    abundDt[, abund := stats::rnorm(.N, mu, sd)]
  } else if (family == 'negbinom') {
    # dispFunc is identical for features of the same group in the same condition
    # this is the way I've figured out how to call functions that are columns
    abundDt[, abund := stats::rnbinom(.N, mu = 2^mu, size = 1 / dispFunc[[1]](2^mu)),
           by = c('cond', 'group')]
  } else if (family == 'bernoulli') {
    # will output NA and a warning for mu < 0 or mu > 1
    abundDt[, abund := stats::rbinom(.N, 1, mu)]
  } else if (family == 'poisson') {
    # will output NA and a warning for mu < 0
    abundDt[, abund := stats::rpois(.N, mu)]}

  data.table::setorderv(abundDt, c('sample', 'feature'))
  features = unique(abundDt$feature)
  samples = unique(abundDt$sample)
  abundMat = matrix(abundDt$abund, nrow = length(features),
                    dimnames = list(features, samples))
  return(abundMat)}


#' Merge abundance data, feature metadata, and sample metadata
#'
#' Merge a simulation's abundance data, feature metadata, and sample metadata
#' into one `data.table`. This function is useful for making plots using
#' ggplot2.
#'
#' @param simData List with the following elements, such as returned by
#' `simphony()`:
#' \describe{
#'   \item{abundData}{Matrix of abundance values, with rownames for features and
#'   colnames for samples.}
#'   \item{sampleMetadata}{`data.table` with columns `sample` and `cond`.}
#'   \item{featureMetadata}{`data.table` with columns `feature` and `cond`.}
#' }
#' @param features Character vector of features for which to get abundance data.
#'   If NULL, then all features.
#'
#' @return `data.table`.
#'
#' @examples
#' library('data.table')
#' featureGroups = data.table(amp = c(0, 1))
#' simData = simphony(featureGroups)
#' mergedSimData = mergeSimData(simData, simData$featureMetadata$feature[1:2])
#'
#' @seealso `\link{simphony}`
#'
#' @export
mergeSimData = function(simData, features = NULL) {
  if (is.null(features)) {
    features = rownames(simData$abundData)}

  d = data.table(simData$abundData[features, , drop = FALSE],
                 keep.rownames = TRUE)
  data.table::setnames(d, 'rn', 'feature')
  d = data.table::melt(d, id.vars = 'feature', variable.name = 'sample',
                       value.name = 'abund')

  d = merge(d, simData$sampleMetadata, by = 'sample')
  d = merge(d, simData$featureMetadata, by = c('feature', 'cond'))
  return(d)}


#' Split differential featureGroups
#'
#' Split a diffFeatureGroups data.frame into a list of two featureGroups
#' data.frames, which can then be passed to `simphony()`.
#'
#' @param diffFeatureGroups `data.frame` with optional columns `meanBase`,
#'   `dBase`, `meanSd`, `dSd`, `meanAmp`, `dAmp`, `meanPhase`, and `dPhase`
#'   describing the changes in abundance between two conditions. Each row
#'   corresponds to a group of features.
#' @param checkValid Logical for whether to only return rows for which both
#'   amplitudes are greater than or equal to zero and both standard deviations
#'   are greater than zero.
#'
#' @return List of two `data.table`s with possible columns `base`, `sd`, `amp`,
#'   and `phase`, depending on the columns in `diffFeatureGroups`.
#'
#' @examples
#' dGroups = data.frame(meanAmp = c(1, 1, 1, 1), dAmp = c(1, 1, 2, 2),
#'                      meanPhase = c(0, 0, 0, 0), dPhase = c(0, 3, 0, 3))
#' featureGroups = splitDiffFeatureGroups(dGroups)
#'
#' @seealso `\link{simphony}`
#'
#' @export
splitDiffFeatureGroups = function(diffFeatureGroups, checkValid = TRUE) {
  dGroups = data.table(diffFeatureGroups)

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

  if (isTRUE(checkValid)) {
    idx = rep(TRUE, nrow(d1))
    if ('amp' %in% colnames(d1)) {
      idx = idx & (d1$amp >= 0) & (d2$amp >= 0)}
    if ('sd' %in% colnames(d1)) {
      idx = idx & (d1$sd > 0) & (d2$sd > 0)}
    d1 = d1[idx]
    d2 = d2[idx]}

  return(list(d1, d2))}
