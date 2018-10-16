#' @importFrom data.table ":="
#' @importFrom foreach "%do%"
globalVariables(c('base', 'amp', 'phase', 'group', 'rhyFunc', 'sd'))


sampleDispersion = function(x) {
  return(3/x)
}

getSingleCondSim = function(exprGroups, timePoints, method) {
  foreach::foreach(group = 1:nrow(exprGroups), .combine = rbind) %do% {
    amp = exprGroups[group, amp]
    phase = exprGroups[group, phase]
    base = exprGroups[group, base]
    rhyFunc = exprGroups[group, rhyFunc][[1]]
    mu = amp * rhyFunc(timePoints + 2 * pi * phase) + base

    if(method == 'gaussian') {
      samples = stats::rnorm(length(mu) * exprGroups[group, geneCount],
                             rep(mu, exprGroups[group, geneCount]),
                             sd = exprGroups[group, sd])
    } else {
      dispersionFunc = exprGroups[group, dispersionFunc][[1]]
      samples = stats::rnbinom(length(mu) * exprGroups[group, geneCount],
                               mu = 2^rep(mu, exprGroups[group, geneCount]),
                               size = 1/sampleDispersion(2^mu)) }
    matrix(samples, nrow = exprGroups[group, geneCount], byrow = TRUE)
  }
}

setOneCondDefault = function(exprGroups, nGenes, randomTimepoints, nSamples,
                             rhyFunc, method) {

  exprGroups = data.table::data.table(exprGroups)

  if(nrow(exprGroups) == 0) {
    stop('No rows in exprGroups. Cannot simulate genes.') }

  if(randomTimepoints && is.null(nSamples)) {
    stop('Number of random timepoint samples not specified.') }

  if(!'geneFrac' %in% colnames(exprGroups)) {
    exprGroups[, geneFrac := 1 / nrow(exprGroups)] }

  if(!'amp' %in% colnames(exprGroups)) {
    exprGroups[, amp := 0] }

  if(!'phase' %in% colnames(exprGroups)) {
    exprGroups[, phase := 0] }

  if(!'rhyFunc' %in% colnames(exprGroups)) {
    exprGroups[, rhyFunc := data.table::data.table(rhyFunc)] }

  if(method == 'negbinom') {
    if(!'dispersionFunc' %in% colnames(exprGroups)) {
      exprGroups[, dispersionFunc := data.table::data.table(sampleDispersion)] }
    if(!'base' %in% colnames(exprGroups)) {
    exprGroups[, base := 7] } }
  else {
    if(!'sd' %in% colnames(exprGroups)) {
      exprGroups[, sd := 1] } 
    if(!'base' %in% colnames(exprGroups)) {
      exprGroups[, base := 0] }}

  if(any(exprGroups[, geneFrac] <= 0)) {
    stop('All groups in exprGroups must have geneFrac > 0.') }

  # Compute a number of genes per group that sum to nGenes.
  exprGroups[, geneFrac := geneFrac / sum(geneFrac)]
  exprGroups[, geneCount := as.integer(geneFrac * nGenes)]
  if(sum(exprGroups[, geneCount]) != nGenes) {
    exprGroups[1L:(nGenes - sum(exprGroups[, geneCount])), geneCount := geneCount + 1]
  }

  return(exprGroups)
}

#' Generate list of two expression groups from a combined differential exprGroup
#' 
#'
#' @param twoCondGroups is the differential exprGroup to convert into two
#'   separate exprGroup data.table objects.
#' @examples
#'   dGroups = data.table::data.table(meanBase = c(0, 0, 1, 1), dBase = c(0, 0, 0.5, 0.5),
#'                                    meanAmp = c(1,2,1,2), dAmp = c(1,1,2,2),
#'                                    meanPhase = c(0, 0, 3, 3), dPhase = c(0, 0, 3, 3),
#'                                    meanSd = c(1, 1, 1, 1), dSd = c(0, 0, 0.5, 0.5))
#'   exprGroups = generateExprGroups(dGroups)
#' @export
generateExprGroups = function(twoCondGroups) {
  exprGroups = list(data.table::data.table(
                      base = twoCondGroups[, meanBase] + twoCondGroups[, dBase],
                      amp = twoCondGroups[, meanAmp] + twoCondGroups[, dAmp],
                      phase = twoCondGroups[, meanPhase] + twoCondGroups[, dPhase],
                      sd = twoCondGroups[, meanSd] + twoCondGroups[, dSd]),
                    data.table::data.table(
                      base = twoCondGroups[, meanBase] - twoCondGroups[, dBase],
                      amp = twoCondGroups[, meanAmp] - twoCondGroups[, dAmp],
                      phase = twoCondGroups[, meanPhase] - twoCondGroups[, dPhase],
                      sd = twoCondGroups[, meanSd] - twoCondGroups[, dSd]))
}

#' Generate simulated gene expression time courses
#'
#' @param exprGroupsList is a list of data.frame or data.table objects with the
#"   following optional columns:
#'   \itemize{
#'     \item{geneFrac}: {Fraction of all simulated genes which fall into this
#'                       group. Defaults to 1/nrow(exprGroups) if not supplied.}
#'     \item{meanBase}: {The mean baseline expression for this group. Defaults 
#'                       to 0 if not supplied.}
#'     \item{dBase}: {The difference in baseline expression across conditions
#'                    for this group. Defaults to 0 if not supplied.}
#'     \item{meanAmp}: {The mean amplitude of the rhythmic component of 
#'                      expression for this group. Defaults to 0 if not
#'                      supplied.}
#'     \item{dAmp}: {The difference in amplitude of the rhythmic component of 
#'                   expression across conditions for this group. Defaults to 0
#'                   if not supplied.}
#'     \item{meanPhase}: {The mean phase of the rhythmic component of expression
#'                        for this group. Defaults to 0 if not supplied.}
#'     \item{dPhase}: {The difference in phase of the rhythmic component of 
#'                     expression across conditions for this group. Defaults to
#'                     0 if not supplied.}
#'     \item{meanSd}: {The mean standard deviation of the sample error for this
#'                     group. Defaults to 1 if not supplied.}
#'     \item{dSd}: {The difference in standard deviation of the sample error for
#'                  this group. Defaults to 0 if not supplied.}
#'     \item{rhyFunc}: {The function used to generate the rhythmic component of
#'                      this group's gene expression. rhyFunc must have a period
#'                      of 2*pi. Defaults to sin if not supplied.}
#'   }
#' @param nGenes is the integer number of total genes to simulate.
#' @param period is the integer number of hours in one rhythmic cycle.
#' @param interval is the integer number of hours between simulated time points.
#' @param nReps is the integer number of replicates per time point.
#' @param nSamples is the integer number of time points to sample, if
#'   randomTimepoints is enabled. This must be supplied if randomTimepoints is
#'   TRUE.
#' @param randomTimepoints is a boolean determining whether to simulate an
#'   experiment with random sample times. Defaults to FALSE.
#' @param rhyFunc is the rhythmic function to set for exprGroups missing a
#'   rhythmic function. Defaults to sin if not supplied.
#' @param method is the data generation method to use. Must be either 'gaussian'
#'   or 'negbinom'.
#' @export
simulateGeneData = function(exprGroupsList, nGenes = 10, period = 24,
                            interval = 4, nReps = 2, nSamples = NULL,
                            randomTimepoints = FALSE, rhyFunc = sin,
                            method = 'gaussian') {
  if (!method %in% c('gaussian', 'negbinom')) {
    stop('Sample method must be either Gaussian or Negative Binomial') }

  exprGroupsList = foreach::foreach(exprGroups = exprGroupsList, .combine = list) %do% {
    setOneCondDefault(exprGroups, nGenes, randomTimepoints, nSamples, rhyFunc, method) }

  if(!is(exprGroupsList, 'list')) {exprGroupsList = list(exprGroupsList)}

  geneData = foreach::foreach(exprGroups = exprGroupsList, cond = 1:length(exprGroupsList), .combine = rbind) %do% {
    data.table::data.table(base = rep(exprGroups[, base], times = exprGroups[, geneCount]),
                           amp = rep(exprGroups[, amp], times = exprGroups[, geneCount]),
                           phase = rep(exprGroups[, phase], times = exprGroups[, geneCount]),
                           sd = ifelse('sd' %in% colnames(exprGroups),
                                       rep(exprGroups[, sd], times = exprGroups[, geneCount]),
                                       NA),
                           rhyFunc = rep(exprGroups[, rhyFunc], times = exprGroups[, geneCount]),
                           group = rep(1:nrow(exprGroups), times = exprGroups[, geneCount]),
                           cond = cond, gene = paste('gene', 1:nGenes, sep = '_'))
  }

  if(!randomTimepoints) {
    tt = (2 * pi / period) * interval * 0:(period %/% interval - (period %% interval == 0))
    tt = rep(tt, each = nReps)
    nSamples = length(tt)
    timePoints = matrix(rep(tt, each = length(exprGroupsList)), ncol = nSamples)
  } else {
    tt = stats::runif(nSamples * length(exprGroupsList), min = 0, max = 2 * pi)
    tt = matrix(tt, nrow = length(exprGroupsList), byrow = TRUE)
    timePoints = t(apply(tt, 1, sort))
    #timePoints = foreach::foreach(cond = 1:length(exprGroupsList), .combine = rbind) %do% {
    #  sort(stats::runif(nSamples, min = 0, max = 2 * pi)) }
  }

  sm = foreach::foreach(cond = 1:length(exprGroupsList), .combine = rbind) %do% {
    data.table::data.table(cond = cond, time = timePoints[cond, ] * period / (2*pi),
                           sample = paste('sample', ((cond-1)*nSamples+1):(cond*nSamples), sep = '_'))
  }

  emat = foreach::foreach(exprGroups = exprGroupsList, cond = 1:length(exprGroupsList), .combine = cbind) %do% {
    getSingleCondSim(exprGroups, timePoints[cond, ], method)
  }

  colnames(emat) = sm$sample
  rownames(emat) = paste('gene', 1:nGenes, sep = '_')

  return(list(emat = emat, sm = sm, gm = geneData))
}