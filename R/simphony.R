#' @importFrom data.table data.table ":="
#' @importFrom foreach foreach "%do%"
globalVariables(c('base', 'amp', 'phase', 'group', 'rhyFunc', 'sd', 'cond',
                  'dAmp', 'dBase', 'dSd', 'dispersionFunc', 'exprGroups',
                  'numGenes', 'fracGenes', 'meanAmp', 'meanBase', 'meanSd',
                  'dPhase', 'meanPhase'))


defaultDispersionFunc = function(x) {
  return(3/x)
}

simulateExprDataOneCond = function(exprGroups, times, method) {
  foreach(group = 1:nrow(exprGroups), .combine = rbind) %do% {
    amp = exprGroups[group, amp]
    phase = exprGroups[group, phase]
    base = exprGroups[group, base]
    rhyFunc = exprGroups[group, rhyFunc][[1]]
    mu = amp * rhyFunc(times + 2 * pi * phase) + base

    if(method == 'gaussian') {
      groupEmat = stats::rnorm(length(mu) * exprGroups[group, numGenes],
                               rep(mu, exprGroups[group, numGenes]),
                               sd = exprGroups[group, sd])
    } else {
      dispersionFunc = exprGroups[group, dispersionFunc][[1]]
      groupEmat = stats::rnbinom(length(mu) * exprGroups[group, numGenes],
                                 mu = 2^rep(mu, exprGroups[group, numGenes]),
                                 size = 1/dispersionFunc(2^mu)) }
    matrix(groupEmat, nrow = exprGroups[group, numGenes], byrow = TRUE)
  }
}

setDefaultExprGroups = function(exprGroups, nGenes, randomTimepoints, nSamples,
                                rhyFunc, method) {

  exprGroups = data.table(exprGroups)
  exprGroups[, group := 1:.N]

  if(nrow(exprGroups) == 0) {
    stop('No rows in exprGroups. Cannot simulate genes.') }

  if(randomTimepoints && is.null(nSamples)) {
    stop('Number of random timepoint samples not specified.') }

  if(!'fracGenes' %in% colnames(exprGroups)) {
    exprGroups[, fracGenes := 1 / nrow(exprGroups)] }

  if(!'amp' %in% colnames(exprGroups)) {
    exprGroups[, amp := 0] }

  if(!'phase' %in% colnames(exprGroups)) {
    exprGroups[, phase := 0] }

  if(!'rhyFunc' %in% colnames(exprGroups)) {
    exprGroups[, rhyFunc := data.table(rhyFunc)] }

  if(method == 'negbinom') {
    if(!'dispersionFunc' %in% colnames(exprGroups)) {
      exprGroups[, dispersionFunc := data.table(defaultDispersionFunc)] }
    if(!'base' %in% colnames(exprGroups)) {
      exprGroups[, base := 7] } }
  else {
    if(!'sd' %in% colnames(exprGroups)) {
      exprGroups[, sd := 1] }
    if(!'base' %in% colnames(exprGroups)) {
      exprGroups[, base := 0] }}

  if(any(exprGroups[, fracGenes] <= 0)) {
    stop('All groups in exprGroups must have fracGenes > 0.') }

  # Compute a number of genes per group that sum to nGenes.
  exprGroups[, fracGenes := fracGenes / sum(fracGenes)]
  exprGroups[, numGenes := as.integer(fracGenes * nGenes)]
  if(sum(exprGroups[, numGenes]) != nGenes) {
    exprGroups[1L:(nGenes - sum(exprGroups[, numGenes])), numGenes := numGenes + 1]
  }

  return(exprGroups)
}

#' Generate list of two expression groups from a combined differential exprGroup
#'
#'
#' @param diffExprGroups is the differential exprGroup to convert into two
#'   separate exprGroup data.table objects.
#' @examples
#'   dGroups = data.table::data.table(meanBase = c(0, 0, 1, 1), dBase = c(0, 0, 0.5, 0.5),
#'                                    meanAmp = c(1,2,1,2), dAmp = c(1,1,2,2),
#'                                    meanPhase = c(0, 0, 3, 3), dPhase = c(0, 0, 3, 3),
#'                                    meanSd = c(1, 1, 1, 1), dSd = c(0, 0, 0.5, 0.5))
#'   exprGroups = splitDiffExprGroups(dGroups)
#' @export
splitDiffExprGroups = function(diffExprGroups) {
  exprGroups = list(
    data.table(
      base = diffExprGroups[, meanBase] + diffExprGroups[, dBase],
      amp = diffExprGroups[, meanAmp] + diffExprGroups[, dAmp],
      phase = diffExprGroups[, meanPhase] + diffExprGroups[, dPhase],
      sd = diffExprGroups[, meanSd] + diffExprGroups[, dSd]),
    data.table(
      base = diffExprGroups[, meanBase] - diffExprGroups[, dBase],
      amp = diffExprGroups[, meanAmp] - diffExprGroups[, dAmp],
      phase = diffExprGroups[, meanPhase] - diffExprGroups[, dPhase],
      sd = diffExprGroups[, meanSd] - diffExprGroups[, dSd]))
}

#' Generate simulated gene expression time courses
#'
#' @param exprGroupsList is a list of data.frame or data.table objects with the
#"   following optional columns:
#'   \itemize{
#'     \item{fracGenes}: {Fraction of all simulated genes which fall into this
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
simulateExprData = function(exprGroupsList, nGenes = 10, period = 24,
                            interval = 4, nReps = 2, nSamples = NULL,
                            randomTimepoints = FALSE, rhyFunc = sin,
                            method = 'gaussian') {
  if (!method %in% c('gaussian', 'negbinom')) {
    stop('Sample method must be either gaussian or negbinom') }

  if(is.data.frame(exprGroupsList)) {
    exprGroupsList = list(setDefaultExprGroups(exprGroupsList, nGenes,
                                            randomTimepoints, nSamples, rhyFunc,
                                            method))
  } else {
    nGroups = sapply(exprGroupsList, nrow)
    if(length(unique(nGroups)) != 1) {
      stop('Number of rows in each exprGroups must be the same for all conditions') }
    exprGroupsList = foreach(exprGroups = exprGroupsList) %do% {
      setDefaultExprGroups(exprGroups, nGenes, randomTimepoints, nSamples, rhyFunc, method) }
  }

  geneNames = sprintf(sprintf('gene_%%0%dd', floor(log10(nGenes)) + 1), 1:nGenes)

  gm = foreach(exprGroups = exprGroupsList, cond = 1:length(exprGroupsList),
               .combine = rbind) %do% {
    gmNow = exprGroups[rep(1:.N, times = numGenes)]
    gmNow[, c('fracGenes', 'numGenes') := NULL]
    gmNow[, cond := ..cond]
    gmNow[, gene := geneNames]
    setcolorder(gmNow, c('cond', 'group', 'gene'))
    gmNow
  }

  if(!randomTimepoints) {
    tt = (2 * pi / period) * interval * 0:(period %/% interval - (period %% interval == 0))
    tt = rep(tt, each = nReps)
    nSamples = length(tt)
    times = matrix(rep(tt, each = length(exprGroupsList)), ncol = nSamples)
  } else {
    tt = stats::runif(nSamples * length(exprGroupsList), min = 0, max = 2 * pi)
    tt = matrix(tt, nrow = length(exprGroupsList), byrow = TRUE)
    times = t(apply(tt, 1, sort))
  }

  sm = foreach(cond = 1:length(exprGroupsList), .combine = rbind) %do% {
    sampleIds = ((cond - 1) * nSamples + 1):(cond * nSamples)
    sampleNames = sprintf(sprintf('sample_%%0%dd', floor(log10(nSamples)) + 1), sampleIds)
    data.table(sample = sampleNames, cond = cond,
               time = times[cond, ] * period / (2*pi))
  }

  emat = foreach(exprGroups = exprGroupsList, cond = 1:length(exprGroupsList),
                 .combine = cbind) %do% {
    simulateExprDataOneCond(exprGroups, times[cond, ], method)
  }

  colnames(emat) = sm$sample
  rownames(emat) = geneNames

  return(list(exprData = emat, sampleMetadata = sm, geneMetadata = gm))
}

#' @export
combineData = function(simData, geneNames) {
  d = data.table(t(simData$exprData[geneNames, ]), keep.rownames = TRUE)
  d = merge(simData$sampleMetadata, d, by.x = 'sample', by.y = 'rn')
  d = melt(d, measure.vars = geneNames, variable.name = 'gene',
           value.name = 'expr')
  d = merge(d, simData$geneMetadata, by = c('gene', 'cond'))
  return(d)
}

#' @export
getExpectedExpr = function(geneMetadata, time, period = 24) {
  d = geneMetadata[rep(1:.N, each = length(time))]
  d[, time := rep(..time, times = nrow(geneMetadata))]
  d[, mu := base + amp * rhyFunc[[1]]((time + phase) * 2 * pi / period)]
  return(d)
}
