#' @importFrom data.table ":="
#' @importFrom foreach "%do%"
globalVariables(c('geneFrac', 'meanExpr', 'dExpr', 'meanPhase', 'group',
                  'geneCount', 'ii', 'dAmp', 'dPhase', 'meanAmp', 'meanSd',
                  'dSd', 'mean'))


sampleDispersion = function(x) {
  return(rep(3, length(x)))
}

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
                               mu = rep(mu, exprGroups[group, geneCount]),
                               size = 1/sampleDispersion(2^mu)) }
    return(matrix(samples, nrow = exprGroups[group, geneCount]))
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

  if(!'base' %in% colnames(exprGroups)) {
    exprGroups[, base := 7] }

  if(!'amp' %in% colnames(exprGroups)) {
    exprGroups[, amp := 0] }

  if(!'phase' %in% colnames(exprGroups)) {
    exprGroups[, phase := 0] }

  if(!'rhyFunc' %in% colnames(exprGroups)) {
    exprGroups[, rhyFunc := data.table::data.table(rhyFunc)] }

  if(method == 'rnbinom') {
    if(!'dispersionFunc' %in% colnames(exprGroups)) {
      exprGroups[, dispersionFunc := data.table::data.table(sampleDispersion)] } }
  else {
    if(!'sd' %in% colnames(exprGroups)) {
      exprGroups[, sd := 1] } }

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


getMultipleCondSimulation = function(exprGroupsList, nGenes = 100, period = 24,
                                     interval = 4, nReps = 2, nSamples = NULL,
                                     randomTimepoints = FALSE, rhyFunc = sin,
                                     method = 'gaussian') {

  exprGroupsList = foreach::foreach(exprGroups = exprGroupsList, .combine = list) %do% {
    setOneCondDefault(exprGroups, nGenes, randomTimepoints, nSamples, rhyFunc, method) }

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
    timePoints = foreach::foreach(cond = 1:length(exprGroupsList), .combine = rbind) %do% {
      sort(stats::runif(nSamples, min = 0, max = 2 * pi)) }
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