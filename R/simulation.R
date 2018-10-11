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
    sd = exprGroups[group, sd]
    rhyFunc = exprGroups[group, rhyFunc][[1]]
    mu = amp * rhyFunc(timePoints + 2 * pi * phase) + base

    foreach::foreach(gene = 1:exprGroups[group, geneCount], .combine = rbind) %do% {
      if(method == 'gaussian'){
        stats::rnorm(length(mu), mu, sd = sd) }
      else {
        stats::rnbinom(length(mu), mu = mu, size = 1/sampleDispersion(2^mu)) }
    }
  }
}

setOneCondDefault = function(exprGroups, nGenes, randomTimepoints, nSamples,
                             defaultExpr, rhyFunc) {

  exprGroups = data.table::data.table(exprGroups)

  if(nrow(exprGroups) == 0) {
    stop('No rows in exprGroups. Cannot simulate genes.') }

  if(randomTimepoints && is.null(nSamples)) {
    stop('Number of random timepoint samples not specified.') }

  if(!'geneFrac' %in% colnames(exprGroups)) {
    exprGroups[, geneFrac := 1 / nrow(exprGroups)] }

  if(!'base' %in% colnames(exprGroups)) {
    exprGroups[, base := defaultExpr] }

  if(!'amp' %in% colnames(exprGroups)) {
    exprGroups[, amp := 0] }

  if(!'phase' %in% colnames(exprGroups)) {
    exprGroups[, phase := 0] }

  if(!'sd' %in% colnames(exprGroups)) {
    exprGroups[, sd := 1] }

  if(!'rhyFunc' %in% colnames(exprGroups)) {
    exprGroups[, rhyFunc := data.table::data.table(rhyFunc)] }

  if(any(exprGroups[, geneFrac] <= 0)) {
    stop('All groups in exprGroups must have geneFrac > 0.') }

  exprGroups[, group := 1:nrow(exprGroups)]

  # Compute a number of genes per group that sum to nGenes.
  exprGroups[, geneFrac := geneFrac / sum(geneFrac)]
  exprGroups[, geneCount := as.integer(geneFrac * nGenes)]
  if(sum(exprGroups[, geneCount]) != nGenes) {
    exprGroups[1L:(nGenes - sum(exprGroups[, geneCount])), geneCount := geneCount + 1]
  }

  return(exprGroups)
}


getMultipleCondSimulation = function(exprGroups, nGenes = 100, period = 24,
                                     interval = 4, nReps = 2, nSamples = NULL,
                                     randomTimepoints = FALSE, rhyFunc = sin,
                                     method = 'gaussian') {

  nCond = length(exprGroups)

  exprGroups = foreach::foreach(group = exprGroups, .combine = list) %do% {
    setOneCondDefault(group, nGenes, randomTimepoints, nSamples,
                      ifelse(method == 'gaussian', 0, 1), rhyFunc) }

  geneData = foreach::foreach(cond = 1:nCond, .combine = rbind) %do% {
    foreach::foreach(group = 1:nrow(exprGroups[[cond]]), .combine = rbind) %do% {
      foreach::foreach(gene = 1:exprGroups[[cond]][group, geneCount], .combine = rbind) %do% {
        data.table::data.table(base = exprGroups[[cond]][group, base],
                               amp = exprGroups[[cond]][group, amp],
                               phase = exprGroups[[cond]][group, phase],
                               sd = exprGroups[[cond]][group, sd],
                               rhyFunc = exprGroups[[cond]][group, rhyFunc],
                               group = group, cond = cond) } } }

  if(!randomTimepoints) {
    tt = (2 * pi / period) * interval * 0:(period %/% interval - (period %% interval == 0))
    tt = rep(tt, each = nReps)
    nSamples = length(tt)
    timePoints = foreach::foreach(cond = 1:nCond, .combine = rbind) %do% { tt }
  } else {
    timePoints = foreach::foreach(cond = 1:nCond, .combine = rbind) %do% {
      sort(stats::runif(nSamples, min = 0, max = 2 * pi)) }
  }

  sm = foreach::foreach(cond = 1:nCond, .combine = rbind) %do% {
    data.table::data.table(cond = cond, time = timePoints[cond, ] * period / (2*pi),
                           sample = paste('sample', ((cond-1)*nSamples+1):(cond*nSamples), sep = '_'))
  }
  sampleNames = sm$sample

  emat = foreach::foreach(group = length(exprGroups), .combine = cbind) %do% {
    getSingleCondSim(exprGroups[[group]], timePoints[group, ], method)
  }

  return(list(emat = emat, sm = sm, gm = geneData))
}