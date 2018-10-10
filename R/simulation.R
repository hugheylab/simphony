#' @importFrom data.table ":="
#' @importFrom foreach "%do%"
globalVariables(c('geneFrac', 'meanExpr', 'dExpr', 'meanPhase', 'group',
                  'geneCount', 'ii', 'dAmp', 'dPhase', 'meanAmp', 'meanSd',
                  'dSd', 'mean'))


sampleDispersion = function(x) {
  return(rep(3, length(x)))
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
        stats::rnbinom(length(mu), mu = mu, size = mu + sampleDispersion(mu) * mu ^ 2) }
    }
  }
}

splitTwoCondExprGroups = function(exprGroups) {
  exprGroups1 = data.table::data.table(
    base = exprGroups[, meanBase] + exprGroups[, dBase],
    phase = exprGroups[, meanPhase] + exprGroups[, dPhase],
    amp = exprGroups[, meanAmp] + exprGroups[, dAmp],
    sd = exprGroups[, meanSd] + exprGroups[, dSd],
    rhyFunc = exprGroups[, rhyFunc],
    group = exprGroups[, group],
    geneCount = exprGroups[, geneCount]
  )
  exprGroups2 = data.table::data.table(
    base = exprGroups[, meanBase] - exprGroups[, dBase],
    phase = exprGroups[, meanPhase] - exprGroups[, dPhase],
    amp = exprGroups[, meanAmp] - exprGroups[, dAmp],
    sd = exprGroups[, meanSd] - exprGroups[, dSd],
    rhyFunc = exprGroups[, rhyFunc],
    group = exprGroups[, group],
    geneCount = exprGroups[, geneCount]
  )

  return(list(exprGroups1 = exprGroups1, exprGroups2 = exprGroups2))
}

setOneCondDefault = function(exprGroups, nGenes, randomTimepoints, nSamples,
                             defaultExpr) {

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
    exprGroups[, rhyFunc := data.table::data.table(sin)] }

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

setTwoCondDefault = function(exprGroups, nGenes, randomTimepoints, nSamples,
                             defaultExpr) {

  exprGroups = data.table::data.table(exprGroups)

  if(nrow(exprGroups) == 0) {
    stop('No rows in exprGroups. Cannot simulate genes.') }

  if(randomTimepoints && is.null(nSamples)) {
    stop('Number of random timepoint samples not specified.') }

  if(!'geneFrac' %in% colnames(exprGroups)) {
    exprGroups[, geneFrac := 1 / nrow(exprGroups)] }

  if(!'meanExpr' %in% colnames(exprGroups)) {
    exprGroups[, meanExpr := defaultExpr] }

  if(!'dExpr' %in% colnames(exprGroups)) {
    exprGroups[, dExpr := 0] }

  if(!'meanAmp' %in% colnames(exprGroups)) {
    exprGroups[, meanAmp := 0] }

  if(!'dAmp' %in% colnames(exprGroups)) {
    exprGroups[, dAmp := 0] }

  if(!'meanPhase' %in% colnames(exprGroups)) {
    exprGroups[, meanPhase := 0] }

  if(!'dPhase' %in% colnames(exprGroups)) {
    exprGroups[, dPhase := 0] }

  if(!'meanSd' %in% colnames(exprGroups)) {
    exprGroups[, meanSd := 1] }

  if(!'dSd' %in% colnames(exprGroups)) {
    exprGroups[, dSd := 0] }

  if(!'rhyFunc' %in% colnames(exprGroups)) {
    exprGroups[, rhyFunc := data.table::data.table(sin)] }

  if(any(exprGroups[, meanSd] - exprGroups[, dSd] / 2 < 0)) {
    stop('Groups cannot have negative standard deviation of sample error.') }

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

getOneMetadata = function(exprGroups, randomTimepoints, period, interval, nReps,
                       nSamples, nTotalGenes) {

  if(!randomTimepoints) {
    timePoints = (2 * pi / period) * interval * 0:(period %/% interval - (period %% interval == 0))
    timePoints = rep(timePoints, each = nReps)
    nSamples = length(timePoints)
  } else {
    timePoints = sort(stats::runif(nSamples, min = 0, max = 2 * pi)) }

  sampleNames = paste('sample', 1:nSamples, sep = '_')
  geneNames = paste('gene', 1:nTotalGenes, sep = '_')

  sampleMetadata = data.table::data.table(cond = 1,
                                          time = timePoints * period / (2 * pi),
                                          sample = sampleNames)

  geneMetadata = data.table::data.table(gene = geneNames,
                                        group = rep(1:nrow(exprGroups),
                                                    times = exprGroups[, geneCount]))

  return(list(timePoints = timePoints,
              sampleNames = sampleNames, geneNames = geneNames,
              sampleMetadata = sampleMetadata, geneMetadata = geneMetadata,
              nSamples = nSamples))
}

getTwoMetadata = function(exprGroups, randomTimepoints, period, interval, nReps,
                       nSamples, nTotalGenes) {

  if(!randomTimepoints) {
    timePoints = (2 * pi / period) * interval * 0:(period %/% interval - (period %% interval == 0))
    timePoints = rep(timePoints, each = nReps)
    nSamples = length(timePoints)
    timePoints = rep(timePoints, 2)
  } else {
    timePoints = foreach::foreach(cond = 1:2, .combine = c) %do% {
      sort(stats::runif(nSamples, min = 0, max = 2 * pi)) }
  }

  sampleNames = paste('sample', 1:(2*nSamples), sep = '_')
  geneNames = paste('gene', 1:nTotalGenes, sep = '_')

  sampleMetadata = data.table::data.table(cond = rep(1:2, each = nSamples),
                                          time = timePoints * period / (2 * pi),
                                          sample = sampleNames)

  geneMetadata = data.table::data.table(gene = geneNames,
                                        group = rep(1:nrow(exprGroups),
                                                    times = exprGroups[, geneCount]))

  return(list(timePoints1 = timePoints[1:(length(timePoints) / 2)],
              timePoints2 = timePoints[(1 + length(timePoints) / 2):length(timePoints)],
              sampleNames = sampleNames, geneNames = geneNames,
              sampleMetadata = sampleMetadata, geneMetadata = geneMetadata,
              nSamples = nSamples))
}

#' @export
getSingleCond = function(exprGroups, nGenes = 100, period = 24, interval = 4,
                         nReps = 2, nSims = 1, nSamples = NULL,
                         randomTimepoints = FALSE, method = 'gaussian') {

  exprGroups = setOneCondDefault(exprGroups, nGenes, randomTimepoints, nSamples,
                                 ifelse(method == 'gaussian', 0, 1))

  metadata = getOneMetadata(exprGroups, randomTimepoints, period, interval, nReps,
                            nSamples, nGenes * nSims)

  emat = foreach::foreach(sim = 1L:nSims, .combine = rbind) %do% {
    getSingleCondSim(exprGroups, metadata$timePoints, method)
  }

  colnames(emat) = metadata$sampleNames
  rownames(emat) = metadata$geneNames

  results = list(emat = emat, sm = metadata$sampleMetadata,
                 gm = metadata$geneMetadata, exprGroups = exprGroups)
  return(results)
}

getTwoCond = function(exprGroups, nGenes = 100, period = 24, interval = 4,
                      nReps = 2, nSims = 1, nSamples = NULL,
                      randomTimepoints = FALSE, method = 'gaussian') {

  exprGroups = setTwoCondDefault(exprGroups, nGenes, randomTimepoints, nSamples,
                                 ifelse(method == 'gaussian', 0, 1))

  metadata = getTwoMetadata(exprGroups, randomTimepoints, period, interval, nReps,
                            nSamples, nGenes * nSims)

  groups = splitTwoCondExprGroups(exprGroups)

  emat = foreach::foreach(sim = 1L:nSims, .combine = rbind) %do% {
    cbind(getSingleCondSim(groups$exprGroups1, metadata$timePoints1, method),
          getSingleCondSim(groups$exprGroups2, metadata$timePoints2, method))
  }

  colnames(emat) = metadata$sampleNames
  rownames(emat) = metadata$geneNames

  results = list(emat = emat, sm = metadata$sampleMetadata,
                 gm = metadata$geneMetadata, exprGroups = exprGroups)
  return(results)
}