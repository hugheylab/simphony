library(data.table)
library(foreach)

getSimulatedExpr = function(exprGroups, nGenes = 100, period = 24, interval = 4,
                            nReps = 2, errSd = 1, nSims = 1,
                            randomTimepoints = FALSE, nSamples = 0) {

  if(nrow(exprGroups) == 0) {
    stop('No rows in exprGroups. Cannot simulate genes.') }

  if(randomTimepoints & nSamples == 0) {
    stop('Number of random timepoint samples not specified (nSamples).') }

  if(!('geneFrac' %in% colnames(exprGroups)) & !('geneCount' %in% colnames(exprGroups))) {
    exprGroups[, geneFrac := 1/nrow(exprGroups)] }

  if(!'meanExpr' %in% colnames(exprGroups)) {
    exprGroups[, meanExpr := 0] }

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

  exprGroups[, index := 1:nrow(exprGroups)]

  # Compute a number of genes per group that sum to nGenes.
  if(!'geneCount' %in% colnames(exprGroups)) {
    exprGroups[, geneCount := as.integer(geneFrac * nGenes)]
    if(sum(exprGroups[, geneCount]) != nGenes) {
      exprGroups[1L:(nGenes - sum(exprGroups[, geneCount])), geneCount := geneCount + 1]
    }
  }

  if(!randomTimepoints) {
    timePoints = (2 * pi / period) * interval * 0:(period / interval - 1)
    nSamples = nReps * period %/% interval
  } else {
    timePoints = sort(runif(nSamples, min = 0, max = 2 * pi))
    nSamples = nSamples * nReps
  }
  timePoints = rep(timePoints, each = nReps)

  sampleNames = paste('sample', 1:(2 * nSamples), sep = '_')
  geneNames = paste('gene', 1:(nGenes * nSims), sep = '_')

  sampleMetadata = data.table::data.table(cond = rep(1:2, each = nSamples),
                                          time = rep(timePoints * period / (2 * pi), 2),
                                          sample = sampleNames)

  geneMetadata = data.table::data.table(gene = geneNames,
                                        index = rep(1:nrow(exprGroups),
                                                    times = exprGroups[, geneCount]))
  

  emat = foreach(sim = 1L:nSims, .combine = rbind) %do% {
    foreach(ii = 1L:nrow(exprGroups), .combine = rbind) %do% {
      expr1 = exprGroups[ii, meanExpr] + exprGroups[ii, dExpr] / 2
      expr2 = exprGroups[ii, meanExpr] - exprGroups[ii, dExpr] / 2

      amp1 = exprGroups[ii, meanAmp] + exprGroups[ii, dAmp] / 2
      amp2 = exprGroups[ii, meanAmp] - exprGroups[ii, dAmp] / 2

      phase1 = exprGroups[ii, meanPhase] + exprGroups[ii, dPhase] / 2
      phase2 = exprGroups[ii, meanPhase] - exprGroups[ii, dPhase] / 2

      # Compute the expression matrix for this exprGroup
      ematNow = foreach(jj = 1L:exprGroups[ii, geneCount], .combine = rbind) %do% {
        timeCourse1 = amp1 * sin(timePoints + 2 * pi * phase1 / period) + expr1
        timeCourse2 = amp2 * sin(timePoints + 2 * pi * phase2 / period) + expr2
        c(timeCourse1, timeCourse2)
      }
    }
  }

  error = matrix(rnorm(nGenes * nSamples, sd = errSd), nGenes, nSamples * 2)
  emat = emat + error

  colnames(emat) = sampleNames
  rownames(emat) = geneNames

  results = list(emat = emat,
                 sm = sampleMetadata,
                 gm = geneMetadata,
                 exprGroups = exprGroups)
  return(results)
}
