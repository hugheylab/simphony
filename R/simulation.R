#' @importFrom data.table ":="
#' @importFrom foreach "%do%"
globalVariables(c('geneFrac', 'meanExpr', 'dExpr', 'meanPhase', 'group',
                  'geneCount', 'ii', 'dAmp', 'dPhase', 'meanAmp'))

#' Generate simulated gene expresion time courses.
#'
#' @param exprGroups is a dataframe of metadata describing the properties of
#'   rhythmic or differentially rhythmic genes.
#' @param nGenes is the integer number of total genes to simulate.
#' @param period is the integer number of hours in one rhythmic cycle.
#' @param interval is the integer number of hours between simulated time points.
#' @param nReps is the integer number of replicates per time point.
#' @param errSd is the standard deviation of the Gaussian sample error.
#' @param nSims is the integer number of simulations to generate.
#' @param randomTimepoints is a boolean determining whether to simulate an
#'   experiment with random sample times. Defaults to FALSE.
#' @param nSamples is the integer number of time points to sample, if
#'   randomTimepoints is enabled. This must be supplied if randomTimepoints is
#'   TRUE.
#' @param rhyFunc is the function defining the rhythmic component of the
#'   simulated gene expression time course. Defaults to sin.
#' @export
getSimulatedExpr = function(exprGroups, nGenes = 100, period = 24,
                                    interval = 4, nReps = 2, errSd = 1,
                                    nSims = 1, randomTimepoints = FALSE,
                                    nSamples = NULL, rhyFunc = sin) {

  exprGroups = data.table::data.table(exprGroups)

  if(nrow(exprGroups) == 0) {
    stop('No rows in exprGroups. Cannot simulate genes.') }

  if(randomTimepoints & is.null(nSamples)) {
    stop('Number of random timepoint samples not specified.') }

  if(!'geneFrac' %in% colnames(exprGroups)) {
    exprGroups[, geneFrac := 1 / nrow(exprGroups)] }

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

  exprGroups[, group := 1:nrow(exprGroups)]

  # Compute a number of genes per group that sum to nGenes.
  exprGroups[, geneFrac := geneFrac / sum(geneFrac)]
  exprGroups[, geneCount := as.integer(geneFrac * nGenes)]
  if(sum(exprGroups[, geneCount]) != nGenes) {
    exprGroups[1L:(nGenes - sum(exprGroups[, geneCount])), geneCount := geneCount + 1]
  }

  if(!randomTimepoints) {
    timePoints = (2 * pi / period) * interval * 0:(period %/% interval - (period %% interval == 0))
    timePoints = rep(timePoints, each = nReps)
    nSamples = length(timePoints)
    timePoints = rep(timePoints, 2)
  } else {
    timePoints = sort(stats::runif(nSamples, min = 0, max = 2 * pi))
    timePoints = c(timePoints, sort(stats::runif(nSamples, min = 0, max = 2 * pi)))
  }

  sampleNames = paste('sample', 1:(nSamples * 2), sep = '_')
  geneNames = paste('gene', 1:(nGenes * nSims), sep = '_')

  sampleMetadata = data.table::data.table(cond = rep(1:2, each = nSamples),
                                          time = timePoints * period / (2 * pi),
                                          sample = sampleNames)

  geneMetadata = data.table::data.table(gene = geneNames,
                                        group = rep(1:nrow(exprGroups),
                                                    times = exprGroups[, geneCount]))
  
  emat = foreach::foreach(sim = 1L:nSims, .combine = rbind) %do% {
    foreach::foreach(ii = 1L:nrow(exprGroups), .combine = rbind) %do% {
      expr1 = exprGroups[ii, meanExpr] + exprGroups[ii, dExpr] / 2
      expr2 = exprGroups[ii, meanExpr] - exprGroups[ii, dExpr] / 2

      amp1 = exprGroups[ii, meanAmp] + exprGroups[ii, dAmp] / 2
      amp2 = exprGroups[ii, meanAmp] - exprGroups[ii, dAmp] / 2

      phase1 = exprGroups[ii, meanPhase] + exprGroups[ii, dPhase] / 2
      phase2 = exprGroups[ii, meanPhase] - exprGroups[ii, dPhase] / 2

      # Compute the expression matrix for this exprGroup
      foreach::foreach(jj = 1L:exprGroups[ii, geneCount], .combine = rbind) %do% {
        timeCourse1 = amp1 * rhyFunc(timePoints[1:(length(timePoints)/2)] + 2 * pi * phase1 / period) + expr1
        timeCourse2 = amp2 * rhyFunc(timePoints[(length(timePoints)/2 + 1):length(timePoints)] + 2 * pi * phase2 / period) + expr2
        c(timeCourse1, timeCourse2)
      }
    }
  }

  error = matrix(stats::rnorm(nGenes * nSamples, sd = errSd), nGenes * nSims, nSamples * 2)
  emat = emat + error

  colnames(emat) = sampleNames
  rownames(emat) = geneNames

  results = list(emat = emat, sm = sampleMetadata, gm = geneMetadata,
                 exprGroups = exprGroups)
  return(results)
}
