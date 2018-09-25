#' @importFrom data.table ":="
#' @importFrom foreach "%do%"
globalVariables(c('geneFrac', 'meanExpr', 'dExpr', 'meanPhase', 'group',
                  'geneCount', 'ii', 'dAmp', 'dPhase', 'meanAmp', 'meanSd',
                  'dSd', 'mean'))

checkExprGroups = function(exprGroups, nGenes, randomTimepoints, nSamples) {

  exprGroups = data.table::data.table(exprGroups)

  if(nrow(exprGroups) == 0) {
    stop('No rows in exprGroups. Cannot simulate genes.') }

  if(randomTimepoints && is.null(nSamples)) {
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

  if(!'meanSd' %in% colnames(exprGroups)) {
    exprGroups[, meanSd := 1] }

  if(!'dSd' %in% colnames(exprGroups)) {
    exprGroups[, dSd := 0] }

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
#' @param useNegBinom is a boolean determining whether to simulate read counts
#'   using a negative binomial distribution. Defaults to FALSE.
#' @param sizeNegBinom is a positive number relating the mean and dispersion factor
#'   of negative binomial with dipersion = sizeNegBinom * mean. Defaults to 0.333.
#' @param logNegBinom is a boolean determining whether to log transform read counts
#'   into expression levels. Defaults to FALSE.
#' @return The simulated expression matrix emat, simulated sample metadata sm,
#'   metadata describing the expression group each gene came from gm, and an
#'   updated exprGroups including missing property columns and an index; group.
#' @details exprGroups must be a data.frame or data.table object, with the
#'   following optional columns:
#'   \itemize{
#'     \item{geneFrac}: {The fraction of all simulated genes which fall into
#'                       this group. Defaults to 1/nrow(exprGroups) if not
#'                       supplied.}
#'     \item{meanExpr}: {The mean baseline expression for this group. Defaults 
#'                       to 0 if not supplied.}
#'     \item{dExpr}: {The difference in baseline expression across conditions
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
#'   }
#' @examples
#'   exprGroups = data.table::data.table(meanAmp = c(1,2,1,2), dAmp = c(1,1,2,2))
#'   gse = getSimulatedExpr(exprGroups, nGenes = 10000, randomTimepoints = TRUE,
#'                          nSamples = 10)
#' @export
getSimulatedExpr = function(exprGroups, nGenes = 100, period = 24, interval = 4,
                            nReps = 2, errSd = 1, nSims = 1, nSamples = NULL,
                            randomTimepoints = FALSE, rhyFunc = sin, useNegBinom = FALSE,
                            sizeNegBinom = 0.333, logNegBinom = FALSE) {

  exprGroups = checkExprGroups(exprGroups, nGenes, randomTimepoints, nSamples)

  if(useNegBinom == FALSE && logNegBinom == TRUE) {
  stop('Log transformed read counts are only supported with useNegBinom = TRUE') }

  if(!randomTimepoints) {
    timePoints = (2 * pi / period) * interval * 0:(period %/% interval - (period %% interval == 0))
    timePoints = rep(timePoints, each = nReps)
    nSamples = length(timePoints)
    timePoints = rep(timePoints, 2) # Number of conditions = 2
  } else {
    timePoints = sort(stats::runif(nSamples, min = 0, max = 2 * pi))
    timePoints = c(timePoints, sort(stats::runif(nSamples, min = 0, max = 2 * pi)))
  }

  timePoints1 = timePoints[1:(length(timePoints)/2)]
  timePoints2 = timePoints[(length(timePoints)/2 + 1):length(timePoints)]

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

      sd1 = exprGroups[ii, meanSd] + exprGroups[ii, dSd] / 2
      sd2 = exprGroups[ii, meanSd] - exprGroups[ii, dSd] / 2

      if(!'rhyFunc' %in% colnames(exprGroups)) {
        rhyFunc = sin }
      else {
        rhyFunc = exprGroups[ii, rhyFunc][[1]] }

      # Compute the expression matrix for this exprGroup
      foreach::foreach(jj = 1L:exprGroups[ii, geneCount], .combine = rbind) %do% {
        timeCourse1 = amp1 * rhyFunc(timePoints1 + 2 * pi * phase1 / period) + expr1
        timeCourse2 = amp2 * rhyFunc(timePoints2 + 2 * pi * phase2 / period) + expr2

        if(!useNegBinom){
          timeCourse1 = timeCourse1 + stats::rnorm(nSamples, sd = sd1)
          timeCourse2 = timeCourse2 + stats::rnorm(nSamples, sd = sd2)
        } else {
          #Default size parameter is 0.333*mean hence variance = 4 * mean (Polyester based)
          timeCourse1 = foreach::foreach(mean = timeCourse1, .combine='c') %do% {
            stats::rnbinom(1, mu = 2^mean - 1, size = sizeNegBinom * (2^mean - 1))
          }
          timeCourse2 = foreach::foreach(mean = timeCourse2, .combine='c') %do% {
            stats::rnbinom(1, mu = 2^mean - 1, size = sizeNegBinom * (2^mean - 1))
          }

          if(logNegBinom){
            timeCourse1 = log2(timeCourse1 + 1) #Prevent taking the log of zero
            timeCourse2 = log2(timeCourse2 + 1)
          }
        }
        c(timeCourse1, timeCourse2)
      }
    }
  }

  colnames(emat) = sampleNames
  rownames(emat) = geneNames

  results = list(emat = emat, sm = sampleMetadata, gm = geneMetadata,
                 exprGroups = exprGroups)
  return(results)
}
