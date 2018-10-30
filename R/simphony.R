#' @importFrom data.table data.table ":="
#' @importFrom foreach foreach "%do%"
globalVariables(c('base', 'amp', 'phase', 'group', 'rhyFunc', 'sd', 'cond',
                  'dispFunc', 'exprGroups', 'numGenes', 'fracGenes', 'time',
                  '..cond', '..time', '.N', '.dummy', 'gene', 'mu'))


#' Generate list of two expression groups from a combined differential exprGroup
#'
#'
#' @param diffExprGroups is the differential exprGroups data.frame to convert
#' into a list of single-condition exprGroups data.frames.
#' @param checkValid stuff
#' @examples
#' dGroups = data.frame(meanAmp = c(1, 1, 1, 1), dAmp = c(1, 1, 2, 2),
#'                      meanPhase = c(0, 0, 0, 0), dPhase = c(0, 3, 0, 3))
#' exprGroups = splitDiffExprGroups(dGroups)
#' @export
splitDiffExprGroups = function(diffExprGroups, checkValid = TRUE) {
  dGroups = data.table(diffExprGroups)

  capCols = c('Base', 'Amp', 'Phase', 'Sd')
  cols = tolower(capCols)
  meanCols = paste0('mean', capCols)
  dCols = paste0('d', capCols)

  d1 = data.table(.dummy = rep(1, nrow(dGroups)))
  d2 = data.table(.dummy = rep(1, nrow(dGroups)))
  for (ii in 1:length(cols)) {
    if (all(c(meanCols[ii], dCols[ii]) %in% colnames(dGroups))) {
      d1[[cols[ii]]] = dGroups[[meanCols[ii]]] - 0.5 * dGroups[[dCols[ii]]]
      d2[[cols[ii]]] = dGroups[[meanCols[ii]]] + 0.5 * dGroups[[dCols[ii]]]
    }
  }
  d1[, .dummy := NULL]
  d2[, .dummy := NULL]

  heldbackCols = setdiff(colnames(dGroups), c(meanCols, dCols))
  if (length(heldbackCols) > 0) {
    dHeldback = dGroups[, heldbackCols, with = FALSE]
    d1 = cbind(d1, dHeldback)
    d2 = cbind(d2, dHeldback)
  }

  if (checkValid) {
    idx = rep(TRUE, nrow(d1))
    if ('amp' %in% colnames(d1)) {
      idx = idx & (d1$amp >= 0) & (d2$amp >= 0)
    }
    if ('sd' %in% colnames(d1)) {
      idx = idx & (d1$sd > 0) & (d2$sd > 0)
    }
    d1 = d1[idx]
    d2 = d2[idx]
  }

  return(list(d1, d2))
}

#' Generate simulated gene expression time courses
#'
#' @param exprGroupsList is a list of data.frame or data.table objects with the
#"   following optional columns:
#'   \itemize{
#'     \item{fracGenes}: {Fraction of all simulated genes which fall into this
#'                       group. Defaults to 1/nrow(exprGroups) if not supplied.}
#'     \item{base}: {The baseline expression for this group. Defaults
#'                       to 0 if not supplied.}
#'     \item{amp}: {The amplitude of the rhythmic component of
#'                      expression for this group. Defaults to 0 if not
#'                      supplied.}
#'     \item{phase}: {The phase of the rhythmic component of expression
#'                        for this group. Defaults to 0 if not supplied.}
#'     \item{sd}: {The standard deviation of the sample error for this
#'                     group. Defaults to 1 if not supplied.}
#'     \item{rhyFunc}: {The function used to generate the rhythmic component of
#'                      this group's gene expression. rhyFunc must have a period
#'                      of 2*pi. Defaults to sine if not supplied.}
#'   }
#' @param fracGenes Fraction of all simulated genes which fall into each group.
#'   Must have length of number of groups in each exprGroups object. Defaults to
#'   1/(number of groups).
#' @param nGenes is the integer number of total genes to simulate.
#' @param period is the integer number of hours in one rhythmic cycle.
#' @param timepointsType stuff
#' @param interval is the integer number of hours between simulated time points.
#' @param nReps is the integer number of replicates per time point.
#' @param timepoints stuff
#' @param nSamplesPerCond stuff
#' @param rhyFunc is the rhythmic function to set for exprGroups missing a
#'   rhythmic function. Defaults to sin if not supplied.
#' @param method is the data generation method to use. Must be either 'gaussian'
#'   or 'negbinom'.
#' @export
simulateExprData = function(exprGroupsList, fracGenes = NULL, nGenes = 10,
                            period = 24, timepointsType = 'auto', interval = 4,
                            nReps = 2, timepoints = NULL, nSamplesPerCond = NULL,
                            rhyFunc = sin, method = 'gaussian') {
  if (!method %in% c('gaussian', 'negbinom')) {
    stop("method must be 'gaussian' or 'negbinom'.")}

  if (is.data.frame(exprGroupsList)) {
    exprGroupsList = list(exprGroupsList)}
  if (length(unique(sapply(exprGroupsList, nrow))) != 1) {
    stop('Each exprGroups data.frame must have the same number of rows.')}

  exprGroupsList = foreach(exprGroups = exprGroupsList) %do% {
    setDefaultExprGroups(exprGroups, nGenes, rhyFunc, method)}
  nCond = length(exprGroupsList)

  if ('fracGenes' %in% colnames(exprGroupsList[[1]])) {
    fracGenes = exprGroupsList[[1]]$fracGenes }
  if (is.null(fracGenes)) {
    fracGenes = rep(1/nrow(exprGroupsList[[1]]), nrow(exprGroupsList[[1]]))}
  numGenes = as.integer(fracGenes * nGenes)
  if (sum(numGenes) != nGenes) {
    numGenes[1L:(nGenes - sum(numGenes))] = numGenes[1L:(nGenes - sum(numGenes))] + 1L}

  if (any(numGenes) == 0) {
    stop(paste(c('At least one group has no genes. Increase nGenes,',
                 'reduce the number of groups, or change fracGenes.'),
               collapse = ' '))}

  times = getTimes(timepointsType, interval, nReps, timepoints,
                   nSamplesPerCond, nCond, period)

  geneNames = sprintf(sprintf('gene_%%0%dd', floor(log10(nGenes)) + 1), 1:nGenes)

  gm = foreach(exprGroups = exprGroupsList, cond = 1:nCond, .combine = rbind) %do% {
    gmNow = exprGroups[rep(1:.N, times = numGenes)]
    gmNow[, cond := ..cond]
    gmNow[, gene := geneNames]
    data.table::setcolorder(gmNow, c('cond', 'group', 'gene'))
    gmNow
  }

  nSamples = prod(dim(times))
  nSamplesPerCond = ncol(times)
  sm = foreach(cond = 1:nCond, .combine = rbind) %do% {
    sampleIds = ((cond - 1) * nSamplesPerCond + 1):(nSamplesPerCond * cond)
    sampleNames = sprintf(sprintf('sample_%%0%dd', floor(log10(nSamples)) + 1),
                          sampleIds)
    data.table(sample = sampleNames, cond = cond,
               time = times[cond, ] * period / (2*pi))
  }

  emat = foreach(exprGroups = exprGroupsList, cond = 1:nCond, .combine = cbind) %do% {
    simulateExprDataOneCond(exprGroups, numGenes, times[cond, ], method)
  }

  colnames(emat) = sm$sample
  rownames(emat) = geneNames

  return(list(exprData = emat, sampleMetadata = sm, geneMetadata = gm))
}

#' @export
combineData = function(simData, geneNames) {
  d = data.table(t(simData$exprData[geneNames, ]), keep.rownames = TRUE)
  d = merge(simData$sampleMetadata, d, by.x = 'sample', by.y = 'rn')
  d = data.table::melt(d, measure.vars = geneNames, variable.name = 'gene',
                       value.name = 'expr')
  d = merge(d, simData$geneMetadata, by = c('gene', 'cond'))
  return(d)
}

#' @export
getExpectedExpr = function(geneMetadata, times, period = 24) {
  d = data.table(geneMetadata)[rep(1:.N, each = length(times))]
  d[, time := rep(times, times = nrow(geneMetadata))]
  d[, mu := base + amp * rhyFunc[[1]]((time + phase) * 2 * pi / period), by = 1:nrow(d)]
  return(data.table::copy(d))
}
