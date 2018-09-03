#' @importFrom data.table ":="
#' @importFrom foreach "%do%"
#' @importFrom stats rnorm time
globalVariables(c('cond', 'dAmp', 'dPhase', 'foreach', 'geneIndex', 'meanAmp',
                  'simGroup'))


getRhythmicExpr = function(timeSteps, phase = 0, amplitude = 1,
                           verticalShift = 0, rhyFunc = sin) {

  exprs = amplitude * rhyFunc(timeSteps + phase) + verticalShift
  return(exprs)
}


getNonRhythmicExpr = function(nSamples, nCond = 1, nRow = 1,
                              verticalShift = 0) {
  exprs = matrix(rep.int(rep.int(verticalShift, nSamples), nCond),
                 nRow, nSamples * nCond)
  return(exprs)
}


getExprAmplitudes = function(nSamples, discrete = FALSE, amps = c()) {
  if(discrete != FALSE) {
    return(rep(amps, nSamples %/% length(amps) + 1))
  }
  return(pmax(rnorm(nSamples, sd = 0.5, mean = 3), 0.5))
}


getPreparedVariables = function(conditions, interval, period, nCond, nReps,
                                nRhyGenes, nGenes, nSims, nDrGenes, sampleNames,
                                geneNames) {

  sampleMetadata = data.table::data.table(cond = conditions,
                              time = rep(interval * 0:(period / interval - 1),
                                         each = nCond * nReps))
  sampleMetadata = sampleMetadata[order(cond, time)]
  sampleMetadata[, sample := sampleNames]

  featureMetadata = data.table::data.table(gene = geneNames,
                               rhy = rep(c(rep(1, nRhyGenes),
                                           rep(0, nGenes - nRhyGenes)), nSims),
                               dRhy = rep(c(rep(0, nRhyGenes - nDrGenes),
                                            rep(1, nDrGenes),
                                            rep(0, nGenes - nRhyGenes)), nSims),
                               geneIndex = rep(rep(as.integer(NA), nGenes),
                                               nSims),
                               simIndex = rep(1:nSims, each = nGenes))

  timeSteps = rep((2 * pi / period) * interval * 0:(period / interval - 1),
                  each = nReps)

  results = list(sampleMetadata = sampleMetadata,
                 featureMetadata = featureMetadata,
                 timeSteps = timeSteps)
  return(results)
}


getRhyExprMatrix = function(nGenes, simGroup, nRhyOnlyGenes, nDrGenes, nSamples,
                            discreteAmps, exprGroups, timeSteps, nCond,
                            geneNames, conditions, nSims, featureMetadata) {

  nRhyGenesPerSim = nRhyOnlyGenes + nDrGenes
  if(nRhyOnlyGenes > 0) {
    nDrGroups = nrow(exprGroups[dAmp != 0 | dPhase != 0, ])
    nRhyGroups = nrow(exprGroups[dAmp == 0 & dPhase == 0, ])

    rhyExprs = matrix(ncol = nSamples, nrow = nRhyOnlyGenes)
    amplitudes = getExprAmplitudes(nRhyOnlyGenes, discrete = discreteAmps,
                                   amps = t(exprGroups[dAmp == 0 & dPhase == 0, meanAmp]))
    for(ii in 1L:nRhyOnlyGenes) {
      condExprs = getRhythmicExpr(timeSteps, amplitude = amplitudes[ii])
      rhyExprs[ii, ] = rep(condExprs, nCond)
      if(discreteAmps) {
        groupIndex = nDrGroups + (ii - 1) %% nRhyGroups + 1
        data.table::set(featureMetadata, i = ii + (simGroup - 1L) * nGenes,
            j = 'geneIndex',
            value = exprGroups[groupIndex, geneIndex])
      }
    }
  } else {
    rhyExprs = matrix(ncol = nSamples, nrow = 0)
  }

  results = list(rhyExprs = rhyExprs,
                 featureMetadata = featureMetadata)
  return(results)
}


getDrExprMatrix = function(simGroup, nDrGenes, exprGroups, nSamples,
                           discreteAmps, period, timeSteps, geneNames,
                           nRhyOnlyGenes, nSims, featureMetadata,
                           nSamplesPerCond, conditions, nGenes) {

  nRhyGenesPerSim = nRhyOnlyGenes + nDrGenes
  if(nDrGenes > 0) {
    nGenesPerGroup = nDrGenes / as.integer(nrow(exprGroups[dAmp != 0 | dPhase != 0, ]))
    dRhyExprs = matrix(nrow = nDrGenes, ncol = nSamples)

    if(discreteAmps) {
      for(row in 1L:nrow(exprGroups[dAmp != 0 | dPhase != 0, ])) {
        condAmps = c(exprGroups[row, meanAmp] -
                     exprGroups[row, dAmp] / 2,
                     exprGroups[row, meanAmp] +
                     exprGroups[row, dAmp] / 2)

        phaseDelta = (2 * pi / period) * exprGroups[row, dPhase] / 2
        condPhases = c(exprGroups[row, meanPhase] + phaseDelta,
                       exprGroups[row, meanPhase] - phaseDelta)

        baseExpr = c(getRhythmicExpr(timeSteps, phase = condPhases[1],
                                     amplitude = condAmps[1]),
                     getRhythmicExpr(timeSteps, phase = condPhases[2],
                                     amplitude = condAmps[2]))

        for(ii in 1L:nGenesPerGroup) {
          rowIndex = (row - 1L) * nGenesPerGroup + ii
          dRhyExprs[rowIndex, ] = baseExpr
          gene = geneNames[rowIndex + nRhyOnlyGenes + (simGroup - 1) * nSims]

          data.table::set(featureMetadata,
              i = rowIndex + nRhyOnlyGenes + (simGroup - 1L) * nGenes,
              j = 'geneIndex', value = exprGroups[row, geneIndex])
        }
      }
    } else {
      amplitudes = getExprAmplitudes(nDrGenes)

      for(ii in 1L:nDrGenes) {
        exprs = getRhythmicExpr(timeSteps, amplitudes[ii])
        exprs = c(exprs, getNonRhythmicExpr(nSamplesPerCond))
        dRhyExprs[ii, ] = exprs
      }
    }
  } else {
    dRhyExprs = matrix(ncol = nSamples, nrow = 0)
  }

  results = list(dRhyExprs = dRhyExprs,
                 featureMetadata = featureMetadata)

  return(results)
}


getNonRhyExprMatrix = function(nNonRhyGenes, nSamplesPerCond, nCond, nGenes,
                               nRhyGenes, nSamples) {
  if(nNonRhyGenes > 0) {
    nonRhyExprs = getNonRhythmicExpr(nSamplesPerCond, nCond, nGenes - nRhyGenes)
  } else {
    nonRhyExprs = matrix(ncol = nSamples, nrow = 0)
  }

  return(nonRhyExprs)
}


#' Generate differentially rhythmic groups from a set of rhythmic conditions.
#'
#' Groups will be generated via a cross product of the supplied condition
#' vectors.
#'
#' @param meanAmps is a list of floats representing the mean amplitudes of
#'   rhythmic groups.
#' @param dAmps is a list of floats representing the difference in amplitude
#'   for a differentially rhythmic group.
#' @param meanPhases is a list of floats representing the mean phases of
#'   rhythmic groups.
#' @param dPhases is a list of floats representing the difference in phase for a
#'   differentially rhythmic group.
getExprGroups = function(meanAmps = 0, dAmps = 0, meanPhases = 0,
                             dPhases = 0) {
  groups = data.table::CJ(meanAmps, dAmps, meanPhases, dPhases)
  colnames(groups) = c('meanAmp', 'dAmp', 'meanPhase', 'dPhase')
  return(groups[(dAmp > 0 | dPhase > 0) &
                ((dAmp < meanAmp * 2) |
                 (dAmp == meanAmp * 2 & dPhase == 0))])
}


#' Generate simulated gene expresion time courses.
#'
#' @param nGenes is the integer number of total genes to simulate.
#' @param nCond is the integer number of conditions to simulate.
#' @param nReps is the integer number of replicates per time point.
#' @param interval is the integer number of hours between simulated time points.
#' @param period is the integer number of hours in one rhythmic cycle.
#' @param errSd is the standard deviation of the Gaussian sample error.
#' @param rhyFrac is the proportion of genes that are generated as rhythmic.
#'   Must be a float between 0 and 1.
#' @param nSims is the integer number of simulations to generate.
#' @param drFrac is the proportion of rhythmic genes that are generated as
#'   differentially rhythmic.
#' @param exprGroups is a dataframe of metadata describing the properties of
#'   rhythmic or differentially rhythmic genes.
#' @export
getSimulatedExpr = function(nGenes = 10000L, nCond = 2, nReps = 2, interval = 4,
                            period = 24, errSd = 1, rhyFrac = 0.25, nSims = 1,
                            drFrac = 0.25,
                            exprGroups = data.table::data.table()) {

  if(period %% interval != 0) stop('period must be divisible by interval')
  if(rhyFrac < 0 | rhyFrac > 1) stop('rhyFrac must be between 0 and 1')
  if(drFrac < 0 | drFrac > 1) stop('drFrac must be between 0 and 1')

  discreteAmps = ('meanAmp' %in% colnames(exprGroups))
  if(!discreteAmps & ncol(exprGroups) > 0) {
    stop('must include meanAmp column in exprGroups')}
  if(!'dAmp' %in% colnames(exprGroups)) {
    if(ncol(exprGroups) == 0){
      exprGroups = data.table::data.table(dAmp = integer())
    } else {
      exprGroups[, dAmp := integer()]}
  }
  if(!'dPhase' %in% colnames(exprGroups)) {
    exprGroups[, dPhase := integer()]}

  exprGroups = data.table::data.table(exprGroups)

  nRhyGenes = as.integer(nGenes * rhyFrac)
  nDrGenes = as.integer(nRhyGenes * drFrac)
  nRhyOnlyGenes = nRhyGenes - nDrGenes
  nNonRhyGenes = nGenes - nRhyGenes
  nSamples = as.integer(nCond * nReps * period %/% interval)
  nSamplesPerCond = nSamples %/% nCond

  if(nrow(exprGroups) > 0) {
    onlyRhythmicGroups = data.table::data.table(
                                    meanAmp = unique(exprGroups$meanAmp),
                                    dAmp = 0,
                                    dPhase = 0)
    exprGroups = rbind(exprGroups, onlyRhythmicGroups)
  }

  if(nrow(exprGroups) > 0) {
    exprGroups[, geneIndex := 1:nrow(exprGroups)]}

  sampleNames = paste('sample', 1:nSamples, sep = '_')
  geneNames = paste('gene', 1:(nGenes * nSims), sep = '_')
  conditions = as.character(rep(1:nCond, nSamplesPerCond))

  preparedVariables = getPreparedVariables(conditions, interval, period, nCond,
                                           nReps, nRhyGenes, nGenes, nSims,
                                           nDrGenes, sampleNames, geneNames)
  sampleMetadata = preparedVariables$sampleMetadata
  featureMetadata = preparedVariables$featureMetadata
  timeSteps = preparedVariables$timeSteps
  rm(preparedVariables)

  expressionMatrix = foreach(simGroup = 1L:nSims, .combine = rbind) %do% {
    rhyResults = getRhyExprMatrix(nGenes, simGroup, nRhyOnlyGenes, nDrGenes,
                                  nSamples, discreteAmps, exprGroups,
                                  timeSteps, nCond, geneNames, conditions,
                                  nSims, featureMetadata)
    rhyExprs = rhyResults$rhyExprs
    featureMetadata = rhyResults$featureMetadata
    rm(rhyResults)

    drResults = getDrExprMatrix(simGroup, nDrGenes, exprGroups, nSamples,
                           discreteAmps, period, timeSteps, geneNames,
                           nRhyOnlyGenes, nSims, featureMetadata,
                           nSamplesPerCond, conditions, nGenes)
    dRhyExprs = drResults$dRhyExprs
    featureMetadata = drResults$featureMetadata
    rm(drResults)

    nonRhyExprs = getNonRhyExprMatrix(nNonRhyGenes, nSamplesPerCond, nCond,
                                      nGenes, nRhyGenes, nSamples)

    return(rbind(rhyExprs, dRhyExprs, nonRhyExprs) +
           matrix(rnorm(nGenes * nSamples, sd = errSd), nGenes, nSamples))
  }

  colnames(expressionMatrix) = sampleNames
  rownames(expressionMatrix) = geneNames

  values = list(
    exprs = expressionMatrix,
    sm = sampleMetadata,
    fm = featureMetadata,
    exprGroups = exprGroups
  )

  return(values)
}
