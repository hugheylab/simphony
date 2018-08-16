#' @importFrom data.table ":="
#' @importFrom foreach "%do%"


getRhythmicExpr = function(timeSteps, phase = 0, amplitude = 1,
                           verticalShift = 0) {

  exprs = amplitude * sin(timeSteps + phase) + verticalShift
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

  sampleMetadata = data.table(cond = conditions,
                              time = rep(interval * 0:(period / interval - 1),
                                         each = nCond * nReps))
  sampleMetadata = sampleMetadata[order(cond, time)]
  sampleMetadata[, sample := sampleNames]

  featureMetadata = data.table(geneId = geneNames,
                               rhy = rep(c(rep(1, nRhyGenes),
                                           rep(0, nGenes - nRhyGenes)), nSims),
                               dRhy = rep(c(rep(0, nRhyGenes - nDrGenes),
                                            rep(1, nDrGenes),
                                            rep(0, nGenes - nRhyGenes)), nSims),
                               rhyIndex = rep(rep(as.integer(NA), nGenes),
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
                            discreteAmps, rhythmicGroups, timeSteps, nCond,
                            geneNames, conditions, nSims, featureMetadata) {

  nRhyGenesPerSim = nRhyOnlyGenes + nDrGenes
  if(nRhyOnlyGenes > 0) {
    nDrGroups = nrow(rhythmicGroups[dAmp != 0 | dPhase != 0, ])
    nRhyGroups = nrow(rhythmicGroups[dAmp == 0 & dPhase == 0, ])

    rhyExprs = matrix(ncol = nSamples, nrow = nRhyOnlyGenes)
    amplitudes = getExprAmplitudes(nRhyOnlyGenes, discrete = discreteAmps,
                                   amps = t(rhythmicGroups[dAmp == 0 & dPhase == 0, meanAmp]))
    for(ii in 1L:nRhyOnlyGenes) {
      condExprs = getRhythmicExpr(timeSteps, amplitude = amplitudes[ii])
      rhyExprs[ii, ] = rep(condExprs, nCond)
      if(discreteAmps) {
        groupIndex = nDrGroups + (ii - 1) %% nRhyGroups + 1
        set(featureMetadata, i = ii + (simGroup - 1L) * nGenes, j = 'rhyIndex',
            value = rhythmicGroups[groupIndex, rhyIndex])
      }
    }
  } else {
    rhyExprs = matrix(ncol = nSamples, nrow = 0)
  }

  results = list(rhyExprs = rhyExprs,
                 featureMetadata = featureMetadata)
  return(results)
}


getDrExprMatrix = function(simGroup, nDrGenes, rhythmicGroups, nSamples,
                           discreteAmps, period, timeSteps, geneNames,
                           nRhyOnlyGenes, nSims, featureMetadata,
                           nSamplesPerCond, conditions, nGenes) {

  nRhyGenesPerSim = nRhyOnlyGenes + nDrGenes
  if(nDrGenes > 0) {
    nGenesPerGroup = nDrGenes / as.integer(nrow(rhythmicGroups[dAmp != 0 | dPhase != 0, ]))
    dRhyExprs = matrix(nrow = nDrGenes, ncol = nSamples)

    if(discreteAmps) {
      for(row in 1L:nrow(rhythmicGroups[dAmp != 0 | dPhase != 0, ])) {
        condAmps = c(rhythmicGroups[row, meanAmp] -
                     rhythmicGroups[row, dAmp] / 2,
                     rhythmicGroups[row, meanAmp] +
                     rhythmicGroups[row, dAmp] / 2)
        condPhases = c(0, (2 * pi / period) * rhythmicGroups[row, dPhase])

        baseExpr = c(getRhythmicExpr(timeSteps, phase = condPhases[1],
                                      amplitude = condAmps[1]),
                      getRhythmicExpr(timeSteps, phase = condPhases[2],
                                      amplitude = condAmps[2]))

        for(ii in 1L:nGenesPerGroup) {
          rowIndex = (row - 1L) * nGenesPerGroup + ii
          dRhyExprs[rowIndex, ] = baseExpr
          gene = geneNames[rowIndex + nRhyOnlyGenes + (simGroup - 1) * nSims]

          set(featureMetadata,
              i = rowIndex + nRhyOnlyGenes + (simGroup - 1L) * nGenes,
              j = 'rhyIndex', value = rhythmicGroups[row, rhyIndex])
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


#' Generate simulated gene expresion time courses.
#'
#' \code{nGenes} is the integer number of total genes to simulate.
#' \code{nCond} is the integer number of conditions to simulate.
#' \code{nReps} is the integer number of replicates per time point.
#' \code{interval} is the integer number of hours between simulated time points.
#' @export
getSimulatedExpr = function(nGenes = 10000L, nCond = 2, nReps = 2, interval = 4,
                            period = 24, errSd = 1, rhyFrac = 0.25, nSims = 1,
                            drFrac = 0.25, rhythmicGroups = data.frame()) {

  if(period %% interval != 0) stop('period must be divisible by interval')
  if(rhyFrac < 0 | rhyFrac > 1) stop('rhyFrac must be between 0 and 1')
  if(drFrac < 0 | drFrac > 1) stop('drFrac must be between 0 and 1')

  discreteAmps = ('meanAmp' %in% colnames(rhythmicGroups))
  if(!discreteAmps & ncol(rhythmicGroups) > 0) {
    stop('must include meanAmp column in rhythmicGroups')}
  if(!'dAmp' %in% colnames(rhythmicGroups)) {
    rhythmicGroups['dAmp'] = rep(0, nrow(rhythmicGroups))}
  if(!'dPhase' %in% colnames(rhythmicGroups)) {
    rhythmicGroups['dPhase'] = rep(0, nrow(rhythmicGroups))}

  rhythmicGroups = data.table(rhythmicGroups)

  nRhyGenes = as.integer(nGenes * rhyFrac)
  nDrGenes = as.integer(nRhyGenes * drFrac)
  nRhyOnlyGenes = nRhyGenes - nDrGenes
  nNonRhyGenes = nGenes - nRhyGenes
  nSamples = as.integer(nCond * nReps * period %/% interval)
  nSamplesPerCond = nSamples %/% nCond

  if(discreteAmps & nDrGenes %% nrow(rhythmicGroups) != 0) {
    stop('number of rows in rhythmicGroups must be divisible by number of DR genes')}

  if(nrow(rhythmicGroups) > 0) {
    onlyRhythmicGroups = data.table(meanAmp = unique(rhythmicGroups$meanAmp),
                                    dAmp = 0,
                                    dPhase = 0)
    rhythmicGroups = rbind(rhythmicGroups, onlyRhythmicGroups)
  }

  if(discreteAmps) {
    rhythmicGroups[, rhyIndex := 1:nrow(rhythmicGroups)]}

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
                                  nSamples, discreteAmps, rhythmicGroups,
                                  timeSteps, nCond, geneNames, conditions,
                                  nSims, featureMetadata)
    rhyExprs = rhyResults$rhyExprs
    featureMetadata = rhyResults$featureMetadata
    rm(rhyResults)

    drResults = getDrExprMatrix(simGroup, nDrGenes, rhythmicGroups, nSamples,
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

  rhythmicGroups = data.frame(rhythmicGroups)

  colnames(expressionMatrix) = sampleNames
  rownames(expressionMatrix) = geneNames

  values = list(
    exprs = expressionMatrix,
    sm = sampleMetadata,
    fm = featureMetadata,
    rhythmicGroups = rhythmicGroups
  )

  return(values)
}
