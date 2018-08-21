context('simulation')
library(data.table)
library(foreach)

test_that('Rhythmic groups created correctly', {
  meanAmps = c(0, 1, 2, 3)
  dAmps = c(0, 0.5, 1, 1.5, 2)
  dPhases = c(0, 1, 2, 3, 4, 5)
  rhythmicGroups = getRhythmicGroups(meanAmps, dAmps, dPhases)

  expect_equal(nrow(rhythmicGroups), 82)

  expect_equal(nrow(rhythmicGroups[meanAmp == 0]), 0)
  expect_equal(nrow(rhythmicGroups[meanAmp == 1]), 23 + 1)
  expect_equal(nrow(rhythmicGroups[meanAmp == 2]), 29)

  expect_equal(nrow(rhythmicGroups[dAmp == 0 & dPhase == 0]), 0)
})

test_that('Simulated GSE classes are correct', {
  simGse = getSimulatedExpr(nGenes = 100)

  expect_is(simGse$sm, 'data.table')
  expect_is(simGse$fm, 'data.table')
  expect_is(simGse$exprs, 'matrix')

  expect_is(simGse$sm$time, 'numeric')
  expect_is(simGse$sm$cond, 'character')
  expect_is(simGse$sm$sample, 'character')
  
  rm(simGse)
})

test_that('getRhythmicExpr yields correct values', {
  timeSteps = seq(0, 2 * pi, 0.1) # length(timeSteps) is 63.
  phase = pi
  amplitude = 3
  verticalShift = 3
  expr = getRhythmicExpr(timeSteps, phase, amplitude, verticalShift)

  expect_equal(length(expr), length(timeSteps))
  expect_equal(sum(expr > amplitude + verticalShift), 0)
  expect_equal(sum(expr < -(amplitude + verticalShift)), 0)

  rm(timeSteps, phase, amplitude, verticalShift, expr)
})

test_that('getNonRhythmicExpr yields correct values', {
  nSamples = 150
  nCond = 2
  nRows = 150
  expression = 5
  expr = getNonRhythmicExpr(nSamples, nCond, nRows, expression)

  expect_equal(ncol(expr), nSamples * nCond)
  expect_equal(nrow(expr), nRows)
  expect_equal(length(unique(c(unique(expr)))), 1)
  expect_equal(unique(c(unique(expr))), expression)

  rm(nSamples, nRows, expression, expr)
})

test_that('Non-discrete getExprAmplitudes yields correct results', {
  simAmps = getExprAmplitudes(1000)

  expect_equal(length(simAmps), 1000)
  expect_equal(sum(simAmps < 0.5), 0)

  rm(simAmps)
})

test_that('Discrete getExprAmplitudes yields correct results', {
  numDistinctAmps = 17
  numExprAmps = 167
  amps = seq(numDistinctAmps)
  simAmps = getExprAmplitudes(numExprAmps, TRUE, amps)

  expect_equal(length(simAmps),
               numDistinctAmps * (numExprAmps %/% numDistinctAmps + 1))
  expect_equal(unique(simAmps), amps)

  rm(numDistinctAmps, numExprAmps, amps, simAmps)
})

test_that('SM times are identical for all conditions', {
  simGse = getSimulatedExpr(nGenes = 100)

  times = as.list(split(simGse$sm$time, simGse$sm$cond))
  expect_equal(length(unique(times)), 1)

  rm(simGse)
})

test_that('Specifying only 1 sim yields expected outputs', {
  simGse = getSimulatedExpr(nGenes = 100, nSims = 1)

  expect_equal(length(unique(simGse$fm$simIndex)), 1)
  expect_equal(nrow(simGse$fm), nrow(simGse$exprs))
  expect_equal(nrow(simGse$fm), length(unique(rownames(simGse$exprs))))

  simGseDefault = getSimulatedExpr(nGenes = 100)

  expect_equal(simGse$fm, simGseDefault$fm)
  expect_equal(rownames(simGse$exprs), rownames(simGseDefault$exprs))

  rm(simGse)
  rm(simGseDefault)
})

test_that('Correct number of simulations are created', {
  simGse = getSimulatedExpr(nGenes = 100, nSims = 5)

  expect_equal(length(unique(simGse$fm$simIndex)), 5)
  expect_equal(nrow(simGse$fm), nrow(simGse$exprs))
  expect_equal(simGse$fm$gene, rownames(simGse$exprs))
  expect_equal(nrow(simGse$fm), 500)

  rm(simGse)
})

test_that('Order of genes is consistent for emat/fm', {
  simGse = getSimulatedExpr(nGenes = 100)

  expect_equal(simGse$fm$gene, rownames(simGse$exprs))

  rm(simGse)
})

test_that('Number of samples is nCond * nReps * period / interval', {
  nCond = 2
  nReps = 3
  period = 24
  interval = 6
  expectedSampleCount = nCond * nReps * period / interval

  simGseSamples = getSimulatedExpr(nGenes = 100, nReps = nReps, nCond = nCond,
                                       period = period, interval = interval)

  expect_equal(ncol(simGseSamples$exprs), expectedSampleCount)
  expect_equal(nrow(simGseSamples$sm), expectedSampleCount)

  rm(simGseSamples)
})

test_that('Correct number of rhy/dr genes is generated', {
  simGse = getSimulatedExpr(nGenes = 100)

  expect_equal(sum(simGse$fm$rhy == 1), 25)
  expect_equal(sum(simGse$fm$dRhy == 1), 6)

  rm(simGse)
})

test_that('No rhy or dr genes when rhyFrac = 0', {
  simGseNoRhy = getSimulatedExpr(nGenes = 100, rhyFrac = 0)

  expect_equal(sum(simGseNoRhy$fm$rhy == 1), 0)
  expect_equal(sum(simGseNoRhy$fm$drhy == 1), 0)

  rm(simGseNoRhy)
})

test_that('No dr genes when drFrac = 0', {
  simGseNoDr = getSimulatedExpr(nGenes = 100, rhyFrac = 0.2, drFrac = 0)

  expect_equal(sum(simGseNoDr$fm$rhy == 1), 20)
  expect_equal(sum(simGseNoDr$fm$drhy == 1), 0)

  rm(simGseNoDr)
})

test_that('Simuation errors are thrown', {
  expect_error(getSimulatedExpr(nGenes = 100, period = 24, interval = 5),
               'period must be divisible by interval')
  expect_error(getSimulatedExpr(nGenes = 100, rhyFrac = -0.5),
               'rhyFrac must be between 0 and 1')
  expect_error(getSimulatedExpr(nGenes = 100, rhyFrac = 1.5),
               'rhyFrac must be between 0 and 1')
  expect_error(getSimulatedExpr(nGenes = 100, drFrac = -0.5),
               'drFrac must be between 0 and 1')
  expect_error(getSimulatedExpr(nGenes = 100, drFrac = 1.5),
               'drFrac must be between 0 and 1')
  expect_error(getSimulatedExpr(nGenes = 100,
                                rhythmicGroups = data.frame(a = c(1))),
               'must include meanAmp column in rhythmicGroups')
})

test_that('Prepared variables are created correctly', {
  nCond = 2
  interval = 3
  period = 24
  nReps = 2
  nGenes = 2560
  nRhyGenes = 160
  nDrGenes = 10
  nSims = 4
  nSamples = as.integer(nCond * nReps * period %/% interval)
  conditions = as.character(rep(1:nCond, nSamples / nCond))
  sampleNames = paste('sample', 1:nSamples, sep = '_')
  geneNames = paste('gene', 1:(nGenes * nSims), sep = '_')

  preparedVariables = getPreparedVariables(conditions, interval, period, nCond,
                                           nReps, nRhyGenes, nGenes, nSims,
                                           nDrGenes, sampleNames, geneNames)

  expect_equal(names(preparedVariables), c('sampleMetadata', 'featureMetadata',
                                           'timeSteps'))

  # Sample Metadata Tests
  expect_equal(nrow(preparedVariables$sampleMetadata), nSamples)
  expect_equal(preparedVariables$sampleMetadata$sample, sampleNames)
  expect_equal(preparedVariables$sampleMetadata,
               preparedVariables$sampleMetadata[order(cond, time)])
  expect_true(!(period %in% preparedVariables$sampleMetadata$time))

  # Feature Metadata Tests
  expect_equal(nrow(preparedVariables$featureMetadata), nGenes * nSims)
  expect_equal(nrow(preparedVariables$featureMetadata[rhy == 1]),
               nRhyGenes * nSims)
  expect_equal(nrow(preparedVariables$featureMetadata[dRhy == 1]),
               nDrGenes * nSims)
  expect_equal(length(unique(preparedVariables$featureMetadata$simIndex)),
               nSims)
  expect_equal(preparedVariables$featureMetadata$gene, geneNames)
})
