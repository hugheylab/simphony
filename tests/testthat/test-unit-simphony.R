context('Simphony Unit Tests')


test_that('Single condition simulation works', {
  featureGroupsList = data.table::data.table(base = 1, amp = 3, phase = 2)
  expect_silent({simData = simphony(featureGroupsList)})

  rm(featureGroupsList, simData)
})

test_that('Multiple condition simulation works', {
  featureGroupsList = list(data.table::data.table(base = c(1, 2)),
                           data.table::data.table(amp = c(1, 2), phase = c(3, 4)))
  expect_silent(simphony(featureGroupsList))
  expect_silent(simphony(featureGroupsList, family = 'negbinom'))

  rm(featureGroupsList)
})

test_that('Number of features and samples simulated are predictable', {
  nFeatures = 100
  sampleInterval = 6
  nReps = 3
  timeRange = c(0, 24)

  featureGroupsList = list(data.table::data.table(base = c(1, 2)),
                           data.table::data.table(amp = c(1, 2), phase = c(3, 4)))
  simData = simphony(featureGroupsList, nFeatures = nFeatures, timeRange = timeRange,
                     interval = sampleInterval, nReps = nReps)

  expect_equal(ncol(simData$abundData), nrow(simData$sampleMetadata))
  expect_equal(ncol(simData$abundData),
               length(featureGroupsList) * nReps * as.integer((timeRange[2] - timeRange[1]) / sampleInterval))
  expect_equal(nrow(simData$abundData), nFeatures)
  expect_equal(nrow(simData$abundData),
               nrow(simData$featureMetadata) / length(featureGroupsList))

  rm(featureGroupsList, simData)
})

test_that('Appropriate errors are thrown', {
  badFeatureGroupsList = list(data.table::data.table(base = c(1, 2)),
                              data.table::data.table(base = 1))
  expect_error(simphony(badFeatureGroupsList),
               'Each featureGroups data.frame must have the same number of rows.')
  rm(badFeatureGroupsList)
})

test_that('Condition names are as expected', {
  featureGroupsList = foreach(cond = 1:20) %do% {
    data.table(amp = cond)
  }
  simData = simphony(featureGroupsList)
  expectedCondNames = sprintf(sprintf('cond_%%0%dd', floor(log10(20)) + 1),
                              1:20)

  expect_equal(unique(simData$sampleMetadata[, cond]), expectedCondNames)
  expect_equal(unique(simData$featureMetadata[, cond]), expectedCondNames)

  rm(featureGroupsList, simData, expectedCondNames)
})

test_that('Sample names are as expected', {
  featureGroups = data.table(amp = 1)
  simData = simphony(featureGroups, timepointsType = 'specified',
                     timepoints = 1:20)
  expectedSampleNames = sprintf(sprintf('sample_%%0%dd', floor(log10(20)) + 1),
                            1:20)

  expect_equal(unique(simData$sampleMetadata[, sample]), expectedSampleNames)

  rm(featureGroups, simData, expectedSampleNames)
})

test_that('Feature names are as expected', {
  featureGroups = data.table(amp = 1)
  simData = simphony(featureGroups, nFeatures = 20)
  expectedFeatureNames = sprintf(sprintf('feature_%%0%dd', floor(log10(20)) + 1),
                                 1:20)

  expect_equal(unique(simData$featureMetadata[, feature]), expectedFeatureNames)

  rm(featureGroups, simData, expectedFeatureNames)
})
