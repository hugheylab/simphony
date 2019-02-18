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
