context('simphony')


test_that('Single condition simulation works', {
  exprGroupsList = data.table::data.table(base = 1, amp = 3, phase = 2)
  expect_silent({simGse = simphony(exprGroupsList)})

  rm(exprGroupsList, simGse)
})

test_that('Multiple condition simulation works', {
  exprGroupsList = list(data.table::data.table(base = c(1, 2)),
                        data.table::data.table(amp = c(1,2), phase = c(3,4)))
  expect_silent(simphony(exprGroupsList))
  expect_silent(simphony(exprGroupsList, family = 'negbinom'))

  rm(exprGroupsList)
})

test_that('Number of features and samples simulated are predictable', {
  nFeatures = 100
  sampleInterval = 6
  nReps = 3

  exprGroupsList = list(data.table::data.table(base = c(1, 2)),
                        data.table::data.table(amp = c(1,2), phase = c(3,4)))
  simGse = simphony(exprGroupsList, nFeatures = nFeatures,
                    interval = sampleInterval, nReps = nReps)

  expect_equal(ncol(simGse$exprData), nrow(simGse$sampleMetadata))
  expect_equal(ncol(simGse$exprData),
               length(exprGroupsList) * nReps * as.integer(24 / sampleInterval))
  expect_equal(nrow(simGse$exprData), nFeatures)
  expect_equal(nrow(simGse$exprData),
               nrow(simGse$featureMetadata) / length(exprGroupsList))

  rm(exprGroupsList, simGse)
})

test_that('Statistics from NBD are as expected', {
  base = c(3, 5, 7)
  expectedVariance = 2^base + defaultDispFunc(2^base) * ((2^base)^2)

  exprGroupsList = data.table::data.table(base = base, amp = 0)
  simGse = simphony(exprGroupsList, nFeatures = 3, nReps = 500, family = 'negbinom')

  exprData = data.table::data.table(expr = c(t(simGse$exprData)),
                                    feature = rep(rownames(simGse$exprData),
                                               each = ncol(simGse$exprData)))
  expect_equal(exprData[, log2(mean(expr)), by = feature][, V1], base, tolerance = 1e-1)
  expect_equal(exprData[, var(expr), by = feature][, V1], expectedVariance, tolerance = 1e-1)

  rm(base, expectedVariance, exprGroupsList, simGse, exprData)
})

test_that('Appropriate errors are thrown', {
  badExprGroupsList = list(data.table::data.table(base = c(1, 2)),
                           data.table::data.table(base = 1))
  expect_error(simphony(badExprGroupsList),
               'Each exprGroups data.frame must have the same number of rows.')
})
