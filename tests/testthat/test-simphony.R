context('simphony')


test_that('Single condition simulation works', {
  exprGroupsList = data.table::data.table(base = 1, amp = 3, phase = 2)
  expect_silent({simGse = simulateExprData(exprGroupsList)})

  rm(exprGroupsList, simGse)
})

test_that('Multiple condition simulation works', {
  exprGroupsList = list(data.table::data.table(base = c(1, 2)),
                        data.table::data.table(amp = c(1,2), phase = c(3,4)))
  expect_silent(simulateExprData(exprGroupsList))
  expect_silent(simulateExprData(exprGroupsList, family = 'negbinom'))

  rm(exprGroupsList)
})

test_that('Number of genes and samples simulated are predictable', {
  nGenes = 100
  sampleInterval = 6
  nReps = 3

  exprGroupsList = list(data.table::data.table(base = c(1, 2)),
                        data.table::data.table(amp = c(1,2), phase = c(3,4)))
  simGse = simulateExprData(exprGroupsList, nGenes = nGenes,
                            interval = sampleInterval, nReps = nReps)

  expect_equal(ncol(simGse$exprData), nrow(simGse$sampleMetadata))
  expect_equal(ncol(simGse$exprData),
               length(exprGroupsList) * nReps * as.integer(24 / sampleInterval))
  expect_equal(nrow(simGse$exprData), nGenes)
  expect_equal(nrow(simGse$exprData), nrow(simGse$geneMetadata) / length(exprGroupsList))

  rm(exprGroupsList, simGse)
})

test_that('Appropriate errors are thrown', {
  badExprGroupsList = list(data.table::data.table(base = c(1, 2)),
                           data.table::data.table(base = 1))
  expect_error(simulateExprData(badExprGroupsList),
               'Number of rows in each exprGroups must be the same for all conditions')

  goodExprGroupsList = list(data.table::data.table(base = c(2)),
                           data.table::data.table(base = 1))
  expect_error(simulateExprData(goodExprGroupsList, family = 'socratic'),
               'family must be \'gaussian\' or \'negbinom\'.')
})
