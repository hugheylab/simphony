context('simulation')

test_that('Outputs are of correct type', {
  exprGroups = data.table::data.table(meanAmp = c(1, 2), dAmp = c(1, 1))
  simGse = getSimulatedExpr(exprGroups)

  expect_is(simGse$sm, 'data.table')
  expect_is(simGse$emat, 'matrix')
  expect_is(simGse$gm, 'data.table')
  expect_is(simGse$exprGroups, 'data.table')
  
  rm(exprGroups, simGse)
})

test_that('All expected columns are accounted for in exprGroups output', {
  expected_columns = c('geneCount', 'meanExpr', 'dExpr', 'meanAmp', 'dAmp',
                       'meanPhase', 'dPhase', 'group')
  exprGroups = data.table::data.table(meanAmp = c(1, 2), dAmp = c(1, 1))
  simGse = getSimulatedExpr(exprGroups)

  expect_true(all(expected_columns %in% colnames(simGse$exprGroups)))
  
  rm(expected_columns, exprGroups, simGse)
})

test_that('exprGroups output has correct group column', {
  exprGroups = data.table::data.table(meanAmp = c(1, 2, 3), dAmp = c(1, 1, 2))
  simGse = getSimulatedExpr(exprGroups)

  expect_equal(simGse$exprGroups$group, c(1, 2, 3))

  rm(exprGroups, simGse)
})

test_that('Correct number of genes are created', {
  exprGroups = data.table::data.table(meanAmp = c(1, 2, 3), dAmp = c(1, 1, 2))
  simGse = getSimulatedExpr(exprGroups, nGenes = 1000)

  expect_equal(simGse$exprGroups[, sum(geneCount)], 1000)
  expect_equal(nrow(simGse$emat), 1000)

  rm(exprGroups, simGse)
})

test_that('exprGroups can be passed as a data.table', {
  exprGroups = data.frame(meanAmp = c(1, 2, 3), dAmps = c(1, 1, 2))

  expect_silent(getSimulatedExpr(exprGroups))

  rm(exprGroups)
})

test_that('Correct number of samples created', {
  exprGroups = data.table::data.table(meanAmp = c(1, 3, 5))

  simGse = getSimulatedExpr(exprGroups, interval = 3)
  expect_equal(ncol(simGse$emat), (24 / 3) * 2 * 2)

  simGse = getSimulatedExpr(exprGroups, interval = 3, nReps = 4)
  expect_equal(ncol(simGse$emat), (24 / 3) * 2 * 4)

  rm(exprGroups, simGse)
})

test_that('Changing nSims yields expected results', {
  exprGroups = data.table::data.table(meanPhase = c(3, 6))
  simGse = getSimulatedExpr(exprGroups, nSims = 10, nGenes = 100)

  expect_equal(nrow(simGse$emat), 1000)

  rm(exprGroups, simGse)
})

test_that('Rhythmic function can be changed', {
  exprGroups = data.table::data.table(meanPhase = c(3, 6))
  expect_silent(getSimulatedExpr(exprGroups, rhyFunc = cos))

  rm(exprGroups)
})

test_that('Random timepoints can be used', {
  exprGroups = data.table::data.table(meanExpr = c(1,2,3))
  expect_silent(getSimulatedExpr(exprGroups, randomTimepoints = TRUE,
                                         nSamples = 10))

  rm(exprGroups)
})

test_that('Appropriate errors are thrown', {
  expect_error(getSimulatedExpr(data.table::data.table()),
               'No rows in exprGroups. Cannot simulate genes.')
  expect_error(getSimulatedExpr(data.table::data.table(meanAmp = c(1,2)),
                                        randomTimepoints = TRUE),
               'Number of random timepoint samples not specified.')
})
