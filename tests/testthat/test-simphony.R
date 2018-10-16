context('simphony')


test_that('Single condition simulation works', {
  exprGroupsList = list(data.table::data.table(base = 1, amp = 3, phase = 2))
  expect_silent({simGse = simulateGeneData(exprGroupsList)})

  rm(exprGroupsList, simGse)
})

test_that('Multiple condition simulation works', {
  exprGroupsList = list(data.table::data.table(base = c(1, 2)),
                        data.table::data.table(amp = c(1,2), phase = c(3,4)))
  expect_silent(simulateGeneData(exprGroupsList))
  expect_silent(simulateGeneData(exprGroupsList, method = 'negbinom'))

  rm(exprGroupsList)
})

test_that('Number of genes and samples simulated are predictable', {
  nGenes = 100
  sampleInterval = 6
  nReps = 3

  exprGroupsList = list(data.table::data.table(base = c(1, 2)),
                        data.table::data.table(amp = c(1,2), phase = c(3,4)))
  simGse = simulateGeneData(exprGroupsList, nGenes = nGenes,
                            interval = sampleInterval, nReps = nReps)

  expect_equal(ncol(simGse$emat), nrow(simGse$sm))
  expect_equal(ncol(simGse$emat),
               length(exprGroupsList) * nReps * as.integer(24 / sampleInterval))
  expect_equal(nrow(simGse$emat), nGenes)
  expect_equal(nrow(simGse$emat), nrow(simGse$gm) / length(exprGroupsList))

  rm(exprGroupsList, simGse)
})
