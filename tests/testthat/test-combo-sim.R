context('simulation')

test_that('Combo sim works', {
  exprGroups = list(data.table::data.table(base = c(1, 2)),
                    data.table::data.table(amp = c(1,2), phase = c(3,4)))
  expect_silent(getMultipleCondSimulation(exprGroups))
})