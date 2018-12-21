context('Simphony Functional Tests')


test_that('Abundances are sampled from the correct trend', {
  featureGroups = data.table::data.table(amp = seq(1, 10, 1), base = 1:10, sd = 0)
  timepoints = seq(0, 22, 0.1)
  simData = simphony(featureGroups, nFeatures = 10, timepoints = timepoints, timepointsType = 'specified')

  usedAbundGroups = simData$experMetadata$featureGroupsList
  usedAmps = usedAbundGroups[[1]][, amp]
  usedRhyFunc = usedAbundGroups[[1]][, rhyFunc]
  usedBase = usedAbundGroups[[1]][, base]
  usedPeriod = simData$experMetadata$period

  expectedAbund = foreach(r = 1:nrow(simData$abundData), .combine = rbind) %do% {
    usedAmps[r] * usedRhyFunc[[r]]((2 * pi) * timepoints / usedPeriod) + usedBase[r]
  }
  diff = abs(simData$abundData - expectedAbund)
  expect_true(all(diff == 0))

  rm(featureGroups, timepoints, simData, usedAbundGroups, usedAmps, usedRhyFunc,
     usedBase, usedPeriod, expectedAbund, diff)
})

test_that('Time-independent statistics from NBD are as expected', {
  base = c(3, 5, 7)
  expectedVariance = 2 ^ base + defaultDispFunc(2 ^ base) * ((2 ^ base) ^ 2)

  featureGroups = data.table::data.table(base = base, amp = 0)
  simData = simphony(featureGroups, nFeatures = 3, nReps = 500, family = 'negbinom')

  abundData = data.table::data.table(abund = c(t(simData$abundData)),
                                     feature = rep(rownames(simData$abundData),
                                                   each = ncol(simData$abundData)))
  expect_equal(abundData[, log2(mean(abund)), by = feature][, V1], base, tolerance = 1e-1)
  expect_equal(abundData[, var(abund), by = feature][, V1], expectedVariance, tolerance = 1e-1)

  rm(base, expectedVariance, featureGroups, simData, abundData)
})

test_that('Time-dependent statistics from NBD are as expected', {
  featureGroups = data.table::data.table(amp = 3, base = 4:8)
  simData = simphony(featureGroups, nFeatures = 5, nReps = 4000, family = 'negbinom')

  for(timeNow in unique(simData$sampleMetadata$time)) {
    samplesNow = simData$sample[, time == timeNow]

    for(groupNow in unique(simData$featureMetadata$group)) {
      featuresNow = simData$featureMetadata[, group == groupNow]
      params = simData$featureMetadata[featuresNow, ][1, ]

      meanExp = 2 ^ (params$amp * params$rhyFunc[[1]](2 * pi * timeNow / simData$experMetadata$period) + params$base)
      varExp = meanExp + params$dispFunc[[1]](meanExp) * (meanExp ^ 2)

      meanObs = mean(simData$abundData[featuresNow, samplesNow])
      varObs = var(simData$abundData[featuresNow, samplesNow])

      expect_equal(meanExp, meanObs, tolerance = 0.1)
      expect_equal(varExp, varObs, tolerance = 0.1)
    }
  }

  rm(featureGroups, simData, samplesNow, featuresNow, params, meanExp, varExp,
     meanObs, varObs)
})