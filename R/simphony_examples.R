library('data.table')

# Simulate data for genes having one of three sets of rhythmic parameters.
exprGroups = data.table(amp = c(0, 1, 1), phase = c(0, 0, 6),
                        rhyFunc = c(cos, cos, sin))
simData = simphony(exprGroups)

# Simulate data for an experiment with specified timepoints and replicates.
exprGroups = data.table(amp = c(0, 1))
simData = simphony(exprGroups, timepointsType = 'specified',
                   timepoints = rep(seq(0, 6, 2), each = 2))

# Simulate data for genes whose rhythmicity varies between two conditions.
exprGroupsList = list(data.table(amp = c(1, 2), phase = c(0, -3)),
                      data.table(amp = c(3, 2), phase = c(0, 3)))
simData = simphony(exprGroupsList)

# Simulate data for 100 genes, half non-rhythmic and half rhythmic, with
# amplitudes for rhythmic genes sampled from a log-normal distribution whose
# parameters were estimated using limma and circadian microarray data from
# mouse lung (GSE54650).
nGenes = 100
rhyFrac = 0.5
nRhyGenes = round(rhyFrac * nGenes)
rhyAmps = 2^rnorm(nRhyGenes, sd = 0.415, mean = 1.32)
fracGenes = c(1 - rhyFrac, rep(rhyFrac / nRhyGenes, nRhyGenes))
exprGroups = data.table(amp = c(0, rhyAmps), fracGenes = fracGenes)
simData = simphony(exprGroups, nGenes = nGenes)

# Simulate data for 100 rhythmic genes, with baseline log2 expected counts and
# residual log dispersion sampled from distributions whose parameters were
# estimated using DESeq2 and circadian RNA-seq data from mouse lung (GSE54651).
nGenes = 100
baseLogCounts = rnorm(nGenes, mean = 8.63, sd = 2.63)
dispFactors = exp(rnorm(nGenes, sd = 1.02))
dispFuncs = sapply(dispFactors, function(z) getDispFunc(x2 = z))
exprGroups = data.table(base = baseLogCounts, dispFunc = dispFuncs, amp = 1)
simData = simphony(exprGroups, nGenes = nGenes, family = 'negbinom')
