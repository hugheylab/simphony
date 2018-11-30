library('data.table')

# Simulate data for genes having one of three sets of rhythmic parameters.
exprGroups = data.table(amp = c(0, 1, 1), phase = c(0, 0, 6),
                        rhyFunc = c(cos, cos, sin))
simData = simphony(exprGroups)

# Simulate data for genes whose rhythmicity varies between two conditions.
exprGroupsList = list(data.table(amp = c(1, 2), phase = c(0, -3)),
                      data.table(amp = c(3, 2), phase = c(0, 3)))
simData = simphony(exprGroupsList)

# Simulate data for 100 genes, half non-rhythmic and half rhythmic, with
# amplitudes for rhythmic genes sampled from a log-normal distribution.
nGenes = 100
rhyFrac = 0.5
nRhyGenes = round(rhyFrac * nGenes)
rhyAmps = 2^rnorm(nRhyGenes, sd = 0.415, mean = 1.32)
fracGenes = c(1 - rhyFrac, rep(rhyFrac / nRhyGenes, nRhyGenes))
exprGroups = data.table(amp = c(0, rhyAmps), fracGenes = fracGenes)
simData = simphony(exprGroups, nGenes = nGenes)

# Simulate data for 100 rhythmic genes, with expression values sampled from
# the negative binomial family and with baseline log2 expected counts sampled
# from a normal distribution whose parameters were estimated using DESeq2 and
# circadian RNA-seq data from mouse lung (NCBI GEO GSE54651).
baseLogCounts = rnorm(100, mean = 8.96, sd = 2.27)
exprGroups = data.table(base = baseLogCounts, amp = 1)
simData = simphony(exprGroups, nGenes = nrow(exprGroups), family = 'negbinom')

# Simulate data for 100 rhythmic genes, with expression values sampled from
# the negative binomial family and with variance of residual log dispersion
# sampled from a normal distribution whose parameters were estimated
# using DESeq2 and circadian RNA-seq data from mouse lung (NCBI GEO GSE54651).
dispFactors = exp(rnorm(100, sd = sqrt(1.02)))
dispFuncs = sapply(dispFactors, function(z) getDispFunc(x2 = z))
exprGroups = data.table(dispFunc = dispFuncs, amp = 1)
simData = simphony(exprGroups, nGenes = nrow(exprGroups), family = 'negbinom')
