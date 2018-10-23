defaultDispFunc = function(x) {
  return(3/x)
}


simulateExprDataOneCond = function(exprGroups, times, method) {
  foreach(group = 1:nrow(exprGroups), .combine = rbind) %do% {
    amp = exprGroups[group, amp]
    phase = exprGroups[group, phase]
    base = exprGroups[group, base]
    rhyFunc = exprGroups[group, rhyFunc][[1]]
    mu = amp * rhyFunc(times + 2 * pi * phase) + base

    if(method == 'gaussian') {
      groupEmat = stats::rnorm(length(mu) * exprGroups[group, numGenes],
                               rep(mu, exprGroups[group, numGenes]),
                               sd = exprGroups[group, sd])
    } else {
      dispFunc = exprGroups[group, dispFunc][[1]]
      groupEmat = stats::rnbinom(length(mu) * exprGroups[group, numGenes],
                                 mu = 2^rep(mu, exprGroups[group, numGenes]),
                                 size = 1/dispFunc(2^mu)) }
    matrix(groupEmat, nrow = exprGroups[group, numGenes], byrow = TRUE)
  }
}


setDefaultExprGroups = function(exprGroups, nGenes, rhyFunc, method) {
  exprGroups = data.table(exprGroups)
  exprGroups[, group := 1:.N]

  if(nrow(exprGroups) == 0) {
    stop('exprGroups must have at least one row.') }

  if(!'fracGenes' %in% colnames(exprGroups)) {
    exprGroups[, fracGenes := 1 / nrow(exprGroups)] }

  if(!'amp' %in% colnames(exprGroups)) {
    exprGroups[, amp := 0] }

  if(!'phase' %in% colnames(exprGroups)) {
    exprGroups[, phase := 0] }

  if(!'rhyFunc' %in% colnames(exprGroups)) {
    exprGroups[, rhyFunc := data.table(rhyFunc)] }

  if(method == 'negbinom') {
    if(!'dispFunc' %in% colnames(exprGroups)) {
      exprGroups[, dispFunc := data.table(defaultDispFunc)] }
    if(!'base' %in% colnames(exprGroups)) {
      exprGroups[, base := 7] } }
  else {
    if(!'sd' %in% colnames(exprGroups)) {
      exprGroups[, sd := 1]
    } else if (!all(exprGroups$sd > 0)) {
      stop('All groups in exprGroups must have standard deviation > 0.')
    }
    if(!'base' %in% colnames(exprGroups)) {
      exprGroups[, base := 0] }}

  if(any(exprGroups[, fracGenes] <= 0)) {
    stop('All groups in exprGroups must have fracGenes > 0.') }

  # Compute a number of genes per group that sum to nGenes.
  exprGroups[, fracGenes := fracGenes / sum(fracGenes)]
  exprGroups[, numGenes := as.integer(fracGenes * nGenes)]
  if(sum(exprGroups[, numGenes]) != nGenes) {
    exprGroups[1L:(nGenes - sum(exprGroups[, numGenes])), numGenes := numGenes + 1]
  }

  return(exprGroups)
}


getTimes = function(timepointsType, interval, nReps, timepoints,
                    nSamplesPerCond, nCond, period) {
  if (timepointsType == 'auto') {
    tt = (2 * pi / period) * interval *
      0:(period %/% interval - (period %% interval == 0))
    tt = rep(tt, each = nReps)
    times = matrix(rep(tt, each = nCond), ncol = length(tt))
  } else if (timepointsType == 'specified') {
    tt = (2 * pi / period) * timepoints # don't do %%, let rhyFunc figure it out
    times = matrix(rep(tt, each = nCond), ncol = length(tt))
  } else if (timepointsType == 'random') {
    tt = stats::runif(nSamplesPerCond * nCond, min = 0, max = 2 * pi)
    tt = matrix(tt, nrow = nCond)
    times = t(apply(tt, 1, sort))
  } else {
    stop("timepointsType must be 'auto', 'specified', or 'random'.")
  }
  return(times)
}
