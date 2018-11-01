defaultDispFunc = function(x) {
  return(3/x)
}

setDefaultExprGroups = function(exprGroups, nGenes, rhyFunc, method) {
  exprGroups = data.table(exprGroups)
  exprGroups[, group := 1:.N]

  if(nrow(exprGroups) == 0) {
    stop('exprGroups must have at least one row.') }

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

  return(exprGroups)
}


getTimes = function(timepointsType, interval, nReps, timepoints,
                    nSamplesPerCond, nCond, period) {
  if (timepointsType == 'auto') {
    tt = interval * 0:(period %/% interval - (period %% interval == 0))
    tt = rep(tt, each = nReps)
    times = matrix(rep(tt, each = nCond), ncol = length(tt))
  } else if (timepointsType == 'specified') {
    # don't do %%, let rhyFunc figure it out
    times = matrix(rep(timepoints, each = nCond), ncol = length(timepoints))
  } else if (timepointsType == 'random') {
    tt = stats::runif(nSamplesPerCond * nCond, min = 0, max = period)
    tt = matrix(tt, nrow = nCond)
    times = t(apply(tt, 1, sort))
  } else {
    stop("timepointsType must be 'auto', 'specified', or 'random'.")
  }
  return(times)
}
