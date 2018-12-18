setDefaultExprGroups = function(exprGroups, nFeatures, dispFunc, rhyFunc, family) {
  if ('group' %in% colnames(exprGroups)) {
    stop("exprGroups must not have a column named 'group'.")}

  exprGroups = data.table(exprGroups)
  exprGroups[, group := 1:.N]

  if (nrow(exprGroups) == 0) {
    stop('exprGroups must have at least one row.') }

  if (!'amp' %in% colnames(exprGroups)) {
    exprGroups[, amp := 0]}

  if (!'phase' %in% colnames(exprGroups)) {
    exprGroups[, phase := 0]}

  if (!'rhyFunc' %in% colnames(exprGroups)) {
    exprGroups[, rhyFunc := data.table(rhyFunc)]}

  if (family == 'negbinom') {
    if (!'dispFunc' %in% colnames(exprGroups)) {
      exprGroups[, dispFunc := data.table(dispFunc)] }
    if (!'base' %in% colnames(exprGroups)) {
      exprGroups[, base := 8]}}
  else {
    if (!'sd' %in% colnames(exprGroups)) {
      exprGroups[, sd := 1]
    } else if (!all(exprGroups$sd > 0)) {
      stop('All groups in exprGroups must have standard deviation > 0.')}
    if (!'base' %in% colnames(exprGroups)) {
      exprGroups[, base := 0]}}

  return(exprGroups)}


getTimes = function(timepointsType, interval, nReps, timepoints,
                    nSamplesPerCond, nConds, period) {
  if (timepointsType == 'auto') {
    tt = interval * 0:(period %/% interval - (period %% interval == 0))
    tt = rep(tt, each = nReps)
    times = matrix(rep(tt, each = nConds), ncol = length(tt))
  } else if (timepointsType == 'specified') {
    if (is.null(timepoints)) {
      stop("timepoints cannot be NULL, if timepointsType is 'specified'.")}
    # don't do %%, let rhyFunc figure it out
    times = matrix(rep(timepoints, each = nConds), ncol = length(timepoints))
  } else if (timepointsType == 'random') {
    if (is.null(nSamplesPerCond)) {
      stop("nSamplesPerCond cannot be NULL, if timepointsType is 'random'.")}
    tt = stats::runif(nSamplesPerCond * nConds, min = 0, max = period)
    tt = matrix(tt, nrow = nConds)
    times = t(apply(tt, 1, sort))}
  return(times)}


getSampleMetadata = function(times) {
  nSamplesPerCond = ncol(times)
  nSamples = prod(dim(times))
  sm = foreach(cond = 1:nrow(times), .combine = rbind) %do% {
    sampleIds = ((cond - 1) * nSamplesPerCond + 1):(nSamplesPerCond * cond)
    sampleNames = sprintf(sprintf('sample_%%0%dd', floor(log10(nSamples)) + 1),
                          sampleIds)
    data.table(sample = sampleNames, cond = cond, time = times[cond, ])}
  return(sm)}


getNFeaturesPerGroup = function(exprGroups, fracFeatures, nFeatures) {
  if ('fracFeatures' %in% colnames(exprGroups)) {
    fracFeatures = exprGroups$fracFeatures
  } else if (is.null(fracFeatures)) {
    fracFeatures = rep(1 / nrow(exprGroups), nrow(exprGroups))
  } else if (length(fracFeatures) != nrow(exprGroups)) {
    stop('Length of fracFeatures must equal number of rows in exprGroups.')}

  nFeaturesPerGroup = as.integer(fracFeatures * nFeatures)
  if (sum(nFeaturesPerGroup) != nFeatures) {
    nFeaturesPerGroup[1L:(nFeatures - sum(nFeaturesPerGroup))] =
      nFeaturesPerGroup[1L:(nFeatures - sum(nFeaturesPerGroup))] + 1L}

  if (any(nFeaturesPerGroup) == 0) {
    stop(paste(c('At least one group has no features. Increase nFeatures,',
                 'reduce the number of groups, or change fracFeatures.'),
               collapse = ' '))}
  return(nFeaturesPerGroup)}


getFeatureMetadata = function(exprGroupsList, fracFeatures, nFeatures) {
  nFeaturesPerGroup = getNFeaturesPerGroup(exprGroupsList[[1]], fracFeatures, nFeatures)
  features = sprintf(sprintf('feature_%%0%dd', floor(log10(nFeatures)) + 1), 1:nFeatures)
  nConds = length(exprGroupsList)

  fm = foreach(exprGroups = exprGroupsList, cond = 1:nConds, .combine = rbind) %do% {
    fmNow = exprGroups[rep(1:.N, times = nFeaturesPerGroup)]
    fmNow[, cond := ..cond]
    fmNow[, feature := features]
    data.table::setcolorder(fmNow, c('cond', 'group', 'feature'))
    fmNow}
  return(fm)}
