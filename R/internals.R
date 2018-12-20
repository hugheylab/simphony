setDefaultExprGroups = function(abundGroups, nFeatures, dispFunc, rhyFunc, family) {
  if ('group' %in% colnames(abundGroups)) {
    stop("abundGroups must not have a column named 'group'.")}

  abundGroups = data.table(abundGroups)
  abundGroups[, group := 1:.N]

  if (nrow(abundGroups) == 0) {
    stop('abundGroups must have at least one row.') }

  if (!'amp' %in% colnames(abundGroups)) {
    abundGroups[, amp := 0]}

  if (!'phase' %in% colnames(abundGroups)) {
    abundGroups[, phase := 0]}

  if (!'rhyFunc' %in% colnames(abundGroups)) {
    abundGroups[, rhyFunc := data.table(rhyFunc)]}

  if (family == 'negbinom') {
    if (!'dispFunc' %in% colnames(abundGroups)) {
      abundGroups[, dispFunc := data.table(dispFunc)] }
    if (!'base' %in% colnames(abundGroups)) {
      abundGroups[, base := 8]}}
  else {
    if (!'sd' %in% colnames(abundGroups)) {
      abundGroups[, sd := 1]
    } else if (!all(abundGroups$sd > 0)) {
      stop('All groups in abundGroups must have standard deviation > 0.')}
    if (!'base' %in% colnames(abundGroups)) {
      abundGroups[, base := 0]}}

  return(abundGroups)}


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


getNFeaturesPerGroup = function(abundGroups, fracFeatures, nFeatures) {
  if ('fracFeatures' %in% colnames(abundGroups)) {
    fracFeatures = abundGroups$fracFeatures
  } else if (is.null(fracFeatures)) {
    fracFeatures = rep(1 / nrow(abundGroups), nrow(abundGroups))
  } else if (length(fracFeatures) != nrow(abundGroups)) {
    stop('Length of fracFeatures must equal number of rows in abundGroups.')}

  nFeaturesPerGroup = as.integer(fracFeatures * nFeatures)
  if (sum(nFeaturesPerGroup) != nFeatures) {
    nFeaturesPerGroup[1L:(nFeatures - sum(nFeaturesPerGroup))] =
      nFeaturesPerGroup[1L:(nFeatures - sum(nFeaturesPerGroup))] + 1L}

  if (any(nFeaturesPerGroup) == 0) {
    stop(paste(c('At least one group has no features. Increase nFeatures,',
                 'reduce the number of groups, or change fracFeatures.'),
               collapse = ' '))}
  return(nFeaturesPerGroup)}


getFeatureMetadata = function(abundGroupsList, fracFeatures, nFeatures) {
  nFeaturesPerGroup = getNFeaturesPerGroup(abundGroupsList[[1]], fracFeatures, nFeatures)
  features = sprintf(sprintf('feature_%%0%dd', floor(log10(nFeatures)) + 1), 1:nFeatures)
  nConds = length(abundGroupsList)

  fm = foreach(abundGroups = abundGroupsList, cond = 1:nConds, .combine = rbind) %do% {
    fmNow = abundGroups[rep(1:.N, times = nFeaturesPerGroup)]
    fmNow[, cond := ..cond]
    fmNow[, feature := features]
    data.table::setcolorder(fmNow, c('cond', 'group', 'feature'))
    fmNow}
  return(fm)}
