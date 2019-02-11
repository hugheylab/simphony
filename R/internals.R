setDefaultFeatureGroups = function(featureGroups, nFeatures, dispFunc, rhyFunc,
                                   family) {
  if ('group' %in% colnames(featureGroups)) {
    stop("featureGroups must not have a column named 'group'.")}

  featureGroups = data.table(featureGroups)
  featureGroups[, group := 1:.N]

  if (nrow(featureGroups) == 0) {
    stop('featureGroups must have at least one row.') }

  if (!'amp' %in% colnames(featureGroups)) {
    featureGroups[, amp := list(list(function(x) 0))]
  } else {
    if (is.numeric(featureGroups[, amp])) {
      makefunc = function(x) { x; function(m) x }
      if (nrow(featureGroups) == 1) {
        featureGroups[, amp := list(list(makefunc(amp)))]
      } else {
        featureGroups[, amp := foreach(a = featureGroups[, amp]) %do% {
                                           makefunc(a) }]
      }
    }
    else if (sum(unlist(lapply(featureGroups[, amp], function(x) !is.function(x)))) == 0) {}
    else { stop('Amplitudes must be either numeric values or functions')}
  }

  featureGroups[, amp0 := sapply(featureGroups[, amp], mapply, 0)]

  if (!'phase' %in% colnames(featureGroups)) {
    featureGroups[, phase := 0]}

  if (!'rhyFunc' %in% colnames(featureGroups)) {
    featureGroups[, rhyFunc := data.table(rhyFunc)]}

  if (family == 'negbinom') {
    if (!'dispFunc' %in% colnames(featureGroups)) {
      featureGroups[, dispFunc := data.table(dispFunc)] }
    if (!'base' %in% colnames(featureGroups)) {
      featureGroups[, base := list(list(function(x) 8))]
    } else {
      if (is.numeric(featureGroups[, base])) {
        makefunc = function(x) { x; function(m) x }
        if (nrow(featureGroups) == 1) {
          featureGroups[, base := list(list(makefunc(base)))]
        } else {
          featureGroups[, base := foreach(b = featureGroups[, base]) %do% {
                                            makefunc(b) }]
        }
      }
      else if (sum(unlist(lapply(featureGroups[, base], function(x) !is.function(x)))) == 0) {}
      else { stop('Bases must be either numeric values or functions')}
    }
  }
  else {
    if (!'sd' %in% colnames(featureGroups)) {
      featureGroups[, sd := 1]
    } else if (!all(featureGroups$sd >= 0)) {
      stop('All groups in featureGroups must have standard deviation >= 0.')}
        if (!'base' %in% colnames(featureGroups)) {
      featureGroups[, base := list(list(function(x) 8))]
    } else {
      if (is.numeric(featureGroups[, base])) {
        makefunc = function(x) { x; function(m) x }
        if (nrow(featureGroups) == 1) {
          featureGroups[, base := list(list(makefunc(base)))]
        } else {
          featureGroups[, base := foreach(b = featureGroups[, base]) %do% {
                                            makefunc(b) }]
        }
      }
      else if (sum(unlist(lapply(featureGroups[, base], function(x) !is.function(x)))) == 0) {}
      else { stop('Bases must be either numeric values or functions')}
    }
  }
  featureGroups[, base0 := sapply(featureGroups[, base], mapply, 0)]

  return(featureGroups)}


getTimes = function(timepointsType, interval, nReps, timepoints, timeRange,
                    nSamplesPerCond, nConds) {
  if (timepointsType == 'auto') {
    tt = interval * (timeRange[1] %/% interval - (timeRange[1] %% interval == 0) + 1):(timeRange[2] %/% interval - (timeRange[2] %% interval == 0))
    tt = rep(tt, each = nReps)
    times = matrix(rep(tt, each = nConds), ncol = length(tt))
  } else if (timepointsType == 'specified') {
    if (is.null(timepoints)) {
      stop("timepoints cannot be NULL, if timepointsType is 'specified'.")}
    times = matrix(rep(timepoints, each = nConds), ncol = length(timepoints))
  } else if (timepointsType == 'random') {
    if (is.null(nSamplesPerCond)) {
      stop("nSamplesPerCond cannot be NULL, if timepointsType is 'random'.")}
    tt = stats::runif(nSamplesPerCond * nConds, min = timeRange[1], max = timeRange[2])
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


getNFeaturesPerGroup = function(featureGroups, fracFeatures, nFeatures) {
  if ('fracFeatures' %in% colnames(featureGroups)) {
    fracFeatures = featureGroups$fracFeatures
  } else if (is.null(fracFeatures)) {
    fracFeatures = rep(1 / nrow(featureGroups), nrow(featureGroups))
  } else if (length(fracFeatures) != nrow(featureGroups)) {
    stop('Length of fracFeatures must equal number of rows in featureGroups.')}

  nFeaturesPerGroup = as.integer(fracFeatures * nFeatures)
  if (sum(nFeaturesPerGroup) != nFeatures) {
    nFeaturesPerGroup[1L:(nFeatures - sum(nFeaturesPerGroup))] =
      nFeaturesPerGroup[1L:(nFeatures - sum(nFeaturesPerGroup))] + 1L}

  if (any(nFeaturesPerGroup) == 0) {
    stop(paste(c('At least one group has no features. Increase nFeatures,',
                 'reduce the number of groups, or change fracFeatures.'),
               collapse = ' '))}
  return(nFeaturesPerGroup)}


getFeatureMetadata = function(featureGroupsList, fracFeatures, nFeatures) {
  nFeaturesPerGroup = getNFeaturesPerGroup(featureGroupsList[[1]], fracFeatures,
                                           nFeatures)
  features = sprintf(sprintf('feature_%%0%dd', floor(log10(nFeatures)) + 1),
                     1:nFeatures)
  nConds = length(featureGroupsList)

  fm = foreach(featureGroups = featureGroupsList, cond = 1:nConds, .combine = rbind) %do% {
    fmNow = featureGroups[rep(1:.N, times = nFeaturesPerGroup)]
    fmNow[, cond := ..cond]
    fmNow[, feature := features]
    data.table::setcolorder(fmNow, c('cond', 'group', 'feature'))
    fmNow}
  return(fm)}
