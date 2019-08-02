setDefaultFeatureGroups = function(featureGroups, nFeatures, dispFunc, rhyFunc,
                                   family, defaultAmp = 0, defaultPhase = 0,
                                   defaultPeriod = 24, defaultSd = 1,
                                   defaultBaseGaussian = 0, defaultBaseNegbinom = 8,
                                   defaultBaseBernoulli = 0.5, defaultBasePoisson = 1) {
  if ('group' %in% colnames(featureGroups)) {
    stop("featureGroups must not have a column named 'group'.")}

  featureGroups = data.table(featureGroups)
  featureGroups[, group := 1:.N]

  if (nrow(featureGroups) == 0) {
    stop('featureGroups must have at least one row.')}

  featureGroups = setFuncs(featureGroups, 'amp', defaultAmp)
  featureGroups[, amp0 := amp[[1]](0), by = 1:nrow(featureGroups)]

  if (!'phase' %in% colnames(featureGroups)) {
    featureGroups[, phase := defaultPhase]}

  if (!'period' %in% colnames(featureGroups)) {
    featureGroups[, period := defaultPeriod]}

  if (!'rhyFunc' %in% colnames(featureGroups)) {
    featureGroups[, rhyFunc := data.table(rhyFunc)]}

  if (family == 'gaussian') {
    if (!'sd' %in% colnames(featureGroups)) {
      featureGroups[, sd := defaultSd]
    } else if (!all(featureGroups$sd >= 0)) {
      stop('All groups in featureGroups must have standard deviation >= 0.')}
    featureGroups = setFuncs(featureGroups, 'base', defaultBaseGaussian)
  } else if (family == 'negbinom') {
    if (!'dispFunc' %in% colnames(featureGroups)) {
      if (is.null(dispFunc)) {
        dispFunc = defaultDispFunc}
      featureGroups[, dispFunc := data.table(dispFunc)]}
    featureGroups = setFuncs(featureGroups, 'base', defaultBaseNegbinom)
  } else if (family == 'bernoulli') {
    featureGroups = setFuncs(featureGroups, 'base', defaultBaseBernoulli)
  } else if (family == 'poisson') {
    featureGroups = setFuncs(featureGroups, 'base', defaultBasePoisson)}

  featureGroups[, base0 := base[[1]](0), by = 1:nrow(featureGroups)]
  return(featureGroups)}


setFuncs = function(featureGroups, varName, defaultValue) {
  if (!varName %in% colnames(featureGroups)) {
    featureGroups[, (varName) := list(list(function(x) defaultValue))]
  } else {
    if (is.numeric(featureGroups[, get(varName)])) {
      makefunc = function(x) {x; function(m) x}
      if (nrow(featureGroups) == 1) {
        featureGroups[, (varName) := list(list(makefunc(featureGroups[1, get(varName)])))]
      } else {
        featureGroups[, (varName) := foreach(v = featureGroups[, get(varName)]) %do% {makefunc(v)}]}
    } else if (!all(sapply(featureGroups[, get(varName)], is.function))) {
      stop(sprintf('%s must be numeric or a list of functions.', varName))}}
  return(featureGroups)}


getTimes = function(timepointsType, interval, nReps, timepoints, timeRange,
                    nSamplesPerCond, nConds) {
  if (timepointsType == 'auto') {
    t1 = timeRange[1] %/% interval - (timeRange[1] %% interval == 0) + 1
    t2 = timeRange[2] %/% interval - (timeRange[2] %% interval == 0)
    tt = interval * t1:t2
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
    condNames = sprintf(sprintf('cond_%%0%dd', floor(log10(nrow(times))) + 1), cond)
    data.table(sample = sampleNames, cond = condNames, time = times[cond, ])}
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
    fmNow[, cond := sprintf(sprintf('cond_%%0%dd', floor(log10(nConds)) + 1), ..cond)]
    fmNow[, feature := features]
    data.table::setcolorder(fmNow, c('cond', 'group', 'feature'))
    fmNow}
  return(fm)}
