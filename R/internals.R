setDefaultExprGroups = function(exprGroups, nGenes, dispFunc, rhyFunc, family) {
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
      exprGroups[, base := 7]}}
  else {
    if (!'sd' %in% colnames(exprGroups)) {
      exprGroups[, sd := 1]
    } else if (!all(exprGroups$sd > 0)) {
      stop('All groups in exprGroups must have standard deviation > 0.')}
    if (!'base' %in% colnames(exprGroups)) {
      exprGroups[, base := 0]}}

  return(exprGroups)}


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
    stop("timepointsType must be 'auto', 'specified', or 'random'.")}
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


getNGenesPerGroup = function(exprGroups, fracGenes, nGenes) {
  if ('fracGenes' %in% colnames(exprGroups)) {
    fracGenes = exprGroups$fracGenes
  } else if (is.null(fracGenes)) {
    fracGenes = rep(1 / nrow(exprGroups), nrow(exprGroups))
  } else if (length(fracGenes) != nrow(exprGroups)) {
    stop('Length of fracGenes must equal number of rows in exprGroups.')}

  nGenesPerGroup = as.integer(fracGenes * nGenes)
  if (sum(nGenesPerGroup) != nGenes) {
    nGenesPerGroup[1L:(nGenes - sum(nGenesPerGroup))] =
      nGenesPerGroup[1L:(nGenes - sum(nGenesPerGroup))] + 1L}

  if (any(nGenesPerGroup) == 0) {
    stop(paste(c('At least one group has no genes. Increase nGenes,',
                 'reduce the number of groups, or change fracGenes.'),
               collapse = ' '))}
  return(nGenesPerGroup)}


getGeneMetadata = function(exprGroupsList, fracGenes, nGenes) {
  nGenesPerGroup = getNGenesPerGroup(exprGroupsList[[1]], fracGenes, nGenes)
  genes = sprintf(sprintf('gene_%%0%dd', floor(log10(nGenes)) + 1), 1:nGenes)
  nCond = length(exprGroupsList)

  gm = foreach(exprGroups = exprGroupsList, cond = 1:nCond, .combine = rbind) %do% {
    gmNow = exprGroups[rep(1:.N, times = nGenesPerGroup)]
    gmNow[, cond := ..cond]
    gmNow[, gene := genes]
    data.table::setcolorder(gmNow, c('cond', 'group', 'gene'))
    gmNow}
  return(gm)}
