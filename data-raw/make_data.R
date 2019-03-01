dispFuncs = readRDS('../../projects/simphony_analysis/results/mouse_disp_funcs.rds')
locfitDispFunc = dispFuncs$local

##### fit a smoothing spline to the locfit curve to speed up load time

makeDefaultDispFunc = function(xLoc, yLoc) {
  ss = smooth.spline(xLoc, yLoc)
  defaultDispFunc = function(x) {
    if (any(x <= 0)) {
      stop('Mean counts must be greater than zero.')}
    y = 10^predict(ss, log10(x))$y
    return(y)}
  return(defaultDispFunc)}

x = seq(-1, 7, 0.01)
y = log10(locfitDispFunc(10^x))
defaultDispFunc = makeDefaultDispFunc(x, y)

plot(x, y, pch = '.')
lines(x, predict(ss)$y, col = 'blue')
lines(x, log10(defaultDispFunc(10^x)), col = 'red')

usethis::use_data(defaultDispFunc, internal = FALSE, overwrite = TRUE)

##### cut a little fat
# env = environment(defaultDispFunc)
#
# ls(env)
# rm('means', 'disps', 'd', 'minDisp', envir = env)
#
# # names(env$fit)
# # env$fit$frame = NULL
# # failed: eva, cell, terms, nvc, mi, , trans
# # ok to remove: box, sty, deriv,
#
# usethis::use_data(defaultDispFunc, internal = FALSE, overwrite = TRUE)
#
##### check that it works

load('data/defaultDispFunc.rda')
x = seq(-1, 7, 0.01)
y = log10(defaultDispFunc(10^x))
plot(x, y, pch = '.')
