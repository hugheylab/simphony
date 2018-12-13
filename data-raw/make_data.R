dispFuncs = readRDS('../../projects/simphony_analysis/results/mouse_disp_funcs.rds')
defaultDispFunc = dispFuncs$local

##### cut a little fat
env = environment(defaultDispFunc)

ls(env)
rm('means', 'disps', 'd', 'minDisp', envir = env)

# names(env$fit)
# env$fit$frame = NULL
# failed: eva, cell, terms, nvc, mi, , trans
# ok to remove: box, sty, deriv,

usethis::use_data(defaultDispFunc, internal = FALSE, overwrite = TRUE)

##### check that it still works

load('data/defaultDispFunc.rda')
defaultDispFunc(2^(0:20))
