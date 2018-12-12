dispFuncs = readRDS('../../projects/simphony_analysis/results/mouse_disp_funcs.rds')
defaultDispFunc = dispFuncs$local
usethis::use_data(defaultDispFunc, internal = FALSE, overwrite = TRUE)
