#
#
# defaultGurobiSolveCarnivalOptions <- function(...){
#   if ( "solver" %in% names(list(...)) ) {
#     stop("Don't try to redefine solver in this function.",
#          " Use other default functions for other solvers.")
#   }
#
#   options <- list(
#     solverPath = '',
#     solver = 'gurobi',
#     mipGap = 0.05,
#     poolrelGap = 0.0001,
#     limitPop = 500,
#     poolCap = 100,
#     threads = 1,
#     outputFolder = getwd(),
#     betaWeight = 0.2,
#     cleanTmpFiles = TRUE,
#     keepLPFiles = TRUE,
#     workdir = getwd())
#
#
# }
#
# carnivalOptions <- default_CARNIVAL_options('cplex')
#
# resultFile <- carnivalOptions$resultFile
#
# carnivalOptions <- CARNIVAL::setCarnivalOptions('cplex')
#
# carnivalOptions
#
# CARNIVAL::getOptionsList()
#
#
# paste0(carnivalOptions$solverPath,
#        " MIPGAP=", carnivalOptions$mipGap,
#        " TimeLimit=", carnivalOptions$timelimit,
#        " PoolGap=", carnivalOptions$poolrelGap,
#        " SolutionLimit=", carnivalOptions$limitPop,
#        " PoolSolutions=", carnivalOptions$poolCap,
#        " Threads=", carnivalOptions$threads,
#        " SolFiles=", resultFile,
#        " PoolSearchMode=2 NumericFocus=2",
#        " ResultFile=", paste0(resultFile, ".sol"),
#        additional_params,
#        " ", lpFile)
#
# carnivalOptions$solverPath <- '/home/lorenzo/opt/gurobi950/linux64/bin/python3.7 /home/lorenzo/opt/gurobi950/linux64/lib/gurobi.py'
#
# carnivalOptions$solver <- 'gurobi'
# additional_params <- NULL
# lpFile <- "lpFile_t11_57_44d14_03_2023n93.lp"
# paste0(carnivalOptions$solverPath,
# " MIPGAP=", carnivalOptions$mipGap,
# " TimeLimit=", carnivalOptions$timelimit,
# " PoolGap=", carnivalOptions$poolrelGap,
# " SolutionLimit=", carnivalOptions$limitPop,
# " PoolSolutions=", carnivalOptions$poolCap,
# " Threads=", carnivalOptions$threads,
# " SolFiles=", resultFile,
# " PoolSearchMode=2 NumericFocus=2",
# " ResultFile=", paste0(resultFile, ".sol"),
# additional_params,
# " ", lpFile)
