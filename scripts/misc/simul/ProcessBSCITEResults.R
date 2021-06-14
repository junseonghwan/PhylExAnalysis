library(Rcpp)

source("R/EvaluationFunctions.R")
source("R/FindClones.R")
sourceCpp("src/rcpp_hello_world.cpp")

simul_path <- "/Users/seonghwanjun/data/simul"
scenarios <- c("binary", "binary_cn", "quadternary_cn")
for (scenario in scenarios) {
    data_path <- paste(simul_path, scenario, sep="/")
    GetBSciteResults(data_path, cases = 1:3, rep_end = 20, is_multiregion = FALSE)
}

scenarios_multiregion <- c("quadternary_multiregion", "quadternary_cn_multiregion")
for (scenario in scenarios_multiregion) {
    data_path <- paste(simul_path, scenario, sep="/")
    GetBSciteResults(data_path, cases = 1:3, rep_end = 20, is_multiregion = TRUE)
}
