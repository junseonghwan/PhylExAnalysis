# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @export
ReadParentVector <- function(file_path, mutation_count) {
    .Call(`_PhylExR_ReadParentVector`, file_path, mutation_count)
}

#' @export
GetChains <- function(parent_vector) {
    .Call(`_PhylExR_GetChains`, parent_vector)
}

get_parent_name <- function(node) {
    .Call(`_PhylExR_get_parent_name`, node)
}

#' @export
GetConfigMatrix <- function(datum2node, ordered_nodes) {
    .Call(`_PhylExR_GetConfigMatrix`, datum2node, ordered_nodes)
}

LogAdd <- function(x, y) {
    .Call(`_PhylExR_LogAdd`, x, y)
}

LogSumExp <- function(x) {
    .Call(`_PhylExR_LogSumExp`, x)
}

IdentifyNodeMutationStatus <- function(datum2node, ordered_nodes, ordered_mutations) {
    .Call(`_PhylExR_IdentifyNodeMutationStatus`, datum2node, ordered_nodes, ordered_mutations)
}

