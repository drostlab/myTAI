

transform_counts <- function(phyex_set, FUN) {
    phyex_set@count_matrix <- FUN(phyex_set@count_matrix)
}


