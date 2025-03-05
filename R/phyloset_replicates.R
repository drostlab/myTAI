
#' @import S7
PhyloExpressionSetReplicates <- new_class("PhyloExpressionSetReplicates",
    parent = PhyloExpressionSet,
    properties = list(
        groups = new_required_property(
            class = class_character,
            name="groups"
        ),
        count_matrix_reps = new_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            getter = function(self) {
                m <- self@data[3:ncol(self@data)] |>
                    as.matrix()
                rownames(m) <- self@gene_ids
                
                return(m)
            },
            validator = \(value) if (any(is.na(value))) "cannot contain NA values. Check data[3:ncol(data)]."
        ),
        count_matrix = new_cached_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            getter = function(self) {
                m <- sapply(unique(self@groups), \(g) rowMeans(self@count_matrix_reps[, self@groups == g]))
                return(m)
            },
            validator = \(value) if (any(is.na(value))) "cannot contain NA values. Check data[3:ncol(data)]."
        ),
        pTXI_reps = new_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            getter = function(self) {
                m <- cpp_pMatrix(self@count_matrix_reps, self@strata_vector)
                rownames(m) <- self@gene_ids
                colnames(m) <- colnames(self@count_matrix_reps)
                return(m)
            }
        ),
        TXI_reps = new_property(
            class = class_double,
            getter = \(self) colSums(self@pTXI_reps)
        ),
        bootstrapped_txis = new_cached_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            getter = .generate_bootstrapped_txis_reps
        )
    )
    )


.generate_bootstrapped_txis_reps <- function(phyex_set,
                                             sample_size=NULL) {
    if (is.null(sample_size))
        sample_size <- phyex_set@null_conservation_sample_size
    
    N <- phyex_set@num_genes
    
    sampled_indices <- matrix(sample(N, size=N*sample_size, replace=TRUE), 
                              nrow=N, ncol=sample_size)
    selection_matrix <- Matrix::sparseMatrix(i = as.vector(sampled_indices),
                                             j = rep(1:sample_size, each = N),
                                             x = 1,
                                             dims= c(N, sample_size))
    
    bootstrap_matrix_reps <- t(phyex_set@pTXI_reps) %*% selection_matrix |> 
        as.matrix() |>
        t() # columns: samples. rows: bootstraps 
    
    random_rep <- function(groups) {
        v <- sapply(unique(groups), \(g) sample(which(groups == g), 1))
        v <- as.integer(1:length(groups) %in% unlist(v))
        return(v)
    }
    rep_mask <- t(replicate(sample_size, random_rep(groups)))

    bootstrap_matrix <- (rep_mask * bootstrap_matrix_reps) |> # different replicates sampled per bootstrap
        (\(m) sapply(unique(phyex_set@groups), \(g) rowSums(m[, phyex_set@groups == g])))() # collapse replicates => columns: conditions. rows: bootstraps

    colnames(bootstrap_matrix) <- phyex_set@conditions

    return(bootstrap_matrix)
}
