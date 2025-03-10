
#' @import S7
PhyloExpressionSetReplicates <- new_class("PhyloExpressionSetReplicates",
    parent = PhyloExpressionSet,
    properties = list(
        ## CONSTRUCTOR PARAMETERS
        groups = new_required_property(
            class = class_character,
            name="groups"
        ),
        count_matrix_reps = new_required_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            validator = \(value) if (any(is.na(value))) "cannot contain NA values. Check data[3:ncol(data)].",
            name="count_matrix_reps"
        ),
        ## FIELDS & PROPERTIES
        sample_names = new_property(
            class = class_character,
            getter = function(self) colnames(self@count_matrix_reps)
        ),
        replicate_map = new_property(
            class = class_list,
            getter = function(self) split(self@sample_names, self@groups)
        ),
        count_matrix = new_cached_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            getter = \(self) .collapse_replicates(self@count_matrix_reps, self@groups),
            validator = \(value) if (any(is.na(value))) "cannot contain NA values. Check data[3:ncol(data)].",
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
            getter = \(self) .generate_bootstrapped_txis_reps(self@pTXI_reps, 
                                                              self@count_matrix_reps,
                                                              self@groups, 
                                                              self@conditions, 
                                                              self@bootstrap_sample_size)
        )
    ),
    constructor = function(data, 
                           groups,
                           index_type = "TXI", 
                           conditions_label = "Ontogeny", 
                           bootstrap_sample_size = 5000L, 
                           null_conservation_sample_size = 5000L) {
        count_matrix_reps <- as.matrix(data[3:ncol(data)])
        count_matrix <- .collapse_replicates(count_matrix_reps, groups)
        new_object(PhyloExpressionSet(data = data.frame(data[1:2], count_matrix),
                                      index_type = index_type, 
                                      conditions_label = conditions_label, 
                                      bootstrap_sample_size = bootstrap_sample_size, 
                                      null_conservation_sample_size = null_conservation_sample_size), 
                   groups = groups, count_matrix_reps=count_matrix_reps)
    }
    )

S7::method(print, PhyloExpressionSetReplicates) <- function(x, ...) {
    cat("\n", "Phylo Expression Set (with replicates)", "\n", sep="")
    cat(x@conditions_label, ":  ", paste(as.character(x@conditions), collapse = ", "), "\n", sep="")
    cat("Replicates: ", paste(as.character(x@sample_names), collapse = ", "), "\n", sep="")
    cat("Number of genes: ", x@num_genes, "\n", sep="")
}

.collapse_replicates <- function(count_matrix, groups) {
    sapply(unique(groups), \(g) rowMeans(count_matrix[, groups == g]))
}

.generate_bootstrapped_txis_reps <- function(pTXI_reps,
                                             count_matrix_reps,
                                             groups,
                                             conditions,
                                             sample_size) {
    N <- nrow(pTXI_reps)
    
    sampled_indices <- matrix(sample(N, size=N*sample_size, replace=TRUE), 
                              nrow=N, ncol=sample_size)
    selection_matrix <- Matrix::sparseMatrix(i = as.vector(sampled_indices),
                                             j = rep(1:sample_size, each = N),
                                             x = 1,
                                             dims= c(N, sample_size))
    
    # this is necessary to renormalise the counts so as to get the average right
    # it's a bit contrived, might be better to just rewrite this in rcpp...
    total_counts_per_bootstrap <- t(count_matrix_reps) %*% selection_matrix |>
        as.matrix() |>
        t()
    total_counts <- colSums(count_matrix_reps)
    
    bootstrap_matrix_reps <- t(pTXI_reps) %*% selection_matrix |> 
        as.matrix() |>
        t() |> # columns: samples. rows: bootstraps
        sweep(2, total_counts, "*")
    
    bootstrap_matrix_reps <- bootstrap_matrix_reps / total_counts_per_bootstrap
    
    random_rep <- function(groups) {
        v <- sapply(unique(groups), \(g) sample(which(groups == g), 1))
        v <- as.integer(1:length(groups) %in% unlist(v))
        return(v)
    }
    rep_mask <- t(replicate(sample_size, random_rep(groups)))

    bootstrap_matrix <- (rep_mask * bootstrap_matrix_reps) |> # different replicates sampled per bootstrap
        (\(m) sapply(unique(groups), \(g) rowSums(m[, groups == g])))() # collapse replicates => columns: conditions. rows: bootstraps

    colnames(bootstrap_matrix) <- conditions
    
    

    return(bootstrap_matrix)
}
