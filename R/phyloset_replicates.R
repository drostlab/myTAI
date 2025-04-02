
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
        )
    ),
    constructor = function(data, 
                           groups = NULL,
                           name = deparse(substitute(data)),
                           index_type = "TXI", 
                           conditions_label = "Ontogeny", 
                           is_time_series = TRUE,
                           bootstrap_sample_size = 5000L, 
                           null_conservation_sample_size = 5000L) {
        count_matrix_reps <- as.matrix(data[3:ncol(data)])
        
        
        if (is.null(groups))
            groups <- colnames(count_matrix_reps)
        
        print(groups)
                
        count_matrix <- .collapse_replicates(count_matrix_reps, groups)
        
        new_object(PhyloExpressionSet(strata_vector = as.numeric(data[[1]]),
                                      gene_ids = as.character(data[[2]]),
                                      count_matrix = count_matrix,
                                      name = name,
                                      index_type = index_type, 
                                      conditions_label = conditions_label,
                                      is_time_series=is_time_series,
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

collapse <- function(phyex_set) {
   S7::convert(phyex_set, PhyloExpressionSet)
}

.collapse_replicates <- function(count_matrix, groups) {
    sapply(unique(groups), \(g) rowMeans(count_matrix[, groups == g, drop = FALSE]))
}

S7::method(transform_counts, PhyloExpressionSetReplicates) <- function(phyex_set, 
                                                                       FUN,
                                                                       FUN_name=deparse(substitute(FUN)),
                                                                       new_name=paste(phyex_set@name, "transformed by", FUN_name)) {
    f <- match.fun(FUN)
    phyex_set@count_matrix_reps <- f(phyex_set@count_matrix_reps)
    phyex_set@name <- new_name
    return(phyex_set)
}

S7::method(select_genes, PhyloExpressionSetReplicates) <- function(phyex_set, 
                                                                   genes) {
    indices <- (phyex_set@gene_ids %in% genes)
    
    phyex_set@strata_vector <- phyex_set@strata_vector[indices]
    phyex_set@gene_ids <- phyex_set@gene_ids[indices]
    phyex_set@count_matrix_reps <- phyex_set@count_matrix_reps[indices, ]
    
    return(phyex_set)
}

