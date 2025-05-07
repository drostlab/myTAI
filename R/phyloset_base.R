
quantile_rank <- function(x) {
    ranks <- base::rank(x, ties.method = "average")
    (ranks - 0.5) / length(x)
}

TI_map <- list(TXI = "Transcriptomic Index",
               TAI = "Transcriptomic Age Index",
               TDI = "Transcriptomic Divergence Index",
               TPI = "Transcriptomic Polymorphism Index",
               TEI = "Transcriptomic Evolutionary Index"
               )


#' @import S7
PhyloExpressionSet <- new_class("PhyloExpressionSet",
    properties = list(
        ## PARAMETERS
        # REQUIRED
        strata_vector = new_required_property(
            class = class_factor,
            validator = function(value) {
                if (any(is.na(value))) "cannot contain NA values. Check data[1]."
                if (length(value) == 0) "cannot be empty. Check data[1]"
                },
            name = "strata_vector"
        ),
        gene_ids = new_required_property(
            class = class_character,
            validator = function(value) {
                if (any(is.na(value))) "cannot contain NA values. Check data[2]."
                if (length(value) == 0) "cannot be empty. Check data[2]"
            },
            name = "gene_ids"
        ),
        count_matrix_reps = new_required_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            validator = function(value) {
                if (any(is.na(value))) "cannot contain NA values. Check data[3:ncol(data)]."
                if (length(value) == 0) "cannot be empty. Check data[3:ncol(data)]"
            },
            name = "count_matrix"
        ),
        rep_groups = new_required_property(
            class = class_character,
            name = "rep_groups"
        ),
        # OPTIONAL
        name = new_property(
            class = class_character,
            default = "Phylo Expression Set"
        ),
        species = new_property(
            class = class_character,
            default = NULL
        ),
        index_type = new_options_property( # the type of transcriptomic index
            class = class_character,
            options = names(TI_map),
            default = "TXI"
            ),
        conditions_label = new_property(
            class = class_character,
            default = "Ontogeny"
            ),
        is_time_series = new_property(
            class = class_logical,
            default = TRUE
            ),
        bootstrap_sample_size = new_property(
            class = class_numeric,
            default = 5000L
            ),
        null_conservation_sample_size = new_property(
            class = class_numeric,
            default = 5000L
            ),
        ## FIELDS & PROPERTIES
        count_matrix = new_cached_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            getter = \(self) .collapse_replicates(self@count_matrix_reps, self@rep_groups),
            validator = \(value) if (any(is.na(value))) "cannot contain NA values. Check data[3:ncol(data)].",
        ),
        data = new_property(
            class = class_data.frame,
            getter <- \(self) tibble::tibble(Stratum=self@strata_vector, 
                                             GeneID=self@gene_ids, 
                                             tibble::as_tibble(self@count_matrix_reps))
        ),
        index_full_name = new_property(
            class = class_character,
            getter <- \(self) TI_map[[self@index_type]]
        ),
        conditions = new_property(
            class = class_factor,
            getter = \(self) factor(colnames(self@count_matrix), 
                                    levels=unique(colnames(self@count_matrix),
                                    ordered=TRUE))
            ),
        num_genes = new_property(
            class = class_integer,
            getter = \(self) nrow(self@count_matrix)
            ),
        num_conditions = new_property(
            class = class_integer,
            getter = \(self) ncol(self@count_matrix)
            ),
        num_strata = new_property(
            class = class_integer,
            getter = \(self) length(unique(self@strata_vector))
        ),
        pTXI = new_cached_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            getter = function(self) {
                m <- cpp_pMatrix(self@count_matrix, self@strata_vector)
                rownames(m) <- self@gene_ids
                colnames(m) <- self@conditions
                return(m)
            }
            ),
        TXI = new_property(
            class = class_double,
            getter = \(self) colSums(self@pTXI)
            ),
        bootstrapped_txis = new_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            getter = \(self) memo_generate_bootstrapped_txis(self@pTXI,
                                                             self@count_matrix,
                                                             self@bootstrap_sample_size)
            ),
        null_conservation_txis = new_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            getter = \(self) memo_generate_conservation_txis(self@strata_vector,
                                                             self@count_matrix,
                                                             self@null_conservation_sample_size)
        ),
        sample_names = new_property(
            class = class_character,
            getter = function(self) colnames(self@count_matrix_reps)
        ),
        replicate_map = new_property(
            class = class_list,
            getter = function(self) split(self@sample_names, self@rep_groups)
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
        )
    )

as_PhyloExpressionSet <- function(data, 
                                  rep_groups = colnames(data[,3:ncol(data)]),
                                  name = deparse(substitute(data)),
                                  ...) {
    gene_ids = as.character(data[[2]])
    strata_vector = factor(as.numeric(data[[1]]), levels=sort(unique(as.numeric(data[[1]]))))
    names(strata_vector) = gene_ids
    
    count_matrix_reps = as.matrix(data[3:ncol(data)])
    rownames(count_matrix_reps) = gene_ids
    
    return(PhyloExpressionSet(
        strata_vector = strata_vector,
        gene_ids = gene_ids,
        count_matrix_reps = count_matrix_reps,
        rep_groups = rep_groups,
        name = name,
    ))
}


S7::method(print, PhyloExpressionSet) <- function(x, ...) {
    cat("\n", "Phylo Expression Set", "\n", sep="")
    cat("\n", x@name, "\n", sep="")
    cat(x@conditions_label, ":  ", paste(as.character(x@conditions), collapse = ", "), "\n", sep="")
    #cat("Replicates: ", paste(as.character(x@sample_names), collapse = ", "), "\n", sep="")
    cat("Number of genes: ", x@num_genes, "\n", sep="")
}


.collapse_replicates <- function(count_matrix, groups) {
    m <- do.call(cbind, lapply(unique(groups), \(g) rowMeans(count_matrix[, groups == g, drop = FALSE])))
    colnames(m) <- unique(groups)
    return(m)
}

collapse <- function(phyex_set) {
    data <- tibble::tibble(Stratum=phyex_set@strata_vector, 
                           GeneID=phyex_set@gene_ids, 
                           tibble::as_tibble(phyex_set@count_matrix))
    as_PhyloExpressionSet(data)
}

transform_counts <- S7::new_generic("transform_counts", "phyex_set")
S7::method(transform_counts, PhyloExpressionSet) <- function(phyex_set, 
                                                             FUN,
                                                             FUN_name=deparse(substitute(FUN)),
                                                             new_name=paste(phyex_set@name, "transformed by", FUN_name)) {
    f <- match.fun(FUN)
    phyex_set@count_matrix_reps <- f(phyex_set@count_matrix_reps)
    #phyex_set@name <- new_name
    return(phyex_set)
}

select_genes <- S7::new_generic("select_genes", "phyex_set")
S7::method(select_genes, PhyloExpressionSet) <- function(phyex_set, 
                                                         genes) {
    indices <- (phyex_set@gene_ids %in% genes)
    
    phyex_set@strata_vector <- phyex_set@strata_vector[indices]
    phyex_set@gene_ids <- phyex_set@gene_ids[indices]
    phyex_set@count_matrix_reps <- phyex_set@count_matrix_reps[indices, ,drop=F]
    
    
    return(phyex_set)
}

#' @import purrr
TXI_conf_int <- function(phyex_set, 
                         low_q=.025,
                         high_q=.975) {
    CIs <- apply(phyex_set@bootstrapped_txis, 2, quantile, probs=c(low_q, high_q))
    # xlist <- map(phyex_set@replicate_map, \(reps) phyex_set@TXI_reps[reps])
    # stderr <- map_dbl(xlist, \(x) sd(x) / sqrt(length(x)))
    # mean <- map_dbl(xlist, \(x) unname(mean(x)))
    # dfs <- map_dbl(xlist, \(x) length(x) - 1)
    # low <- mean - map2_dbl(stderr, dfs, \(s, d) qt(1 - low_q, df=d) * s)
    # high <- mean + map2_dbl(stderr, dfs, \(s, d) qt(high_q, df=d) * s)
    return(list(low=CIs[1, ], high=CIs[2, ]))
}




sTXI <- function(phyex_set,
                 option="identity") {
    mat <- rowsum(phyex_set@pTXI, phyex_set@strata_vector)
    if (option == "add") {
        mat <- apply(mat,2,cumsum)
    }
    if (! option %in% c("identity", "add"))
        stop("'option' must be one of 'identity' or 'add",
             call. = FALSE
        )
    return(mat)
}

remove_genes <- function(phyex_set, genes, new_name = paste(phyex_set@name, "perturbed")) {
    selected_genes <- setdiff(phyex_set@gene_ids, genes)
    s <- select_genes(phyex_set, selected_genes)
    s@name <- new_name
    return(s)
}







