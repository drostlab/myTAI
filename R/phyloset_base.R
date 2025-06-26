
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
        stratas = new_required_property(
            class = class_factor,
            validator = function(value) {
                if (any(is.na(value))) "cannot contain NA values. Check data[1]."
                if (length(value) == 0) "cannot be empty. Check data[1]"
                },
            name = "stratas"
        ),
        gene_ids = new_required_property(
            class = class_character,
            validator = function(value) {
                if (any(is.na(value))) "cannot contain NA values. Check data[2]."
                if (length(value) == 0) "cannot be empty. Check data[2]"
            },
            name = "gene_ids"
        ),
        counts = new_required_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            validator = function(value) {
                if (any(is.na(value))) "cannot contain NA values. Check data[3:ncol(data)]."
                if (length(value) == 0) "cannot be empty. Check data[3:ncol(data)]"
            },
            name = "counts"
        ),
        groups = new_required_property(
            class = class_factor,
            name = "groups"
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
        counts_collapsed = new_cached_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            getter = \(self) .collapse_replicates(self@counts, self@groups),
            validator = \(value) if (any(is.na(value))) "cannot contain NA values. Check data[3:ncol(data)].",
        ),
        data = new_property(
            class = class_data.frame,
            getter <- \(self) tibble::tibble(Stratum=self@stratas, 
                                             GeneID=self@gene_ids, 
                                             tibble::as_tibble(self@counts))
        ),
        data_collapsed = new_property(
            class = class_data.frame,
            getter <- \(self) tibble::tibble(Stratum=self@stratas, 
                                             GeneID=self@gene_ids, 
                                             tibble::as_tibble(self@counts_collapsed))
        ),
        index_full_name = new_property(
            class = class_character,
            getter <- \(self) TI_map[[self@index_type]]
        ),
        conditions = new_property(
            class = class_factor,
            getter = \(self) factor(colnames(self@counts_collapsed), 
                                    levels=unique(colnames(self@counts_collapsed),
                                    ordered=TRUE))
            ),
        num_genes = new_property(
            class = class_integer,
            getter = \(self) nrow(self@counts)
            ),
        num_conditions = new_property(
            class = class_integer,
            getter = \(self) length(self@conditions)
            ),
        num_stratas = new_property(
            class = class_integer,
            getter = \(self) length(unique(self@stratas))
        ),
        pTXI = new_cached_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            getter = function(self) {
                m <- pTXI(self@counts_collapsed, self@stratas)
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
                                                             self@counts_collapsed,
                                                             self@bootstrap_sample_size)
            ),
        precomputed_null_conservation_txis = new_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            default = NULL
        ),
        null_conservation_txis = new_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            getter = function(self) {
                if(is.null(self@precomputed_null_conservation_txis))
                   memo_generate_conservation_txis(self@stratas,
                                                   self@counts_collapsed,
                                                   self@null_conservation_sample_size)
                else
                   self@precomputed_null_conservation_txis
            }
        ),
        sample_names = new_property(
            class = class_character,
            getter = function(self) colnames(self@counts)
        ),
        group_map = new_property(
            class = class_list,
            getter = function(self) split(self@sample_names, self@groups)
        ),
        pTXI_reps = new_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            getter = function(self) {
                m <- pTXI(self@counts, self@stratas)
                rownames(m) <- self@gene_ids
                colnames(m) <- self@sample_names
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
                                  groups = colnames(data[,3:ncol(data)]),
                                  name = deparse(substitute(data)),
                                  ...) {
    gene_ids = as.character(data[[2]])
    stratas = factor(as.numeric(data[[1]]), levels=sort(unique(as.numeric(data[[1]]))))
    names(stratas) = gene_ids
    
    groups = factor(groups, levels=unique(groups))
    
    counts = as.matrix(data[3:ncol(data)])
    rownames(counts) = gene_ids
    
    return(PhyloExpressionSet(
        stratas = stratas,
        gene_ids = gene_ids,
        counts = counts,
        groups = groups,
        name = name,
        ...
    ))
}

#' @import dplyr
match_map <- function(data, 
                      phylomap,
                      groups = colnames(data[,2:ncol(data)]),
                      name = NULL,
                      ...) {
    if (is.null(name)) name <- deparse(substitute(data))
    colnames(phylomap) <- c("Stratum", "GeneID")
    
    data <- data |>
        inner_join(phylomap, by="GeneID") |>
        relocate(Stratum, .before = 1)
    
    return(as_PhyloExpressionSet(data, groups=groups, name=name, ...))
}




S7::method(print, PhyloExpressionSet) <- function(x, ...) {
    cat("\n", "Phylo Expression Set", "\n", sep="")
    cat("\n", x@name, "\n", sep="")
    cat(x@conditions_label, ":  ", paste(as.character(x@conditions), collapse = ", "), "\n", sep="")
    #cat("Replicates: ", paste(as.character(x@sample_names), collapse = ", "), "\n", sep="")
    cat("Number of genes: ", x@num_genes, "\n", sep="")
}


.collapse_replicates <- function(counts, groups) {
    m <- do.call(cbind, lapply(unique(groups), \(g) rowMeans(counts[, groups == g, drop = FALSE])))
    colnames(m) <- unique(groups)
    return(m)
}

collapse <- function(phyex_set) {
    data <- tibble::tibble(Stratum=phyex_set@stratas, 
                           GeneID=phyex_set@gene_ids, 
                           tibble::as_tibble(phyex_set@counts_collapsed))
    as_PhyloExpressionSet(data)
}

normalise_stage_expression <- function(phyex_set, total=1e6) {
    phyex_set@counts <- sweep(phyex_set@counts, 2, colSums(phyex_set@counts), FUN="/") * total
    phyex_set
}

transform_counts <- S7::new_generic("transform_counts", "phyex_set")
S7::method(transform_counts, PhyloExpressionSet) <- function(phyex_set, 
                                                             FUN,
                                                             FUN_name=deparse(substitute(FUN)),
                                                             new_name=paste(phyex_set@name, "transformed by", FUN_name)) {
    f <- match.fun(FUN)
    phyex_set@counts <- f(phyex_set@counts)
    #phyex_set@name <- new_name
    return(phyex_set)
}

#' @export
tf <- transform_counts


select_genes <- S7::new_generic("select_genes", "phyex_set")
S7::method(select_genes, PhyloExpressionSet) <- function(phyex_set, 
                                                         genes) {
    indices <- (phyex_set@gene_ids %in% genes)
    
    phyex_set@stratas <- phyex_set@stratas[indices]
    phyex_set@gene_ids <- phyex_set@gene_ids[indices]
    phyex_set@counts <- phyex_set@counts[indices, ,drop=F]
    
    
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


pTXI <- function(counts, stratas) {
    
    sweep(counts, 2, colSums(counts), "/") * as.numeric(stratas)
}


sTXI <- function(phyex_set,
                 option="identity") {
    mat <- rowsum(phyex_set@pTXI, phyex_set@stratas)
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
    s@precomputed_null_conservation_txis <- phyex_set@null_conservation_txis
    return(s)
}








