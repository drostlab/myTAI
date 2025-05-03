#' 
#' 
#' #' @import S7
#' PhyloExpressionSetMasked <- new_class("PhyloExpressionSetMasked",
#'     parent=PhyloExpressionSet,
#'     properties = list(
#'         ## CONSTRUCTOR PARAMETERS
#'         full_set = new_required_property( # the full phyex set 
#'             class = PhyloExpressionSet,
#'             validator = \(value) if(!S7_inherits(value, PhyloExpressionSet)) "must be a PhyloExpressionSet object",
#'             name = "full_set",
#'         ),
#'         removed_genes = new_property(
#'             class = class_character,
#'             default = character(0)
#'         ),
#'         ## FIELDS & PROPERTIES
#'         genes_bitvector = new_property(
#'             class = class_logical,
#'             getter = \(self) ! self@gene_ids %in% self@removed_genes
#'         ),
#'         # reuse from full set
#'         # bootstrapped_txis = new_property(
#'         #     #class = class_matrix, # S7 doesn't support class_matrix yet
#'         #     getter = \(self) self@full_set@bootstrapped_txis
#'         # ),
#'         null_conservation_txis = new_cached_property(
#'             #class = class_matrix, # S7 doesn't support class_matrix yet
#'             getter = \(self) self@full_set@null_conservation_txis
#'         )
#'         
#'     ),
#'     validator = function(self) {
#'         if (any(! self@removed_genes %in% self@full_set@gene_ids))
#'             "@removed_genes must be a subset of @full_set@gene_ids"
#'     },
#'     constructor = function (full_set,
#'                             removed_genes,
#'                             name=paste(full_set@name, "masked"),
#'                             ...) {
#'         new_object(remove_genes(PhyloExpressionSet, removed_genes, ...),
#'                    full_set = full_set, removed_genes = removed_genes)
#'     }
#' )
#' 
#' 
#' S7::method(print, PhyloExpressionSetMasked) <- function(x, ...) {
#'     cat("\n", "Phylo Expression Set (masked)", "\n", sep="")
#'     cat(x@conditions_label, ":  ", paste(as.character(x@conditions), collapse = ", "), "\n", sep="")
#'     cat("Number of genes: ", x@full_set@num_genes, ", out of which ", length(x@removed_genes), " removed\n" ,sep="")
#'     cat("Removed genes: ", paste(x@removed_genes, collapse = ", "), "\n", sep="")
#' }
