#' @title Compute Partial Transcriptome Evolutionary Index (TEI) Strata Values
#' @description
#' This function computes the partial transcriptome evolutionary
#' index (TEI) values combined for each strata.
#'
#' In detail, each gene gets a \emph{TEI contribution profile}.
#'
#' \deqn{TEI_{is} = f_{is} * ps_i}
#'
#' where \eqn{TEI_{is}} is the partial TEI value of gene i,
#' \eqn{f_{is} = e_{is} / \sum e_{is}} and \eqn{ps_i} is the phylostratum of gene i.
#'
#' \eqn{TEI_{is}} values are combined per \eqn{ps}.
#' @param ExpressionSet expression object with rownames as GeneID (dgCMatrix)
#' or standard PhyloExpressionSet object.
#' @param Phylostratum a named vector representing phylostratum per GeneID with
#' names as GeneID (not used if Expression is PhyloExpressionSet).
#' @param split specify number of columns to split
#' @param showprogress boolean if progressbar should be shown
#' @param threads specify number of threads
#' @details The partial TEI values combined per strata give an overall
#' impression of the contribution of each
#' strata to the global \code{\link{TEI}} pattern.
#' @return a numeric matrix storing the summed partial TEI values for each
#' strata.
#' @references
#' Domazet-Loso T. and Tautz D. (2010).
#' \emph{A phylogenetically based transcriptome age index mirrors ontogenetic
#' divergence patterns}. Nature (468): 815-818.
#'
#' Quint M et al. (2012).
#' \emph{A transcriptomic hourglass in plant embryogenesis}.
#' Nature (490): 98-101.
#'
#' Drost HG et al. (2015)
#' Mol Biol Evol. 32 (5): 1221-1231 doi:10.1093/molbev/msv012
#'
#' @examples
#'
#' # reading a standard PhyloExpressionSet
#' data(PhyloExpressionSetExample, package = "myTAI")
#'
#' # computing partial TEI contribution per gene
#' pS <- pStrataTEI(PhyloExpressionSetExample)
#'
#' @export pStrataTEI
#' @author Kristian K Ullrich

pStrataTEI <- function(ExpressionSet,
                       Phylostratum = NULL,
                       split = 100000,
                       showprogress = TRUE,
                       threads = 1) {
  if (methods::is(ExpressionSet, "Matrix")) {
    common_ids <- sort(Reduce(intersect, list(
      rownames(ExpressionSet),
      names(Phylostratum)
    )))
    Phylostratum <-
      Phylostratum[names(Phylostratum) %in% common_ids]
    Phylostratum <- Phylostratum[order(names(Phylostratum))]
    PhylostratumGroups <- sort(unique(Phylostratum))
    if (ncol(ExpressionSet) > split) {
      split_start <- seq(from = 1,
                         to = ncol(ExpressionSet),
                         by = split)
      if (ncol(ExpressionSet) %% split != 0) {
        split_end <- c(seq(
          from = split,
          to = ncol(ExpressionSet),
          by = split
        ),
        ncol(ExpressionSet))
      } else{
        split_end <- seq(from = split,
                         to = ncol(ExpressionSet),
                         by = split)
      }
      if (rev(split_end - split_start)[1] == 0) {
        split_start <- split_start[-length(split_start)]
        split_end[length(split_end) - 1] <-
          rev(split_end)[1]
        split_end <- split_end[-length(split_end)]
      }
      if (showprogress) {
        pb <- utils::txtProgressBar(min = 1,
                             max = length(split_start),
                             style = 3)
      }
      OUT <- NULL
      es <- NULL
      es_ps <- NULL
      for (i in seq_along(split_start)) {
        es <- ExpressionSet[, split_start[i]:split_end[i]]
        es <- es[rownames(es) %in% common_ids, , drop = FALSE]
        es <- es[order(rownames(es)), ]
        es_ps <- methods::as(
          rcpp_pStrata_parallel(es, Phylostratum,
                                PhylostratumGroups, threads),
          "sparseMatrix"
        )
        colnames(es_ps) <- colnames(es)
        rownames(es_ps) <- PhylostratumGroups
        if (showprogress) {
          utils::setTxtProgressBar(pb, i)
        }
        es <- NULL
        OUT <- cbind(OUT, es_ps)
      }
      ps <- OUT
    } else{
      ExpressionSet <- ExpressionSet[rownames(ExpressionSet) %in%
                                       common_ids, , drop = FALSE]
      ExpressionSet <-
        ExpressionSet[order(rownames(ExpressionSet)), ]
      ps <- methods::as(
        rcpp_pStrata_parallel(
          ExpressionSet,
          Phylostratum,
          PhylostratumGroups,
          threads
        ),
        "sparseMatrix"
      )
      colnames(ps) <- colnames(ExpressionSet)
      rownames(ps) <- PhylostratumGroups
    }
  }
  if (methods::is(ExpressionSet, "data.frame") |
      methods::is(ExpressionSet, "tibble")) {
    if (is.ExpressionSet(ExpressionSet)) {
      Phylostratum <- stats::setNames(ExpressionSet$Phylostratum,
                               ExpressionSet$GeneID)
      PhylostratumGroups <- sort(unique(Phylostratum))
      ExpressionSet <-
        methods::as(data.matrix(ExpressionSet[, 3:ncol(ExpressionSet)]), "sparseMatrix")
      rownames(ExpressionSet) <- names(Phylostratum)
    }
    if (ncol(ExpressionSet) > split) {
      split_start <- seq(from = 1,
                         to = ncol(ExpressionSet),
                         by = split)
      if (ncol(ExpressionSet) %% split != 0) {
        split_end <- c(seq(
          from = split,
          to = ncol(ExpressionSet),
          by = split
        ),
        ncol(ExpressionSet))
      } else{
        split_end <- seq(from = split,
                         to = ncol(ExpressionSet),
                         by = split)
      }
      if (rev(split_end - split_start)[1] == 0) {
        split_start <- split_start[-length(split_start)]
        split_end[length(split_end) - 1] <-
          rev(split_end)[1]
        split_end <- split_end[-length(split_end)]
      }
      if (showprogress) {
        pb <- utils::txtProgressBar(min = 1,
                             max = length(split_start),
                             style = 3)
      }
      OUT <- NULL
      es <- NULL
      es_ps <- NULL
      for (i in seq_along(split_start)) {
        es <- ExpressionSet[, split_start[i]:split_end[i]]
        es_ps <- methods::as(
          rcpp_pStrata_parallel(es, Phylostratum,
                                PhylostratumGroups, threads),
          "sparseMatrix"
        )
        colnames(es_ps) <- colnames(es)
        rownames(es_ps) <- PhylostratumGroups
        if (showprogress) {
          utils::setTxtProgressBar(pb, i)
        }
        es <- NULL
        OUT <- cbind(OUT, es_ps)
      }
      ps <- OUT
    } else{
      ps <- methods::as(
        rcpp_pStrata_parallel(
          ExpressionSet,
          Phylostratum,
          PhylostratumGroups,
          threads
        ),
        "sparseMatrix"
      )
      colnames(ps) <- colnames(ExpressionSet)
      rownames(ps) <- PhylostratumGroups
    }
  }
  return(ps)
}
