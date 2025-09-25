#' @title Early Expression Pattern
#' @description Generate an ideal early expression pattern for S developmental stages.
#' 
#' @param S Number of developmental stages
#' 
#' @return Numeric vector representing early expression pattern
#' 
#' @details
#' Creates a pattern where expression starts low, increases in early stages,
#' and remains high in later stages.
#' 
#' @keywords internal
early_gene <- function(S) {
    c(rep(-1, length.out=S%/%4),
      seq(from=-1, to=1, length.out=S-2*(S%/%4)),
      rep(1, length.out=S%/%4))
}

#' @title Mid Expression Pattern
#' @description Generate an ideal mid expression pattern for S developmental stages.
#' 
#' @param S Number of developmental stages
#' 
#' @return Numeric vector representing mid expression pattern
#' 
#' @details
#' Creates a pattern where expression starts low, peaks in middle stages,
#' and returns to low in later stages.
#' 
#' @keywords internal
mid_gene <- function(S) {
    c(rep(-1, length.out=S%/%6),
      seq(from=-1, to=1, length.out=S%/%4),
      rep(1, length.out=S-2*(S%/%6)-2*(S%/%4)),
      seq(from=1, to=-1, length.out=S%/%4),
      rep(-1, length.out=S%/%6))
}

#' @title Late Expression Pattern
#' @description Generate an ideal late expression pattern for S developmental stages.
#' 
#' @param S Number of developmental stages
#' 
#' @return Numeric vector representing late expression pattern
#' 
#' @details
#' Creates a pattern where expression starts high and decreases in later stages.
#' This is the inverse of the early pattern.
#' 
#' @keywords internal
late_gene <- function(S) {
    -early_gene(S)
}

#' @title Reverse Mid Expression Pattern
#' @description Generate an ideal reverse mid expression pattern for S developmental stages.
#' 
#' @param S Number of developmental stages
#' 
#' @return Numeric vector representing reverse mid expression pattern
#' 
#' @details
#' Creates a pattern where expression starts high, drops in middle stages,
#' and returns to high in later stages. This is the inverse of the mid pattern.
#' 
#' @keywords internal
rev_mid_gene <- function(S) {
    -mid_gene(S)
}

#' @title Modulo Pi Function
#' @description Wrap angles to the range (-pi, pi].
#' 
#' @param x Numeric vector of angles
#' 
#' @return Numeric vector of wrapped angles
#' 
#' @details
#' Ensures angles are in the standard range (-pi, pi] by wrapping around 2pi.
#' 
#' @keywords internal
mod_pi <- function(x) {
    (x + pi) %% (2*pi) - pi
}

#' @title Calculate Gene Expression Angles
#' @description Calculate developmental stage angles for genes based on their expression patterns.
#' 
#' @param e Matrix of expression values with genes as rows and stages as columns
#' 
#' @return Numeric vector of angles representing each gene's expression pattern
#' 
#' @details
#' This function uses PCA to project genes and ideal expression patterns into 2D space,
#' then calculates the angle of each gene relative to ideal patterns. The angles
#' represent the peak expression timing during development.
#' 
#' @keywords internal
get_angles <- function(e) {
    N <- nrow(e)
    S <- ncol(e)
    ideal_genes <- rbind(early_gene(S), mid_gene(S), 
                         late_gene(S), rev_mid_gene(S))
    rownames(ideal_genes) <- c("early", "mid", "late", "rev_mid")
    e_aug <- rbind(e, ideal_genes)
    
    pca <- stats::prcomp(e_aug, scale. = FALSE)
    coords <- pca$x[, 1:2]
    
    df <- data.frame(PC1 = coords[,1], PC2 = coords[,2])
    
    df$angle <- atan2(df$PC2, df$PC1)
    
    df$angle <- mod_pi(df$angle - df["rev_mid", "angle"] + pi)
    
    if (df["early", "angle"] < df["late", "angle"])
        df$angle <- mod_pi(-df$angle)
    
    df[1:N, "angle"]
}
