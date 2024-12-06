library(ggplot2)
library(tibble)
library(dplyr)

data("PhyloExpressionSetExample")

phylo_expr_set <- PhyloExpressionSetExample

S <- 7

expr_mat.raw <- phylo_expr_set[ , 3:(S+2)]

# optionally split dataset into early/mid/late
modules = list(early=1:2, mid=3:5, late=6:7)
expr_mat.3D <- data.frame(Early=rowMeans(expr_mat.raw[modules$early]),
                          Mid=rowMeans(expr_mat.raw[modules$mid]),
                          Late=rowMeans(expr_mat.raw[modules$late]))

expr_mat <- expr_mat.3D
#expr_mat <- expr_mat.raw

# Plot count distribution
df <- as_tibble(expr_mat.raw) |>
    rowid_to_column("Id") |>
    tidyr::pivot_longer(cols=colnames(expr_mat.raw), names_to="Stage", values_to="Count")
ggplot(data=df, aes(x=Count, fill=Stage)) +
    geom_histogram(bins=100)

pMatrixResiduals <- function(expr_mat,
                             transformation = "none",
                             mask = NULL) {
    N <- nrow(expr_mat)
    if(is.null(mask))
        mask = rep(1, N)
    
    # Apply transformation
    transf = switch(transformation,
        "none" = identity,
        "sqrt" = sqrt,
        "log" = \(x) log2(x+1)
    )
    expr_mat[] <- lapply(expr_mat, transf)
    
    expr_mat <- expr_mat * mask
    
    # Normalise the counts, stage/column wise
    expr_mat.normalised <- mapply('/', expr_mat, colSums(expr_mat))
    
    # Centralise the gene profiles/row wise 
    expr_mat.centralised <- expr_mat.normalised - rowMeans(expr_mat.normalised)
    
    # Multiply by gene age
    ps_vec <- phylo_expr_set[, 1]
    expr_mat.aged <- expr_mat.normalised * ps_vec
    
    return(expr_mat.aged)
}

rho_mat <- pMatrixResiduals(expr_mat, transformation="log")


# Compute TAI
TAI.residuals <- colSums(rho_mat)



# Plot all the genes
df <- as_tibble(rho_mat) |> 
    add_column(Group = "Gene residuals", 
               Phylostratum = phylo_expr_set[[1]],
               GeneID = phylo_expr_set[[2]]) |>
    tidyr::pivot_longer(cols=names(TAI.residuals), names_to="Stage", values_to="Value")

ggplot(data=df, 
       aes(x=factor(Stage, levels=unique(Stage)), 
           y=Value, 
           colour=Phylostratum, 
           group=GeneID)) + 
    geom_line(alpha=0.2)

# Plot the genes as vectors in 3D
library(plotly)

rho_mat.3D <- rho_mat

df <- as_tibble(rho_mat.3D) |> 
    add_column(Group = "Gene residuals", 
               Phylostratum = phylo_expr_set[[1]],
               ID = phylo_expr_set[[2]])


plot_ly(df, x=~Early, y=~Mid, z=~Late,
        marker=list(color = ~Phylostratum, colorscale='Blues', size=1.5, alpha=0.5))

# Plot the genes as vectors in 2D
projectResiduals <- function(expr_mat, tai) {
    tai <- unname(tai)
    proj_mat <- rbind(TAI.residuals, pracma::cross(c(1,1,1), TAI.residuals))
    expr_mat.2D <- expr_mat |>
        data.matrix() |>
        t() |>
        (\(x) proj_mat %*% x)() |>
        t() |>
        unname()
    return(expr_mat.2D)
}

rho_mat.2D <- projectResiduals(rho_mat.3D, TAI.residuals)



df <- as_tibble(rho_mat.2D) |>
    add_column(Group = "Gene residuals", 
               Phylostratum = phylo_expr_set[[1]],
               ID = phylo_expr_set[[2]])

ggplot(df, aes(x=V1, y=V2, size=Group)) +
    geom_point(alpha=0.5)
    
# See how much gene values get changed
randomMask <- function(N, k) {
    mask <- c(rep(1, N-k), rep(0, k)) |>
        sample()
    return(mask)
}
N <- nrow(expr_mat)
k <- 1000
set.seed(32)
mask <- randomMask(N, k)
rho_mat.removed_genes <- pMatrixResiduals(expr_mat, transformation="log", mask=mask)

rho_mat.removed_genes.2D <- projectResiduals(rho_mat.removed_genes, TAI.residuals)

df <- as_tibble(rho_mat.2D) |>
    add_column(U1 = rho_mat.removed_genes.2D[, 1], U2 = rho_mat.removed_genes.2D[, 2]) |>
    filter(mask == 1)



p <- ggplot(df, aes(x=V1, y=V2)) +
    geom_segment(aes(xend=U1, yend=U2), arrow=arrow(length=unit(0.1, "cm")), alpha=0.1)
ggplotly(p)


normalisingFactor <- function(expr_mat, 
                              transformation = "none",
                              mask = NULL) {
    N <- nrow(expr_mat)
    if(is.null(mask))
        mask = rep(1, N)
    
    # Apply transformation
    transf = switch(transformation,
                    "none" = identity,
                    "sqrt" = sqrt,
                    "log" = \(x) log2(x+1)
    )
    expr_mat[] <- lapply(expr_mat, transf)
    expr_mat <- expr_mat * mask
    
    return(colSums(expr_mat))
}

genesToMask <- function(all_genes, subset) {
    mask <- as.integer(!(all_genes %in% subset))
    return(mask)
}

genes <- TopVarianceGenes(phylo_expr_set, p=0.8)

mask <- genesToMask(phylo_expr_set$GeneID, genes)

scaling_factor <- normalisingFactor(expr_mat, transformation="log") / normalisingFactor(expr_mat, transformation="log", mask=mask) 
scaling_factor



