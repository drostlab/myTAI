#' @title Generate Phylostratum Colors
#' @description Generate a color palette for phylostrata visualization using a log-scaled transformation.
#' @param n number of colors to generate
#' @return A character vector of color codes
#' @examples
#' # Generate colors for 5 phylostrata
#' colors <- PS_colours(5)
#' @importFrom grDevices colorRampPalette
#' @export
PS_colours <- function(n) {
    vals <- 1:n |>
        log() |>
        scales::rescale()
   
    pal <- grDevices::colorRampPalette(c("black", "#AD6F3B", "lightgreen"))
   
    pal(100)[floor(vals * 99) +1]
}
