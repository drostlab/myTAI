
re.colors <- function(n)
{
        
#         colos <- c("black","red","green","brown","darkmagenta",
#                    "blue","darkred","darkblue","darkgreen", "orange",
#                    "azure4","gold4","greenyellow","hotpink4",
#                    "mediumorchid3","mediumorchid3","peachpuff4",
#                    "hotpink","lightgoldenrod", "peru", "slateblue3", "yellow4", "yellowgreen")
        
        return(grDevices::colorRampPalette(RColorBrewer::brewer.pal(8,"Dark2"))(n) )
        
}



PS_colours <- function(n) {
    vals <- 1:n |>
        log() |>
        scales::rescale()
   
    pal <- colorRampPalette(c("black", "#AD6F3B", "lightgreen"))
   
    pal(100)[floor(vals * 99) +1]
}   