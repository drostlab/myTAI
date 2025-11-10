## ----echo = FALSE, message = FALSE--------------------------------------------
library(myTAI)
knitr::opts_chunk$set(
  comment = "#>",
  error = FALSE,
  tidy = FALSE,
  width = 750)

## ----message = FALSE, results = FALSE, fig.height=4, fig.width=6, fig.alt="plot_signature function output", dev.args = list(bg = 'transparent'), fig.align='center'----
library(myTAI)
# obtain an example phylo-expression object
data("example_phyex_set")
# plot away!
myTAI::plot_signature(example_phyex_set)  

