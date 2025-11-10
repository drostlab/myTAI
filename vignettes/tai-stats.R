## ----message = FALSE, warning = FALSE-----------------------------------------
library(myTAI)

## ----message = FALSE, results = FALSE, fig.height=4, fig.width=6, fig.alt="plot_signature function output", dev.args = list(bg = 'transparent'), fig.align='center'----
data("example_phyex_set_old")
myTAI::plot_signature(example_phyex_set_old, show_p_val = FALSE)

## ----message = FALSE, eval = FALSE--------------------------------------------
# myTAI::stat_flatline_test(example_phyex_set_old, plot_result = TRUE)

## ----message = FALSE, warning = FALSE, echo = FALSE, fig.height=3, fig.width=5, fig.alt="stat_flatline_test function output", dev.args = list(bg = 'transparent'), fig.align='center'----
myTAI::stat_flatline_test(example_phyex_set_old, plot_result = TRUE)

## ----message = FALSE, warning = FALSE-----------------------------------------
res_flt <- myTAI::stat_flatline_test(example_phyex_set_old, plot_result = FALSE)

## ----message = FALSE, eval = FALSE--------------------------------------------
# myTAI::plot_cullen_frey(res_flt)

## ----message = FALSE, warning = FALSE, fig.height=5, fig.width=5, fig.alt="plot_cullen_frey function output for stat_flatline_test", dev.args = list(bg = 'transparent'), fig.align='center'----
myTAI::plot_cullen_frey(res_flt)

## ----message = FALSE, eval = FALSE--------------------------------------------
# myTAI::plot_null_txi_sample(res_flt)

## ----message = FALSE, warning = FALSE, fig.height=4, fig.width=6, fig.alt="plot_null_txi_sample function output", dev.args = list(bg = 'transparent'), fig.align='center'----
myTAI::plot_null_txi_sample(res_flt)

## ----message = FALSE, eval = FALSE--------------------------------------------
# modules <- list(early = 1:3, mid = 4:5, late = 6:7)
# myTAI::stat_reductive_hourglass_test(
#   example_phyex_set_old, plot_result = TRUE,
#   modules = modules)

## ----message = FALSE, warning = FALSE, echo = FALSE, fig.height=3, fig.width=5, fig.alt="stat_reductive_hourglass_test function output", dev.args = list(bg = 'transparent'), fig.align='center'----
modules <- list(early = 1:3, mid = 4:5, late = 6:7)
myTAI::stat_reductive_hourglass_test(
  example_phyex_set_old, plot_result = TRUE,
  modules = modules)

## ----message = FALSE, warning = FALSE-----------------------------------------
res_flt <- myTAI::stat_reductive_hourglass_test(
  example_phyex_set_old, plot_result = FALSE,
  modules = modules)

## ----message = FALSE, warning = FALSE, fig.height=5, fig.width=5, fig.alt="plot_cullen_frey function output for stat_reductive_hourglass_test", dev.args = list(bg = 'transparent'), fig.align='center'----
myTAI::plot_cullen_frey(res_flt)

## ----message = FALSE, eval = FALSE--------------------------------------------
# modules <- list(early = 1:3, mid = 4:5, late = 6:7)
# myTAI::stat_reverse_hourglass_test(
#   example_phyex_set_old, plot_result = TRUE,
#   modules = modules)

## ----message = FALSE, warning = FALSE, echo = FALSE, fig.height=3, fig.width=5, fig.alt="stat_reverse_hourglass_test function output", dev.args = list(bg = 'transparent'), fig.align='center'----
modules <- list(early = 1:3, mid = 4:5, late = 6:7)
myTAI::stat_reverse_hourglass_test(
  example_phyex_set_old, plot_result = TRUE,
  modules = modules)

## ----message = FALSE, eval = FALSE--------------------------------------------
# modules <- list(early = 1:3, mid = 4:5, late = 6:7)
# myTAI::stat_early_conservation_test(
#   example_phyex_set_old, plot_result = TRUE,
#   modules = modules)

## ----message = FALSE, warning = FALSE, echo = FALSE, fig.height=3, fig.width=5, fig.alt="stat_early_conservation_test function output", dev.args = list(bg = 'transparent'), fig.align='center'----
modules <- list(early = 1:3, mid = 4:5, late = 6:7)
myTAI::stat_early_conservation_test(
  example_phyex_set_old, plot_result = TRUE,
  modules = modules)

## ----message = FALSE, eval = FALSE--------------------------------------------
# modules <- list(early = 1:3, mid = 4:5, late = 6:7)
# myTAI::stat_late_conservation_test(
#   example_phyex_set_old, plot_result = TRUE,
#   modules = modules)

## ----message = FALSE, warning = FALSE, echo = FALSE, fig.height=3, fig.width=5, fig.alt="stat_late_conservation_test function output", dev.args = list(bg = 'transparent'), fig.align='center'----
modules <- list(early = 1:3, mid = 4:5, late = 6:7)
myTAI::stat_late_conservation_test(
  example_phyex_set_old, plot_result = TRUE,
  modules = modules)

## ----message = FALSE, eval = FALSE--------------------------------------------
# modules <- list(contrast1 = 1:4, contrast2 = 5:7)
# myTAI::stat_pairwise_test(
#   example_phyex_set_old, plot_result = TRUE,
#   modules = modules,
#   alternative = "greater")

## ----message = FALSE, warning = FALSE, echo = FALSE, fig.height=3, fig.width=5, fig.alt="stat_pairwise_test function output", dev.args = list(bg = 'transparent'), fig.align='center'----
modules <- list(contrast1 = 1:4, contrast2 = 5:7)
myTAI::stat_pairwise_test(
  example_phyex_set_old, plot_result = TRUE,
  modules = modules,
  alternative = "greater")

