#' @title Plot evolutionary signatures across transcriptomes and RNA-seq transformations
#' @description \emph{PlotSignatureTransformed} aims to statistically evaluate the
#' stability of \code{\link{ReductiveHourglassTest}}, \code{\link{FlatLineTest}}, 
#' \code{\link{ReverseHourglassTest}}, \code{\link{EarlyConservationTest}}, or 
#' \code{\link{LateConservationTest}} 
#' (all based on \code{\link{TAI}} or \code{\link{TDI}} computations) against different
#' data transformations AND plot the resulting TAI profiles using \code{\link{PlotSignature}}.
#' The corresponding p-value quantifies the probability that a given TAI or TDI pattern (or any phylotranscriptomics pattern) 
#' does not support the chosen test. A p-value < 0.05 indicates that the corresponding phylotranscriptomics pattern does
#' indeed support the chosen test.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param measure type of transcriptome index that shall be computed. E.g. 
#' \itemize{
#' \item \code{measure = "TAI"} (Transcriptome Age Index)
#' \item \code{measure = "TDI"} (Transcriptome Divergence Index)
#' \item \code{measure = "TPI"} (Transcriptome Polymorphism Index)
#' }
#' @param TestStatistic a string defining the type of test statistics to be used to quantify the statistical significance the present phylotranscriptomics pattern.
#' Possible values can be: 
#' \itemize{
#' \item \code{TestStatistic} = \code{"FlatLineTest"} : Statistical test for the deviation from a flat line
#' \item \code{TestStatistic} = \code{"ReductiveHourglassTest"} : Statistical test for the existence of a hourglass shape (high-low-high pattern)
#' \item \code{TestStatistic} = \code{"ReverseHourglassTest"} : Statistical test for the existence of a reverse hourglass pattern (low-high-low pattern)
#' \item \code{TestStatistic} = \code{"EarlyConservationTest"} : Statistical test for the existence of a early conservation pattern (low-high-high pattern)
#' \item \code{TestStatistic} = \code{"LateConservationTest"} : Statistical test for the existence of a late conservation pattern (high-high-low pattern)
#' }
#' @param transforms a character vector of any valid function that transforms gene expression levels.
#' @param modules a list storing three elements for the \code{\link{ReductiveHourglassTest}}, \code{\link{EarlyConservationTest}}, \code{\link{LateConservationTest}}, 
#' or \code{\link{ReverseHourglassTest}}: early, mid, and late. 
#' Each element expects a numeric vector specifying the developmental stages 
#' or experiments that correspond to each module. For example:
#' \itemize{
#' \item \code{module} = \code{list(early = 1:2, mid = 3:5, late = 6:7)} divides a dataset storing seven developmental stages into 3 modules.
#' } 
#' @param permutations a numeric value specifying the number of permutations to be performed for the \code{\link{FlatLineTest}}, \code{\link{EarlyConservationTest}}, \code{\link{LateConservationTest}}, \code{\link{ReductiveHourglassTest}} or \code{\link{ReverseHourglassTest}}.
#' @param pseudocount any valid number to add to the expression matrix prior to transformation.
#' @param p.value a boolean value specifying whether the p-value of the test statistic shall be printed within the plot area.
#' @param shaded.area a boolean value specifying whether a shaded area shall 
#' be drawn for the developmental stages defined to be the presumptive phylotypic period.
#' @param xlab label of x-axis.
#' @param ylab label of y-axis.
#' @param main figure title.
#' @param lwd line width.
#' @param alpha transparency of the shaded area (between [0,1]). Default is \code{alpha = 0.1}.
#' @param y.ticks number of ticks on the y-axis. Default is \code{ticks = 3}.
#' @details 
#' Visualisation and assessment for the stability of data transforms on the permutation test of choice.
#' For details, please consult the main function \code{\link{PlotSignature}},
#' as well as \code{\link{tf}}, \code{\link{ReductiveHourglassTest}}, 
#' \code{\link{FlatLineTest}}, \code{\link{ReverseHourglassTest}},
#' \code{\link{LateConservationTest}} or \code{\link{EarlyConservationTest}}.
#' 
#' In most cases, users can replace \code{\link{PlotSignature}} simply with 
#' \code{PlotSignatureTransformed} to obtain the multi-panel plot with different
#' transformations to visualise the stability of the pattern observed with \code{\link{PlotSignature}}.
#' 
#' @return a ggplot object containing the following information 
#' for each visualised transcriptome indices:
#' \code{p.value} : the p-value quantifying the statistical significance (depending on the chosen test) of the given phylotranscriptomics pattern.
#' @references 
#' 
#' Lotharukpong JS et al. (2023) bioRxiv
#'
#'   
#' @author Jaruwatana Sodai Lotharukpong
#' @seealso \code{\link{PlotSignature}},\code{\link{tfStability}}, 
#' \code{\link{FlatLineTest}}, \code{\link{ReverseHourglassTest}}, 
#' \code{\link{EarlyConservationTest}}, \code{\link{ReductiveHourglassTest}}, 
#' \code{\link{LateConservationTest}}
#' @examples
#' data(PhyloExpressionSetExample)
#' 
#' # Flat line test
#' PlotSignatureTransformed(ExpressionSet = PhyloExpressionSetExample,
#'                      TestStatistic = "FlatLineTest",
#'                      transforms = c("log2", "sqrt", "none"))
#' 
#' # Reductive hourglass test
#' PlotSignatureTransformed(ExpressionSet = PhyloExpressionSetExample,
#'                      TestStatistic = "ReductiveHourglassTest",
#'                      transforms = c("log2", "sqrt", "none"),
#'                      modules = list(early = 1:2, mid = 3:5, late = 6:7))
#' 
#' library(DESeq2)
#' PlotSignatureTransformed(ExpressionSet = PhyloExpressionSetExample,
#'                      TestStatistic = "ReductiveHourglassTest",
#'                      transforms = c("log2", "sqrt", "none", "vst"),
#'                      modules = list(early = 1:2, mid = 3:5, late = 6:7))
#'                      
#' @export
PlotSignatureTransformed <- 
  function(ExpressionSet,
           measure            = "TAI",
           TestStatistic      = "FlatLineTest",
           transforms         = c("none", "sqrt", "log2", "rank", "squared"),
           modules            = NULL,
           permutations       = 1000,
           pseudocount        = 1,
           p.value            = FALSE,
           shaded.area        = FALSE,
           xlab               = "Ontogeny",
           ylab               = "Transcriptome Index",
           main               = "",
           lwd                = 4,
           alpha              = 0.1,
           y.ticks            = 3){
    
    # Check input data and output error messages when 
    myTAI::is.ExpressionSet(ExpressionSet)
    
    if(!TestStatistic %in% c("FlatLineTest", "ReductiveHourglassTest", "ReverseHourglassTest", "EarlyConservationTest", "LateConservationTest"))
      stop("Please select the available test: 'FlatLineTest', 'ReductiveHourglassTest', 'ReverseHourglassTest', 'EarlyConservationTest' or 'LateConservationTest' using the argument test = 'FlatLineTest'", call. = FALSE)
    
    if(TestStatistic %in% c("ReductiveHourglassTest", "ReverseHourglassTest", "EarlyConservationTest", "LateConservationTest") & is.null(modules))
      stop("Please specify the three modules: early, mid, and late using the argument 'module = list(early = ..., mid = ..., late = ...)'.", call. = FALSE)
    
    # message(paste("Proceeding with the", TestStatistic))
    
    # However, this step repeats the permutation test. 
    # In future releases, I recommend removing it and redirecting the p-value
    # from the PlotSignature via
    # ggplot2::ggplot_build(p[[i]])$data[[3]][["label"]]
    # but the existing annotation cannot be removed post facto.
    p_value <- 
      tfStability(
        ExpressionSet,
        TestStatistic      = TestStatistic,
        transforms         = transforms,
        modules            = modules,
        permutations       = permutations,
        pseudocount        = pseudocount)
    p_value <- 
      scales::pvalue(p_value, accuracy= 0.0001, add_p = TRUE)
    
    # Set the ggplot-space to include x number of plots (from the number of transformations chosen)
    p <- list()
    
    # output transforms in a vector
    for (i in transforms) {
      if(i == "none")
        tfExpressionSet <- ExpressionSet
      else if(i %in% c('log2', 'log', 'log10'))
        tfExpressionSet <- tf(ExpressionSet, FUN = i, pseudocount = pseudocount)
      else if(i == "squared")
        tfExpressionSet <- tf(ExpressionSet, FUN = function(x) x*x)
      else if(i %in% c('vst', 'rlog'))
        tfExpressionSet <- tf(ExpressionSet, FUN = i, integerise = TRUE)
      else if(i == "rank")
        tfExpressionSet <- tf(ExpressionSet, FUN = function(x) apply(x, 2, base::rank))
      else
        tfExpressionSet <- tf(ExpressionSet, FUN = i)
      
      tfExpressionSet <- stats::na.omit(tfExpressionSet)
      
      # Here we make the plots using PlotSignature
      
      p[[i]] <- 
        PlotSignature(
          tfExpressionSet,
          measure = measure,
          TestStatistic = TestStatistic,
          modules = modules,
          permutations = permutations,
          p.value = p.value,
          shaded.area  = shaded.area,
          xlab = NULL,
          ylab = NULL,
          main = i,
          lwd = lwd,
          alpha = alpha,
          y.ticks = y.ticks) + 
        ggplot2::scale_x_discrete(labels = NULL) +
        ggplot2::theme(
          axis.text.y = ggplot2::element_text(size = 15),
          plot.subtitle = ggplot2::element_text(size=10)
        )
      p[[i]] <- p[[i]] +
        ggplot2::labs(
          subtitle = p_value[[i]]
        )
      # I tried to use the precomputed p-values from PlotSignature() and
      # replace the subtitle with it. However, I cannot remove the annotations.
      # I have opened an issue in the github.
      # https://github.com/drostlab/myTAI/issues/16
      # p[[i]] <- p[[i]] +
      #   ggplot2::labs(
      #     subtitle = ggplot2::ggplot_build(p[[i]])$data[[3]][["label"]]
      #   )
    } 
    
    p_out <- 
      base::do.call(
        ggpubr::ggarrange,
        p
      )
    
    p_out <-
      ggpubr::annotate_figure(
        p_out,
        top = ggpubr::text_grob(
          TestStatistic, 
          face = "bold", 
          size = 20),
        left = ggpubr::text_grob(
          ylab, 
          rot = 90,
          face = "bold", 
          size = 20),
        bottom = ggpubr::text_grob(
          xlab, 
          face = "bold", 
          size = 20)
      )
    
    return(p_out)
  }