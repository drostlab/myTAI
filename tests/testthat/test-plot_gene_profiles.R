gene_options <- list(
    NULL,
    example_phyex_set@gene_ids[1],
    example_phyex_set@gene_ids[1:10],
    example_phyex_set@gene_ids[1:100]
)

params <- expand.grid(
    transformation = c("log", "std_log", "none"),
    colour_by = c("manual", "stage", "strata"),
    show_bg = c(TRUE, FALSE),
    show_set_mean = c(TRUE, FALSE),
    show_reps = c(TRUE, FALSE),
    show_labels = c(TRUE, FALSE),
    show_legend = c(TRUE, FALSE),
    stringsAsFactors = FALSE
)

params <- params[rep(seq_len(nrow(params)), each = length(gene_options)), ]
params$genes <- rep(gene_options, times = nrow(params) / length(gene_options))

test_that("plot_gene_profiles handles all options", {
    for (i in seq_len(nrow(params))) {
        p <- plot_gene_profiles(
            example_phyex_set,
            genes = params$genes[[i]],
            transformation = params$transformation[i],
            colour_by = params$colour_by[i],
            show_bg = params$show_bg[i],
            show_set_mean = params$show_set_mean[i],
            show_reps = params$show_reps[i],
            show_labels = params$show_labels[i],
            show_legend = params$show_legend[i]
        )
        label <- paste0("plot_", params$transformation[i], "_", params$colour_by[i], 
                        "_bg", params$show_bg[i], "_mean", params$show_set_mean[i],
                        "_reps", params$show_reps[i], "_labels", params$show_labels[i], 
                        "_legend", params$show_legend[i], "_genes", length(params$genes[[i]]))
        vdiffr::expect_doppelganger(label, p)
        expect_s3_class(p, "ggplot")
    }
})