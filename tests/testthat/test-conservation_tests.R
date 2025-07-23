test_that("Flatline test works", {
    # Basic flatline test
    result <- stat_flatline_test(example_phyex_set, plot_result = FALSE)
    expect_s7_class(result, myTAI::ConservationTestResult)
    expect_true(is.numeric(result@p_value))
    expect_true(result@p_value >= 0 && result@p_value <= 1)
    expect_true(is.numeric(result@test_stat))
    expect_equal(result@method_name, "Flat Line Test")
})

test_that("Reductive hourglass test works", {
    # Test with early-mid-late modules
    modules <- list(
        early = c(1, 2),
        mid = c(3, 4), 
        late = c(5, 6)
    )
    
    result <- stat_reductive_hourglass_test(example_phyex_set, modules, plot_result = FALSE)
    expect_s7_class(result, myTAI::ConservationTestResult)
    expect_true(is.numeric(result@p_value))
    expect_true(result@p_value >= 0 && result@p_value <= 1)
    expect_equal(result@method_name, "Reductive Hourglass Test")
    expect_equal(result@modules, modules)
})

test_that("Reverse hourglass test works", {
    # Test with early-mid-late modules
    modules <- list(
        early = c(1, 2),
        mid = c(3, 4),
        late = c(5, 6)
    )
    
    result <- stat_reverse_hourglass_test(example_phyex_set, modules, plot_result = FALSE)
    expect_s7_class(result, myTAI::ConservationTestResult)
    expect_true(is.numeric(result@p_value))
    expect_true(result@p_value >= 0 && result@p_value <= 1)
    expect_equal(result@method_name, "Reverse Hourglass Test")
})

test_that("Early conservation test works", {
    modules <- list(
        early = c(1, 2, 3),
        mid = c(4, 5),
        late = c(6)
    )
    
    result <- stat_early_conservation_test(example_phyex_set, modules, plot_result = FALSE)
    expect_s7_class(result, myTAI::ConservationTestResult)
    expect_true(is.numeric(result@p_value))
    expect_true(result@p_value >= 0 && result@p_value <= 1)
    expect_equal(result@method_name, "Early Conservation Test")
})

test_that("Late conservation test works", {
    modules <- list(
        early = c(1, 2),
        mid = c(3, 4),
        late = c(5, 6)
    )
    
    result <- stat_late_conservation_test(example_phyex_set, modules, plot_result = FALSE)
    expect_s7_class(result, myTAI::ConservationTestResult)
    expect_true(is.numeric(result@p_value))
    expect_true(result@p_value >= 0 && result@p_value <= 1)
    expect_equal(result@method_name, "Late Conservation Test")
})

test_that("Pairwise test works", {
    # Test with stage pairs
    stage_pairs <- list(
        contrast1 = c(1, 2, 3),
        contrast2 = c(4, 5, 6)
    )
    
    result <- stat_pairwise_test(example_phyex_set, stage_pairs, plot_result = FALSE)
    expect_s7_class(result, myTAI::ConservationTestResult)
    expect_true(is.numeric(result@p_value))
    expect_true(result@p_value >= 0 && result@p_value <= 1)
    expect_equal(result@method_name, "Pairwise Test")
})



test_that("Test result helper functions work", {
    result <- stat_flatline_test(example_phyex_set, plot_result = FALSE)
    
    # Test confidence intervals (internal function)
    ci <- myTAI:::conf_int(result)
    expect_true(is.numeric(ci))
    expect_equal(length(ci), 2)
    expect_true(ci[1] <= ci[2])
    
    # Test goodness of fit (internal function)
    gof <- myTAI:::goodness_of_fit(result)
    expect_true(is.list(gof) || inherits(gof, "htest"))
    
    # Test exp_p (internal function) - returns call object for ggplot2
    exp_p_val <- myTAI:::exp_p(result@p_value)
    expect_true(is.call(exp_p_val))
    
    # Test that it can be converted to character for ggplot2
    expect_true(is.character(as.character(exp_p_val)))
})

test_that("tf_stability works", {
    # Test with a simple transformation set
    simple_transforms <- list(
        identity = function(x) x,
        log1p = function(x) log1p(x)
    )
    
    stability_result <- tf_stability(
        example_phyex_set, 
        conservation_test = function(x, ...) flatline_test(x, ...),
        transforms = simple_transforms
    )
    
    expect_true(is.numeric(stability_result))
    expect_equal(length(stability_result), length(simple_transforms))
    expect_equal(names(stability_result), names(simple_transforms))
    expect_true(all(stability_result >= 0 & stability_result <= 1))
})
