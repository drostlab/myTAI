#' @title Create S7 Options Property
#' @description Helper function to create an S7 property with a limited set of valid options.
#' 
#' @param class The S7 class for the property (default: S7::class_any)
#' @param options Character vector of valid options for the property
#' @param default Default value for the property (default: first option)
#' @param name Optional name for the property
#' 
#' @return An S7 property object with validation for allowed options
#' 
#' @details
#' # Create a property that only accepts specific values
#' # color_prop <- new_options_property(options = c("red", "blue", "green"))
#' 
#' @keywords internal
new_options_property <- function(class = S7::class_any,
                                 options,
                                 default = options[[1]],
                                 name = NULL) {
    return(S7::new_property(
        class = class,
        default = default,
        validator = \(value) if (!value %in% options) paste("must be one of", paste(options, collapse = ", ")),
        name = name
    ))
}

#' @title Create S7 Required Property
#' @description Helper function to create an S7 property that is required (throws error if not provided).
#' 
#' @param class The S7 class for the property (default: S7::class_any)
#' @param name Name of the property (used in error messages)
#' @param validator Optional validation function for the property
#' @param getter Optional getter function for the property
#' @param setter Optional setter function for the property
#' 
#' @return An S7 property object that is required
#' 
#' @details
#' # Create a required property
#' # required_prop <- new_required_property(class = class_character, name = "gene_ids")
#' 
#' @keywords internal
new_required_property <- function(class = S7::class_any,
                                  name,
                                  validator = NULL,
                                  getter = NULL,
                                  setter = NULL) {
    msg <- sprintf("@%s is required", name)
    return(S7::new_property(
        class = class,
        default = substitute(stop(msg), list(msg = msg)),
        validator = validator,
        name = name,
        getter = getter,
        setter = setter
    ))
}
