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

new_required_property <- function(class = S7::class_any,
                                  name,
                                  validator = NULL) {
    msg <- sprintf("@%s is required", name)
    return(S7::new_property(
        class = class,
        default = substitute(stop(msg), list(msg = msg)),
        validator = validator,
        name = name
    ))
}
