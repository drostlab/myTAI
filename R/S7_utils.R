
# adapted from https://github.com/RConsortium/S7/issues/418
#' @import memoise cachem
new_cached_property <- function(class = S7::class_any, 
                                getter,
                                setter = NULL,
                                validator = NULL,
                                default = NULL,
                                name = NULL) {
    mc <<- cachem::cache_mem()
    
    cached_getter <- function(self) {
        f <- function(self, getter_f) {
            getter_f(self)
        }
        mf <- memoise::memoise(f, cache = mc)
        
        out <- mf(self, getter)
        return(out)
    }
    
    return(S7::new_property(class = class,
                            getter = cached_getter,
                            setter = setter,
                            validator = validator,
                            default = default,
                            name = name))
}

new_options_property <- function(class = S7::class_any,
                                 options,
                                 default = options[[1]],
                                 name = NULL) {
    return(S7::new_property(
        class = class,
        default = default,
        validator = \(value) if (! value %in% options) paste("must be one of", paste(options, collapse=", ")),
        name = name
    ))
}

new_required_property <- function(class = S7::class_any,
                                  name,
                                  validator = NULL
                                  ) {
    msg <- sprintf("@%s is required", name)
    return(S7::new_property(
        class = class,
        default = substitute(stop(msg), list(msg=msg)),
        validator = validator,
        name = name
    ))
}


