# Create S7 Options Property

Helper function to create an S7 property with a limited set of valid
options.

## Usage

``` r
new_options_property(
  class = S7::class_any,
  options,
  default = options[[1]],
  name = NULL
)
```

## Arguments

- class:

  The S7 class for the property (default: S7::class_any)

- options:

  Character vector of valid options for the property

- default:

  Default value for the property (default: first option)

- name:

  Optional name for the property

## Value

An S7 property object with validation for allowed options

## Details

\# Create a property that only accepts specific values \# color_prop \<-
new_options_property(options = c("red", "blue", "green"))
