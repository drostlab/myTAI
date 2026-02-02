# Create S7 Required Property

Helper function to create an S7 property that is required (throws error
if not provided).

## Usage

``` r
new_required_property(
  class = S7::class_any,
  name,
  validator = NULL,
  getter = NULL,
  setter = NULL
)
```

## Arguments

- class:

  The S7 class for the property (default: S7::class_any)

- name:

  Name of the property (used in error messages)

- validator:

  Optional validation function for the property

- getter:

  Optional getter function for the property

- setter:

  Optional setter function for the property

## Value

An S7 property object that is required

## Details

\# Create a required property \# required_prop \<-
new_required_property(class = class_character, name = "gene_ids")
