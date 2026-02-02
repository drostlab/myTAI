# Retrieve taxonomy categories from NCBI Taxonomy

This function retrieves category information from NCBI Taxonomy and is
able to filter kingdom specific taxids.

## Usage

``` r
taxid(db.path, download = FALSE, update = FALSE, filter = NULL)
```

## Arguments

- db.path:

  path to download and store the NCBI Taxonomy `categories.dmp` file.
  Default is the [`tempdir()`](https://rdrr.io/r/base/tempfile.html)
  directory.

- download:

  a logical value specifying whether or not the `categories.dmp` shall
  be downloaded (`download = TRUE`) or whether a local version already
  exists on the users machine (`download = TRUE` - in this case please
  specify the `db.path` argument to target the local `categories.dmp`
  file).

- update:

  should the local file be updated? Please specify the `db.path`
  argument to target the local `categories.dmp` file.

- filter:

  a character string specifying the kingdom of life for which taxids
  shall be returned. Options are `"Archea"`, `"Bacteria"`,
  `"Eukaryota"`, `"Viruses"`, `"Unclassified"`.

## Value

A tibble object containing taxonomy category information

## Author

Hajk-Georg Drost

## Examples

``` r
if (FALSE) {
# download categories.dmp file to current working directory 
# and filter for 'Archea' taxids
Archea.taxids <- taxid(db.path = getwd(), filter = "Archea", download = TRUE)

# Once the NCBI Taxonomy 'categories.dmp' file is downloaded to your machine ('download = TRUE')
# the 'taxid()' function can be proceed on the local 'categories.dmp' file
# e.g. filter for Virus taxids
Virus.taxids <- taxid(db.path = getwd(), filter = "Viruses")
}
```
