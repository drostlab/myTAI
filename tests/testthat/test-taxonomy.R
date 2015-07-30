context("Test: taxonomy() ")


# test_that("db = eol interface works properly...",{
#         
#         test.tax <- taxonomy("Arabidopsis thaliana",db = "eol")
#         expect_equal(test.tax[nrow(test.tax), 1],"Arabidopsis thaliana")
#         
#         
#         
#         test.tax2 <- taxonomy("Arabidopsis",db = "eol")
#         expect_equal(test.tax2[nrow(test.tax2), 1],"Arabidopsis")
#         
#         
#         
#         test.tax3 <- taxonomy("Homo sapiens",db = "eol")
#         expect_equal(test.tax3[nrow(test.tax3), 1],"Homo sapiens")
#         
#         
#         
#         test.tax4 <- taxonomy("Mus musculus",db = "eol")
#         expect_equal(test.tax4[nrow(test.tax4), 1],"Mus musculus")
#         
#         
#         
#         test.tax5 <- taxonomy("Danio rerio",db = "eol")
#         expect_equal(test.tax5[nrow(test.tax5), 1],"Danio rerio")
#         
#         
#         expect_equal(as.numeric(taxonomy("Arabidopsis thaliana",db = "eol", output = "taxid")),3702)
#         
# })

# test_that("db = tnrs interface works properly...",{
#         
#         test.tax <- taxonomy("Arabidopsis thaliana",db = "tnrs")
#         expect_equal(test.tax[nrow(test.tax), 1],"Arabidopsis thaliana")
#         
#         
#         
#         test.tax2 <- taxonomy("Arabidopsis",db = "tnrs")
#         expect_equal(test.tax2[nrow(test.tax2), 1],"Arabidopsis")
#         
#         
#         
#         test.tax3 <- taxonomy("Homo sapiens",db = "tnrs")
#         expect_equal(test.tax3[nrow(test.tax3), 1],"Homo sapiens")
#         
#         
#         
#         test.tax4 <- taxonomy("Mus musculus",db = "tnrs")
#         expect_equal(test.tax4[nrow(test.tax4), 1],"Mus musculus")
#         
#         
#         
#         test.tax5 <- taxonomy("Danio rerio",db = "tnrs")
#         expect_equal(test.tax5[nrow(test.tax5), 1],"Danio rerio")
#         
#         expect_equal(as.numeric(taxonomy("Arabidopsis thaliana",db = "tnrs", output = "taxid")),3702)
# })

# test_that("db = itis interface works properly...",{
#         
#         skip_on_cran()
#         
#         test.tax <- taxonomy("Arabidopsis thaliana",db = "itis")
#         expect_equal(test.tax[nrow(test.tax), 1],"Arabidopsis thaliana")
#         
#         
#         
#         test.tax2 <- taxonomy("Arabidopsis",db = "itis")
#         expect_equal(test.tax2[nrow(test.tax2), 1],"Arabidopsis")
#         
#         
#         
#         test.tax3 <- taxonomy("Homo sapiens",db = "itis")
#         expect_equal(test.tax3[nrow(test.tax3), 1],"Homo sapiens")
#         
#         
#         
#         test.tax4 <- taxonomy("Mus musculus",db = "itis")
#         expect_equal(test.tax4[nrow(test.tax4), 1],"Mus musculus")
#         
#         
#         
#         test.tax5 <- taxonomy("Danio rerio",db = "itis")
#         expect_equal(test.tax5[nrow(test.tax5), 1],"Danio rerio")
#         
#         expect_equal(as.numeric(taxonomy("Arabidopsis thaliana",db = "itis", output = "taxid")),23041)
# })


# test_that("db = phylomatic interface works properly...",{
#         
#         test.tax <- taxonomy("Arabidopsis thaliana",db = "phylomatic")
#         expect_equal(test.tax[nrow(test.tax), 1],"Arabidopsis thaliana")
#         
#         
#         
#         test.tax2 <- taxonomy("Arabidopsis",db = "phylomatic")
#         expect_equal(test.tax2[nrow(test.tax2), 1],"Arabidopsis")
#         
#         
#         
#         test.tax3 <- taxonomy("Homo sapiens",db = "phylomatic")
#         expect_equal(test.tax3[nrow(test.tax3), 1],"Homo sapiens")
#         
#         
#         
#         test.tax4 <- taxonomy("Mus musculus",db = "phylomatic")
#         expect_equal(test.tax4[nrow(test.tax4), 1],"Mus musculus")
#         
#         
#         
#         test.tax5 <- taxonomy("Danio rerio",db = "phylomatic")
#         expect_equal(test.tax5[nrow(test.tax5), 1],"Danio rerio")
#         
#         
# })

# test_that("db = ubio interface works properly...",{
#         
#         test.tax <- taxonomy("Arabidopsis thaliana",db = "ubio")
#         expect_equal(test.tax[nrow(test.tax), 1],"Arabidopsis thaliana")
#         
#         
#         
#         test.tax2 <- taxonomy("Arabidopsis",db = "ubio")
#         expect_equal(test.tax2[nrow(test.tax2), 1],"Arabidopsis")
#         
#         
#         
#         test.tax3 <- taxonomy("Homo sapiens",db = "ubio")
#         expect_equal(test.tax3[nrow(test.tax3), 1],"Homo sapiens")
#         
#         
#         
#         test.tax4 <- taxonomy("Mus musculus",db = "ubio")
#         expect_equal(test.tax4[nrow(test.tax4), 1],"Mus musculus")
#         
#         
#         
#         test.tax5 <- taxonomy("Danio rerio",db = "ubio")
#         expect_equal(test.tax5[nrow(test.tax5), 1],"Danio rerio")
#         
#         
# })
# 
# test_that("db = gnr interface works properly...",{
#         
#         test.tax <- taxonomy("Arabidopsis thaliana",db = "gnr")
#         expect_equal(test.tax[nrow(test.tax), 1],"Arabidopsis thaliana")
#         
#         
#         
#         test.tax2 <- taxonomy("Arabidopsis",db = "gnr")
#         expect_equal(test.tax2[nrow(test.tax2), 1],"Arabidopsis")
#         
#         
#         
#         test.tax3 <- taxonomy("Homo sapiens",db = "gnr")
#         expect_equal(test.tax3[nrow(test.tax3), 1],"Homo sapiens")
#         
#         
#         
#         test.tax4 <- taxonomy("Mus musculus",db = "gnr")
#         expect_equal(test.tax4[nrow(test.tax4), 1],"Mus musculus")
#         
#         
#         
#         test.tax5 <- taxonomy("Danio rerio",db = "gnr")
#         expect_equal(test.tax5[nrow(test.tax5), 1],"Danio rerio")
#         
#         
# })
# 
# test_that("db = gni interface works properly...",{
#         
#         test.tax <- taxonomy("Arabidopsis thaliana",db = "gni")
#         expect_equal(test.tax[nrow(test.tax), 1],"Arabidopsis thaliana")
#         
#         
#         
#         test.tax2 <- taxonomy("Arabidopsis",db = "gni")
#         expect_equal(test.tax2[nrow(test.tax2), 1],"Arabidopsis")
#         
#         
#         
#         test.tax3 <- taxonomy("Homo sapiens",db = "gni")
#         expect_equal(test.tax3[nrow(test.tax3), 1],"Homo sapiens")
#         
#         
#         
#         test.tax4 <- taxonomy("Mus musculus",db = "gni")
#         expect_equal(test.tax4[nrow(test.tax4), 1],"Mus musculus")
#         
#         
#         
#         test.tax5 <- taxonomy("Danio rerio",db = "gni")
#         expect_equal(test.tax5[nrow(test.tax5), 1],"Danio rerio")
#         
#         
# })
# 
# test_that("db = iucn interface works properly...",{
#         
#         test.tax <- taxonomy("Arabidopsis thaliana",db = "iucn")
#         expect_equal(test.tax[nrow(test.tax), 1],"Arabidopsis thaliana")
#         
#         
#         
#         test.tax2 <- taxonomy("Arabidopsis",db = "iucn")
#         expect_equal(test.tax2[nrow(test.tax2), 1],"Arabidopsis")
#         
#         
#         
#         test.tax3 <- taxonomy("Homo sapiens",db = "iucn")
#         expect_equal(test.tax3[nrow(test.tax3), 1],"Homo sapiens")
#         
#         
#         
#         test.tax4 <- taxonomy("Mus musculus",db = "iucn")
#         expect_equal(test.tax4[nrow(test.tax4), 1],"Mus musculus")
#         
#         
#         
#         test.tax5 <- taxonomy("Danio rerio",db = "iucn")
#         expect_equal(test.tax5[nrow(test.tax5), 1],"Danio rerio")
#         
#         
# })

# test_that("db = tp interface works properly...",{
#         
#         test.tax <- taxonomy("Arabidopsis thaliana",db = "tp")
#         expect_equal(test.tax[nrow(test.tax), 1],"Arabidopsis thaliana")
#         
#         
#         
#         test.tax2 <- taxonomy("Arabidopsis",db = "tp")
#         expect_equal(test.tax2[nrow(test.tax2), 1],"Arabidopsis")
#         
#         
#         
#         test.tax3 <- taxonomy("Homo sapiens",db = "tp")
#         expect_equal(test.tax3[nrow(test.tax3), 1],"Homo sapiens")
#         
#         
#         
#         test.tax4 <- taxonomy("Mus musculus",db = "tp")
#         expect_equal(test.tax4[nrow(test.tax4), 1],"Mus musculus")
#         
#         
#         
#         test.tax5 <- taxonomy("Danio rerio",db = "tp")
#         expect_equal(test.tax5[nrow(test.tax5), 1],"Danio rerio")
#         
#         
# })

# test_that("db = plantminer interface works properly...",{
#         
#         test.tax <- taxonomy("Arabidopsis thaliana",db = "plantminer")
#         expect_equal(test.tax[nrow(test.tax), 1],"Arabidopsis thaliana")
#         
#         
#         
#         test.tax2 <- taxonomy("Arabidopsis",db = "plantminer")
#         expect_equal(test.tax2[nrow(test.tax2), 1],"Arabidopsis")
#         
#         
#         
# })

# test_that("db = col interface works properly...",{
#         
#         test.tax <- taxonomy("Arabidopsis thaliana",db = "col")
#         expect_equal(test.tax[nrow(test.tax), 1],"Arabidopsis thaliana")
#         
#         
#         
#         test.tax2 <- taxonomy("Arabidopsis",db = "col")
#         expect_equal(test.tax2[nrow(test.tax2), 1],"Arabidopsis")
#         
#         
#         
#         test.tax3 <- taxonomy("Homo sapiens",db = "col")
#         expect_equal(test.tax3[nrow(test.tax3), 1],"Homo sapiens")
#         
#         
#         
#         test.tax4 <- taxonomy("Mus musculus",db = "col")
#         expect_equal(test.tax4[nrow(test.tax4), 1],"Mus musculus")
#         
#         
#         
#         test.tax5 <- taxonomy("Danio rerio",db = "col")
#         expect_equal(test.tax5[nrow(test.tax5), 1],"Danio rerio")
#         
#         
# })

# test_that("db = ncbi interface works properly...",{
#         
#         skip_on_cran()
#         
#         test.tax <- taxonomy("Arabidopsis thaliana",db = "ncbi")
#         expect_equal(test.tax[nrow(test.tax), 1],"Arabidopsis thaliana")
#         
#         
#         
#         test.tax2 <- taxonomy("Arabidopsis",db = "ncbi")
#         expect_equal(test.tax2[nrow(test.tax2), 1],"Arabidopsis")
#         
#         
#         
#         test.tax3 <- taxonomy("Homo sapiens",db = "ncbi")
#         expect_equal(test.tax3[nrow(test.tax3), 1],"Homo sapiens")
#         
#         
#         
#         test.tax4 <- taxonomy("Mus musculus",db = "ncbi")
#         expect_equal(test.tax4[nrow(test.tax4), 1],"Mus musculus")
#         
#         
#         
#         test.tax5 <- taxonomy("Danio rerio",db = "ncbi")
#         expect_equal(test.tax5[nrow(test.tax5), 1],"Danio rerio")
#         
#         
# })

# test_that("db = vascan interface works properly...",{
#         
#         test.tax <- taxonomy("Arabidopsis thaliana",db = "vascan")
#         expect_equal(test.tax[nrow(test.tax), 1],"Arabidopsis thaliana")
#         
#         
#         
#         test.tax2 <- taxonomy("Arabidopsis",db = "vascan")
#         expect_equal(test.tax2[nrow(test.tax2), 1],"Arabidopsis")
#         
#         
#         
#         test.tax3 <- taxonomy("Homo sapiens",db = "vascan")
#         expect_equal(test.tax3[nrow(test.tax3), 1],"Homo sapiens")
#         
#         
#         
#         test.tax4 <- taxonomy("Mus musculus",db = "vascan")
#         expect_equal(test.tax4[nrow(test.tax4), 1],"Mus musculus")
#         
#         
#         
#         test.tax5 <- taxonomy("Danio rerio",db = "vascan")
#         expect_equal(test.tax5[nrow(test.tax5), 1],"Danio rerio")
#         
#         
# })
# 
# test_that("db = ipni interface works properly...",{
#         
#         test.tax <- taxonomy("Arabidopsis thaliana",db = "ipni")
#         expect_equal(test.tax[nrow(test.tax), 1],"Arabidopsis thaliana")
#         
#         
#         
#         test.tax2 <- taxonomy("Arabidopsis",db = "ipni")
#         expect_equal(test.tax2[nrow(test.tax2), 1],"Arabidopsis")
#         
#         
#         
#         test.tax3 <- taxonomy("Homo sapiens",db = "ipni")
#         expect_equal(test.tax3[nrow(test.tax3), 1],"Homo sapiens")
#         
#         
#         
#         test.tax4 <- taxonomy("Mus musculus",db = "ipni")
#         expect_equal(test.tax4[nrow(test.tax4), 1],"Mus musculus")
#         
#         
#         
#         test.tax5 <- taxonomy("Danio rerio",db = "ipni")
#         expect_equal(test.tax5[nrow(test.tax5), 1],"Danio rerio")
#         
#         
# })
# 
# test_that("db = bold interface works properly...",{
#         
#         test.tax <- taxonomy("Arabidopsis thaliana",db = "bold")
#         expect_equal(test.tax[nrow(test.tax), 1],"Arabidopsis thaliana")
#         
#         
#         
#         test.tax2 <- taxonomy("Arabidopsis",db = "bold")
#         expect_equal(test.tax2[nrow(test.tax2), 1],"Arabidopsis")
#         
#         
#         
#         test.tax3 <- taxonomy("Homo sapiens",db = "bold")
#         expect_equal(test.tax3[nrow(test.tax3), 1],"Homo sapiens")
#         
#         
#         
#         test.tax4 <- taxonomy("Mus musculus",db = "bold")
#         expect_equal(test.tax4[nrow(test.tax4), 1],"Mus musculus")
#         
#         
#         
#         test.tax5 <- taxonomy("Danio rerio",db = "bold")
#         expect_equal(test.tax5[nrow(test.tax5), 1],"Danio rerio")
#         
#         
# })

# test_that("db = nbn interface works properly...",{
#         
#         test.tax <- taxonomy("Arabidopsis thaliana",db = "nbn")
#         expect_equal(test.tax[nrow(test.tax), 1],"Arabidopsis thaliana")
#         
#         
#         
#         test.tax2 <- taxonomy("Arabidopsis",db = "nbn")
#         expect_equal(test.tax2[nrow(test.tax2), 1],"Arabidopsis")
#         
#         
#         
#         test.tax3 <- taxonomy("Homo sapiens",db = "nbn")
#         expect_equal(test.tax3[nrow(test.tax3), 1],"Homo sapiens")
#         
#         
#         
#         test.tax4 <- taxonomy("Mus musculus",db = "nbn")
#         expect_equal(test.tax4[nrow(test.tax4), 1],"Mus musculus")
#         
#         
#         
#         test.tax5 <- taxonomy("Danio rerio",db = "nbn")
#         expect_equal(test.tax5[nrow(test.tax5), 1],"Danio rerio")
# })



