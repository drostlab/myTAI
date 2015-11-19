context("Test: taxid() ")

test_that("taxid interface works properly... ", {
        
        skip_on_cran()
        skip_on_travis()
#         tax.category <- taxid(db.path = getwd(), download = TRUE)
#         tax.category.Archea <- taxid(db.path = getwd(), filter = "Archea")
#         tax.category.Bacteria <- taxid(db.path = getwd(), filter = "Bacteria")
#         tax.category.Eukaryota <- taxid(db.path = getwd(), filter = "Eukaryota")
#         tax.category.Viruses <- taxid(db.path = getwd(), filter = "Viruses")
#         tax.category.Unclassified <- taxid(db.path = getwd(), filter = "Unclassified")
#         
#         # test general interface
#         expect_equal(tax.category[1 , "tax_id"], 7)
#         # test proper filter 'Archea'
#         expect_equal(tax.category.Archea[1 , "tax_id"], 2161)
#         # test proper filter 'Bacteria'
#         expect_equal(tax.category.Bacteria[1 , "tax_id"], 7)
#         # test proper filter 'Eukaryota'
#         expect_equal(tax.category.Eukaryota[1 , "tax_id"], 2708)
#         # test proper filter 'Viruses'
#         expect_equal(tax.category.Viruses[1 , "tax_id"], 10243)
#         # test proper filter 'Unclassified'
#         expect_equal(tax.category.Unclassified[1 , "tax_id"], 32644)
})