# 
# 
# phylo_set <- S7::new_class("phylo_set",
#                            properties=list(
#                                data = S7::class_data.frame,
#                                rna_transf="sqrt",
#                                type="tai"
#                            ),
#                            validator = function(self) {
#                                # check data is of decent quality and shoot warnings if not!
#                            }
#                        )
# 
# 
# # add null hypothesis sample once computed