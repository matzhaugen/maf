setwd("~/cloud/MafPackageR/")
current.code = as.package("maf")
devtools::load_all(current.code)
devtools::document(current.code)