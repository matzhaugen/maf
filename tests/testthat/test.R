mafObject = maf(treeringTimeseries)
expect_equal(21, dim(mafObject$mafs)[2])