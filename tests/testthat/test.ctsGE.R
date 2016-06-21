
data_dir <- system.file("extdata", package = "ctsGE")
files <- dir(path=data_dir,pattern = "\\.xls$")
rts <- readTSGE(files, path = data_dir, skip = 10625 )
prts <- PreparingTheIndexes(rts)
tsCI <- ClustIndexes(prts)

expect_is(rts,"ctsGEList")
expect_is(prts,"ctsGEList")
expect_is(tsCI,"ctsGEList")
expect_equal(length(rts),4)
expect_equal(length(prts),7)
expect_equal(length(tsCI),9)
expect_equal(sum(is.na(tsCI$optimalK[,2])),0)
expect_error(PlotIndexesClust(prts,idx="lkj"))
expect_error(PlotIndexesClust(prts,idx=000000))
expect_error(PlotIndexesClust(rts,idx="000000"),
             "Please run PreparingTheIndexes first")
