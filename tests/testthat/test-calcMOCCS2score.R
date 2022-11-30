test_that("test of calcMOCCS2score", {
    load("trueResult.rda")
    fastaPath <- system.file("count.fasta", package = "moccsR")
    testResult <- calcMOCCS2score(fastaPath, 6)
    expect_identical(testResult, trueResult)
})
