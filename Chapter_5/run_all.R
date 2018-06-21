allFiles <- list.files(pattern = glob2rx("Section_5*nimble.R"))
comparisonDir <- Sys.getenv("AHMnimble_comparisons")
if(identical(comparisonDir, ""))
    compariseonDir <- "test_comparison_pages"
for(f in allFiles) {
    args <- paste0("-e \"print(getwd()); library(methods); ",
                  "outputDirectory  <- '", comparisonDir, "';",
                  "source('",
                  f,
                  "')\"")
    print(args)
    system2("Rscript", args)
}
q('no')
