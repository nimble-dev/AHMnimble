allFiles <- list.files(pattern = glob2rx("Section_7*nimble.R"))
comparisonDir <- Sys.getenv("AHMnimble_comparisons")
if(identical(comparisonDir, ""))
    comparisonDir <- "test_comparison_pages"
for(f in allFiles) {
    args <- paste0("-e \"print(getwd()); library(methods); ",
                  "outputDirectory  <- '", comparisonDir, "';",
                  "source('",
                  f,
                  "')\"")
    print(args)
    if(.Platform$OS.type == "unix")
      system2("caffeinate", c("Rscript", args))
    else
    system2("Rscript", args)
}
q('no')
