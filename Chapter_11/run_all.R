filter <- commandArgs(trailingOnly = TRUE)
if(length(filter) > 0) {
    if(length(filter) > 1)
        warning(paste0("Ignoring all elements of file filter except ", filter[1]))
    allFiles <- list.files(pattern = glob2rx(filter[1]))
} else 
    allFiles <- list.files(pattern = glob2rx("Section_11*nimble.R"))
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
