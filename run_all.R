curdir <- getwd()
Sys.setenv(AHMnimble_comparisons = "test_comparison_pages")
chapters_to_run <- c(
#    "Chapter_5"
#   ,
#    "Chapter_6"
#,"Chapter_7"
#,"Chapter_8"
    "Chapter_9"
    ,"Chapter_10"
)

run_chapter <- function(chapter) {
    setwd(chapter)
    on.exit( setwd(curdir) )
    system2("caffeinate", c("Rscript", "run_all.R"))
}

for(chapter in chapters_to_run) {
    run_chapter(chapter)
}
