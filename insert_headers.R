## This code inserts the contents of header.txt at the beginning of every
## R file (except a few) in the Chapter directories.
## It also removes any old version of the header content, if present.
## This allows the header to be edited and re-inserted everywhere, if needed.
## This is dangerous, since a bug could damage or erase contents of code files.
## We expect this only to be used in a github repository, allowing any
## mistakes to be discarded.
chapters <- 5:11
chapterDirs <- paste0("Chapter_",chapters)
headerLines <- readLines("header.txt")
for(chapterDir in chapterDirs) {
    allRfiles <- list.files(chapterDir)
    allRfiles <- allRfiles[grepl(".R", allRfiles)]
    allRfiles <- allRfiles[!grepl("run_all.R", allRfiles)]
    for(oneRfile in allRfiles) {
        pathedOneRfile <- file.path(chapterDir, oneRfile)
        allLines <- readLines(pathedOneRfile)
        if(!(isTRUE(any(grepl("your_comparison_pages", allLines))))) {
            startOldHeader <- grep("(begin AHMnimble header)", allLines)
            endOldHeader <- grep("(end AHMnimble header)", allLines)
            if(length(startOldHeader) > 1 | length(endOldHeader) > 1)
                stop(paste0("Problem picking old header out of ", oneRfile))
            if(length(startOldHeader) == 1 & length(endOldHeader) == 1) {
                allLines <- allLines[-(startOldHeader:endOldHeader)]
            }
            writeLines(paste0("updating header in ", pathedOneRfile))
            writeLines(c(headerLines, allLines),
                       con = pathedOneRfile)
        }
    }
}
