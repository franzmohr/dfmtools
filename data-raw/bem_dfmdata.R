## How data/bem_dfmdata.rda is produced.
##
## inst/COPYRIGHTS points at this file as the record of where the object came
## from, so it is kept even though the two inputs are not redistributed: they
## are third-party data this package has no licence to pass on, and .gitignore
## keeps them out while keeping this script in.
##
##   bem_dfmdata.csv  the quarterly panel of the data sets that accompany Chan,
##                    Koop, Poirier and Tobias (2019), headerless.
##   fred_qd.csv      FRED-QD as those data sets ship it. Only its header is
##                    used, for the series names; its first column is the date
##                    and its first two rows are the transformation codes.
##
## The two describe the same columns in the same order, which is the one thing
## here that cannot be checked from the files themselves and is asserted below:
## a mismatch would name all 196 series wrongly and nothing downstream would
## notice.

x <- read.csv(file = "data-raw/bem_dfmdata.csv", header = FALSE)
fred <- read.csv(file = "data-raw/fred_qd.csv")[-c(1, 2), -1]

stopifnot(ncol(x) == ncol(fred))

## Columns with a complete history. The panel starts in 1959Q3 and several
## series begin later, so they are dropped rather than carried as gaps.
pos <- which(vapply(x, function(column) !anyNA(column), logical(1)))

bem_dfmdata <- ts(x[, pos], start = c(1959, 3), frequency = 4,
                  names = names(fred)[pos])

stopifnot(ncol(bem_dfmdata) == length(pos), !anyNA(bem_dfmdata))

usethis::use_data(bem_dfmdata, overwrite = TRUE)
