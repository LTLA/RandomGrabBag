#' Summarize clonal expansion
#'
#' Generate some summary statistics for clonal expansion based on the number of cells per clonotype.
#'
#' @param counts A list of integer vectors such as that produced by \code{\link{countCellsPerClonotype}}.
#' Each vector corresponds to a group of cells and contains the number of cells for each clonotype in that group.
#' @param use.mean Logical scalar indicating whether to report the mean number of cells per clonotype.
#' @param use.gini Logical scalar indicating whether to report the Gini index.
#' @param use.top Integer vector specifying the number of clonotypes to use to compute the top percentage.
#' @param use.hill Integer scalar specifying the orders to use to compute Hill numbers.
#' 
#' @return A \linkS4class{DataFrame} with one row per group in \code{counts},
#' containing summary statistics on the clonal expansion in that group.
#'
#' @author Aaron Lun
#'
#' @details
#' If \code{use.mean=TRUE}, the output will contain the numeric \code{"mean"} column,
#' containing the average number of cells per clonotype computed across all clonotypes for a given group.
#' Larger values indicate that cells are concentrated within a small number of clonotypes.
#'
#' If \code{use.gini=TRUE}, the output will contain the numeric \code{"gini"} column,
#' containing the Ginit index for clonotype diversity in each group.
#' Larger values indicate that cells are concentrated within a small number of clonotypes.
#'
#' If \code{use.top} is specified, the output will contain the numeric \code{"topX"} columns,
#' containing the proportion of cells assigned to the top X clonotypes (where X is the value of each \code{use.top} entry).
#' Larger values indicate that cells are concentrated within a small number of clonotypes.
#'
#' If \code{use.hill} is specified, the output will contain numeric \code{"hillX"} columns containing the Hill numbers.
#' \code{"hill0"} is simply the number of observed clonotypes (i.e., species richness)
#' while \code{"hill1"} and \code{"hill2"} quantify the evenness of the distribution of cells across clonotypes.
#' (\code{"hill2"} differs from \code{"hill1"} in that the former gives more weight to dominant clonotypes.)
#' 
#' @examples
#' df <- data.frame(
#'     cell.id=sample(LETTERS, 30, replace=TRUE),
#'     clonotype=sample(paste0("clonotype_", 1:5), 30, replace=TRUE)
#' )
#' 
#' y <- splitToCells(df, field="cell.id")
#' out <- countCellsPerClonotype(y, "clonotype",
#'    group=sample(3, length(y), replace=TRUE))
#'
#' summarizeClonalExpansion(out)
#' 
#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom alakazam calcDiversity
summarizeClonalExpansion <- function(counts, use.mean=TRUE, use.gini=TRUE, use.top=c(5, 20, 100), use.hill=0:2) {
    stats <- DataFrame(row.names=names(counts))

    if (use.mean) {
        stats[["mean"]] <- vapply(counts, FUN=mean, FUN.VALUE=0)
    }

    if (use.gini) {
        # Discrete version of the gini index.
        stats[["gini"]] <- vapply(counts, FUN.VALUE=0, FUN=function(y) {
            total <- sum(y)
            1- sum(cumsum(sort(y))) / (total * (total + 1)/2)
        })
    }

    for (i in use.top) {
        stats[[paste0("top", i)]] <- vapply(counts, FUN.VALUE=0, FUN=function(y) {
            sum(head(sort(y, decreasing=TRUE), i))/sum(y) 
        })
    }

    if (length(use.hill)) {
        hill <- lapply(counts, calcDiversity, q=use.hill)
        hill <- do.call(rbind, hill)
        colnames(hill) <- sprintf("hill%i", use.hill)
        stats <- cbind(stats, hill)
    }

    stats
}
