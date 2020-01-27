#' Summarize clonotype counts 
#'
#' Generate some summary statistics for clonal expansion based on the number of cells per clonotype.
#'
#' @param counts A list of integer vectors such as that produced by \code{\link{countCellsPerClonotype}}.
#' Each vector corresponds to a group of cells and contains the number of cells for each clonotype in that group.
#' @param use.mean Logical scalar indicating whether to report the mean number of cells per clonotype.
#' @param use.gini Logical scalar indicating whether to report the Gini index.
#' @param use.top Integer vector specifying the number of clonotypes to use to compute the top percentage.
#' @param use.hill Integer scalar specifying the orders to use to compute Hill numbers.
#' @inheritParams countCellsPerClonotype
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
#' containing the Gini index for clonotype diversity in each group.
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
#' See \code{\link{countClonotypesPerCell}} for why downsampling is turned on by default.
#' @examples
#' df <- data.frame(
#'     cell.id=sample(LETTERS, 30, replace=TRUE),
#'     clonotype=sample(paste0("clonotype_", 1:5), 30, replace=TRUE)
#' )
#' 
#' y <- splitDataFrameByCell(df, field="cell.id")
#' out <- countCellsPerClonotype(y, "clonotype",
#'    group=sample(3, length(y), replace=TRUE))
#'
#' summarizeClonotypeCounts(out)
#' 
#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom alakazam calcDiversity
summarizeClonotypeCounts <- function(counts, use.mean=TRUE, use.gini=TRUE, 
    use.top=c(5, 20, 100), use.hill=0:2,
    downsample=TRUE, down.ncells=NULL)
{
    if (downsample) {
        counts <- .downsample_list(counts, down.ncells)
    }

    stats <- DataFrame(row.names=names(counts))

    if (use.mean) {
        stats[["mean"]] <- vapply(counts, FUN=mean, FUN.VALUE=0)
    }

    if (use.gini) {
        stats[["gini"]] <- vapply(counts, FUN=.compute_gini, FUN.VALUE=0)
    }

    for (n in use.top) {
        stats[[paste0("top", n)]] <- vapply(counts, FUN=.compute_top_prop, n=n, FUN.VALUE=0)
    }

    if (length(use.hill)) {
        hill <- lapply(counts, calcDiversity, q=use.hill)
        hill <- do.call(rbind, hill)
        colnames(hill) <- sprintf("hill%i", use.hill)
        stats <- cbind(stats, hill)
    }

    stats
}

# Discrete version of the gini index.
.compute_gini <- function(y) {
    total <- sum(y)
    1- sum(cumsum(sort(y))) / (total * (total + 1)/2)
}

#' @importFrom utils head
.compute_top_prop <- function(y, n) {
    sum(head(sort(y, decreasing=TRUE), n))/sum(y) 
}
