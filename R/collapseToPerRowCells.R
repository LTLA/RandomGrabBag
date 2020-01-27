#' Collapse to one cell per row
#'
#' Collapse a \linkS4class{SplitDataFrameList} representation into a DataFrame with one row per cell.
#'
#' @param x A \linkS4class{SplitDataFrameList} object containing one DataFrame per cell,
#' where each row of each DataFrame contains information for one sequence in that cell.
#' @param cov.field String specifying the column of \code{x} containing the UMI/read count per sequence.
#' @param fill Logical scalar indicating whether cells with no sequences should be filled in with \code{NA} in the output.
#'
#' @return A \linkS4class{DataFrame} with one row for each cell.
#'
#' @details 
#' This function collapses the SplitDataFrameList into a DataFrame such that each cell is represented by exactly one row.
#' If a cell has multiple sequences, one representative sequence is chosen:
#' \itemize{
#' \item If \code{cov.field} is specified, the sequence with the largest count is selected.
#' This favors the dominant sequence with the highest number of captured molecules.
#' \item Otherwise, the first sequence for each cell is selected.
#' This is effectively an arbitrary choice as the ordering of sequences has no meaning.
#' }
#' 
#' If a cell has no sequences, the output is filled in with \code{NA} if \code{fill=TRUE}.
#' Otherwise, it is simply not reported.
#'
#' @author Aaron Lun
#'
#' @examples
#' df <- data.frame(
#'     cell.id=sample(LETTERS, 30, replace=TRUE),
#'     v_gene=sample(c("TRAV1", "TRAV2", "TRAV3"), 30, replace=TRUE),
#'     j_gene=sample(c("TRAJ4", "TRAJ5", "TRAV6"), 30, replace=TRUE),
#'     umi=pmax(1, rpois(30, 1))
#' )
#'
#' Y <- splitToCells(df, "cell.id")
#' collapseToPerCellRows(Y, "umi")
#'
#' collapseToPerCellRows(Y)
#'
#' collapseToPerCellRows(Y, fill=FALSE)
#' 
#' @export
#' @importFrom utils head
collapseToPerCellRows <- function(x, cov.field, fill=TRUE) {
    if (!missing(cov.field)) {
        cov <- x[,cov.field]
        y <- unlist(x[cov==max(cov)])
    } else {
        # Getting the first element for each cell.
        keep <- c(1L, head(cumsum(lengths(x)), -1L)+1L)
        y <- unlist(x)[keep,,drop=FALSE]
    }

    if (fill) {
        expander <- match(seq_along(x), which(lengths(x)!=0L))
        y <- y[expander,,drop=FALSE]
        rownames(y) <- names(x)
    } else {
        rownames(y) <- names(x)[lengths(x)!=0L]
    }

    y
}
