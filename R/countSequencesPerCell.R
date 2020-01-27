#' Count sequences per cell
#'
#' Count the number of sequences per cell, possibly after filtering.
#'
#' @param x A \linkS4class{SplitDataFrameList} where each \linkS4class{DataFrame} is a cell and each row is a sequence.
#' @param filter.field Character vector specifying the columns on which to filter sequences prior to counting.
#' @param filter.value Character vector of length equal to \code{filter.field},
#' specifying the values to retain for each filter field.
#'
#' @details
#' The number of sequences per cell is often a useful diagnostic.
#' At its simplest, we can use it to determine whether a particular cell contributes to the immune repertoire at all,
#' e.g., to verify clusters that are B or T cells.
#'
#' A more complex use case is to identify cells that express multiple sequences.
#' This is generally a minority occurrence due to allelic exclusion in most cells (see also \code{\link{topCoveragePropPerCell}})
#' but can be inflated by technical artifacts such as doublets or contamination from ambient noise.
#'
#' The filtering enables us to perform more complex diagnostics,
#' e.g., count the number of productive, full-length, high-quality sequences in each cell.
#'
#' @return An integer scalar containing the number of (filtered) sequences per cell.
#'
#' @author Aaron Lun
#' 
#' @examples
#' df <- data.frame(
#'     cell.id=sample(LETTERS, 30, replace=TRUE),
#'     v_gene=sample(c("TRAV1", "TRAV2", "TRAV3"), 30, replace=TRUE),
#'     j_gene=sample(c("TRAJ4", "TRAJ5", "TRAV6"), 30, replace=TRUE),
#'     productive=sample(c("True", "False"), 30, replace=TRUE)
#' )
#'
#' y <- splitDataFrameByCell(df, field="cell.id")
#' countSequencesPerCell(y)
#' countSequencesPerCell(y, filter.field="productive", filter.value="True")
#' 
#' @export
countSequencesPerCell <- function(x, filter.field=NULL, filter.value=NULL) {
    if (length(filter.field)) {
        keep <- TRUE
        for (i in seq_along(filter.field)) {
            keep <- keep & x[,filter.field[[i]]]==filter.value[[i]]
        }
        sum(keep)
    } else {
        lengths(x) 
    }
}

