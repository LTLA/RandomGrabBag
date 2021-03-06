#' Create a SplitDataFrameList
#'
#' Create a \linkS4class{SplitDataFrameList} to store repertoire data.
#' Each cell corresponds to one \link{DataFrame}, while each row of that DataFrame corresponds to a sequence detected in the cell.
#' 
#' @param df A data frame of repertoire data where each row corresponds to a detected sequence.
#' @param field A string specifying the field of \code{df} that contains the identity of the cell.
#' @param ids A character vector of length equal to \code{nrow(df)} specifying the identity of the cell for each sequence.
#' @param universe Character vector of all possible cell identities.
#'
#' @details
#' A cell may have anywhere from zero, one or multiple sequences in a repertoire sequencing experiment.
#' This makes repertoire data rather inconvenient to store and manipulate;
#' at some points, we would like to operate on cells, while at other points, we would like to operate on the raw sequences.
#'
#' Our solution is to use a SplitDataFrameList where each entry corresponds to a single cell.
#' However, each entry is also a DataFrame where each row corresponds to a sequence in that cell.
#' This achieves a per-cell representation without discarding per-sequence information.
#'
#' We can further take advantage of powerful \linkS4class{List} grammar to perform per-sequence operations on the SplitDataFrameList; see the relevant documentation from the \pkg{IRanges} package for more details.
#'
#' Setting \code{universe} is useful to ensure that the output object is of the same length as,
#' say, the number of columns in a \linkS4class{SingleCellExperiment} object containing expression data.
#'
#' @return
#' A \linkS4class{SplitDataFrameList} containing a per-cell perspective of repertoire data.
#'
#' @author Aaron Lun
#'
#' @examples
#' df <- data.frame(
#'     cell.id=sample(LETTERS, 30, replace=TRUE),
#'     v_gene=sample(c("TRAV1", "TRAV2", "TRAV3"), 30, replace=TRUE),
#'     j_gene=sample(c("TRAJ4", "TRAJ5", "TRAV6"), 30, replace=TRUE)
#' )
#'
#' splitDataFrameByCell(df, field="cell.id")
#'
#' splitDataFrameByCell(df, field="cell.id", universe=c(LETTERS, "AA"))
#'     
#' @export
#' @importFrom S4Vectors split DataFrame
splitDataFrameByCell <- function(df, field, ids=df[[field]], universe=sort(unique(ids))) {
    split(DataFrame(df), factor(ids, universe))
}
