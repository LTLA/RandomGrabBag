#' Summarize sequence status
#'
#' Obtain a quick summary of the status of the sequences for a particular component.
#'
#' @param x A \linkS4class{SplitDataFrameList} where each \linkS4class{DataFrame} corresponds to a cell 
#' and each row in that DataFrame is a sequence in that cell.
#' @param group Factor of length equal to \code{x} indicating the group to which each cell belongs.
#'
#' @return If \code{group=NULL}, a \linkS4class{DataFrame} is returned with one row per cell,
#' indicating whether that cell has:
#' \itemize{
#' \item any sequence
#' \item multiple sequences
#' \item any productive sequence
#' \item any full-length sequence
#' \item any high-confidence sequence
#' \item any awesome (productive, full-length and high-confidence) sequence
#' }
#'
#' If \code{group} is specified, the DataFrame instead contains one row per level of \code{group}.
#' Each value then represents the proportion of cells in that group with any sequence, multiple sequences, etc.
#'
#' @details
#' By default, this assumes that the fields are named in the same manner as the annotation files produced by CellRanger.
#' Future iterations will provide support for more standardized formats like AIRR.
#' 
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{countSequencesPerCell}}, for which this function is a wrapper.
#'
#' @examples
#' df <- data.frame(
#'     cell.id=sample(LETTERS, 30, replace=TRUE),
#'     clonotype=sample(paste0("clonotype_", 1:5), 30, replace=TRUE),
#'     full_length=sample(c("True", "False"), 30, replace=TRUE),
#'     high_confidence=sample(c("True", "False"), 30, replace=TRUE),
#'     productive=sample(c("True", "False"), 30, replace=TRUE),
#'     umi=pmax(1, rpois(30, 5))
#' )
#'
#' Y <- splitDataFrameByCell(df, "cell.id")
#' summarizeSequenceStatus(Y)
#'
#' summarizeSequenceStatus(Y, group=sample(1:3, nrow(df), replace=TRUE))
#'
#' @export
#' @importFrom S4Vectors DataFrame
summarizeSequenceStatus <- function(x, group=NULL) {
    output <- DataFrame(
        Any=countSequencesPerCell(x) > 0,
        Multi=countSequencesPerCell(x) > 1,
        Productive=countSequencesPerCell(x,
            filter.field="productive", filter.value="True") > 0,
        FullLength=countSequencesPerCell(x,
            filter.field="full_length", filter.value="True") > 0,
        HighConf=countSequencesPerCell(x,
            filter.field="high_confidence", filter.value="True") > 0,
        Awesome=countSequencesPerCell(x,
            filter.field=c("productive", "full_length", "high_confidence"),
            filter.value=c("True", "True", "True")) > 0
    )

    if (!is.null(group)) {
        by.cluster <- split(output, group)
        props <- lapply(by.cluster, function(x) colMeans(as.matrix(x)))
        props <- do.call(rbind, props)
        DataFrame(props)
    } else {
        rownames(output) <- names(x)
        output
    }
}
