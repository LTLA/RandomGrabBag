#' Read to UMI per cell
#'
#' Compute the read-to-UMI ratio for each cell.
#'
#' @param x A \linkS4class{SplitDataFrameList} where each DataFrame is a cell and each row is a sequence.
#' @param read.field String containing the name of the column containing the read count data.
#' @param umi.field String containing the name of the column containing the UMI count data.
#'
#' @author Aaron Lun
#'
#' @return A numeric vector containing the ratio of reads to UMI for each cell.
#'
#' @details
#' This function is designed to evaluate the degree of redundancy in the read coverage of each UMI.
#' High values indicate that the reads are highly redundant such that little can be gained from further sequencing.
#'
#' Note that, in repertoire data, the definition of \dQuote{high} is somewhat different from usual.
#' This is because only deeply sequenced transcripts will survive the assembly and annotation process,
#' such that the reported sequences are likely to be biased towards very high read-to-UMI ratios.
#' Values around 1000 seem to be typical.
#'
#' If a cell has multiple sequences, their counts are simply added together across sequences to compute the per-cell ratio.
#' 
#' @examples
#' df <- data.frame(
#'     cell.id=sample(LETTERS, 30, replace=TRUE),
#'     v_gene=sample(c("TRAV1", "TRAV2", "TRAV3"), 30, replace=TRUE),
#'     j_gene=sample(c("TRAJ4", "TRAJ5", "TRAV6"), 30, replace=TRUE),
#'     reads=rnbinom(30, mu=20, size=0.5),
#'     umis=rnbinom(30, mu=2, size=1)
#' )
#'
#' y <- splitDataFrameByCell(df, field="cell.id")
#' readToUmiPerCell(y, "reads", "umis")
#' 
#' @export
readToUmiPerCell <- function(x, read.field, umi.field) {
    sum(x[,read.field])/sum(x[,umi.field])
}
