#' Compute the top coverage proportion per cell
#'
#' For each cell, compute the proportion of reads/UMIs that are taken up by the most abundant sequence.
#'
#' @param x A \linkS4class{SplitDataFrameList} where each DataFrame is a cell and each row is a sequence.
#' @param cov.field String specifying the column containing the coverage data, usually read or UMI counts.
#'
#' @return A numeric vector containing the per-cell proportion of reads/UMIs taken up by that cell's most abundant sequence.
#'
#' @details
#' This function is designed to summarize the distribution of sequence abundances within a cell.
#' A proportion close to unity indicates that there is only one dominant sequence for that cell,
#' while lower proportions suggest that the secondary sequences have comparable expression to the top sequence.
#' The presence of low proportions may be interesting as most cells should undergo allelic exclusion for TCR and Ig components. 
#' 
#' The UMI count is usually the more relevant metric as it avoids amplification biases
#' that add noise to comparisons between counts of different sequences.
#'
#' @author Aaron Lun
#'
#' @examples
#' df <- data.frame(
#'     cell.id=sample(LETTERS, 30, replace=TRUE),
#'     v_gene=sample(c("TRAV1", "TRAV2", "TRAV3"), 30, replace=TRUE),
#'     j_gene=sample(c("TRAJ4", "TRAJ5", "TRAV6"), 30, replace=TRUE),
#'     umis=rnbinom(30, mu=2, size=1)
#' )
#'
#' y <- splitToCells(df, field="cell.id")
#' topCoveragePropPerCell(y, "umis")
#' @export
topCoveragePropPerCell <- function(x, cov.field) {
    counts <- x[,cov.field]
    totals <- sum(counts)
    max <- max(counts)

    out <- max/totals
    out[totals==0] <- NA_real_
    out
}


