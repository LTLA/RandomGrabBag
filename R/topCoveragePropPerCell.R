#' Compute the top coverage proportion per cell
#'
#' For each cell, compute the proportion of reads/UMIs that are taken up by the most abundant sequence.
#'
#' @param x A \linkS4class{SplitDataFrameList} where each DataFrame is a cell and each row is a sequence.
#' @param cov.field String specifying the column containing the coverage data, usually read or UMI counts.
#' @param second.ratio Logical scalar indicating whether the ratio of coverage of the second-most-abundance sequence
#' to the most abundant sequence should instead be computed.
#'
#' @return A numeric vector containing the per-cell proportion of reads/UMIs taken up by that cell's most abundant sequence.
#'
#' If \code{second.ratio=TRUE}, the vector instead contains the ratio of the reads/UMIs taken up by that cell's second-most
#' abundant sequence, compared to that taken by the most abundant sequence.
#'
#' @details
#' This function is designed to summarize the distribution of sequence abundances within a cell.
#' A proportion close to unity indicates that there is only one dominant sequence for that cell,
#' while lower proportions suggest that the secondary sequences have comparable expression to the top sequence.
#' The presence of low proportions may be interesting as most cells should undergo allelic exclusion for TCR and Ig components. 
#'
#' If \code{second.ratio=TRUE}, the ratio of the second-most-abundant to the most-abundant sequence is returned instead.
#' This may provide a more intuitive representation of the relevance of the secondary sequences,
#' with larger values indicating that the secondary sequence is expressed at a comparable level to the dominant sequence.
#' It is also unaffected by the presence of further sequences, e.g., due to ambient contamination.
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
#'     umis=pmin(1, rnbinom(30, mu=2, size=1))
#' )
#'
#' y <- splitToCells(df, field="cell.id")
#' topCoveragePropPerCell(y, "umis")
#' @export
topCoveragePropPerCell <- function(x, cov.field, second.ratio=FALSE) {
    counts <- x[,cov.field]
    max <- max(counts)

    if (second.ratio) {
        second <- vapply(counts, function(y) sort(y, decreasing=TRUE)[2], 0)
        out <- second/max
        out[max==0] <- NA_real_

    } else {
        totals <- sum(counts)
        out <- max/totals
        out[totals==0] <- NA_real_
    }

    out
}


