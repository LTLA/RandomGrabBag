#' Summarize gene combination counts 
#'
#' Generate some summary statistics for the gene diversity, based on the number of cells expressing each gene combination.
#'
#' @param counts A \linkS4class{SummarizedExperiment} containing cell counts for each gene combination (row) and group (column),
#' such as that produced by \code{\link{countCellsPerGeneCombo}}.
#' Alternatively, a count matrix containing the same information.
#' @inheritParams summarizeClonotypeCounts
#' 
#' @return A \linkS4class{DataFrame} with one row per group in \code{counts},
#' containing summary statistics on the diversity of gene expression in that group.
#'
#' @author Aaron Lun
#'
#' @details
#' If \code{use.gini=TRUE}, the output will contain the numeric \code{"gini"} column,
#' containing the Gini index for gene combination diversity in each group.
#' Larger values indicate that a small number of gene combinations are expressed in many cells.
#'
#' If \code{use.top} is specified, the output will contain the numeric \code{"topX"} columns,
#' containing the proportion of cells expressing the top X gene combinations (where X is the value of each \code{use.top} entry).
#' Larger values indicate that a small number of gene combinations are expressed in many cells.
#'
#' If \code{use.hill} is specified, the output will contain numeric \code{"hillX"} columns containing the Hill numbers.
#' \code{"hill0"} is simply the number of observed gene combinations (i.e., species richness)
#' while \code{"hill1"} and \code{"hill2"} quantify the evenness of the distribution of cells across gene combinations.
#' (\code{"hill2"} differs from \code{"hill1"} in that the former gives more weight to dominant gene combinations.)
#'
#' See \code{\link{countCellsPerGeneCombo}} for an explanation of why downsampling is turned on by default.
#' @examples
#' df <- data.frame(
#'     cell.id=sample(LETTERS, 30, replace=TRUE),
#'     v_gene=sample(c("TRAV1", "TRAV2", "TRAV3"), 30, replace=TRUE),
#'     j_gene=sample(c("TRAJ4", "TRAJ5", "TRAV6"), 30, replace=TRUE)
#' )
#'
#' y <- splitDataFrameByCell(df, field="cell.id")
#' out <- countCellsPerGeneCombo(y, c("v_gene", "j_gene"),
#'    group=sample(10, length(y), replace=TRUE))
#'
#' summarizeGeneComboCounts(out)
#'
#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom alakazam calcDiversity
#' @importFrom SummarizedExperiment assay
summarizeGeneComboCounts <- function(counts, use.gini=TRUE, use.top=c(5, 20, 100), use.hill=0:2,
    downsample=TRUE, down.ncells=NULL)
{
    if (is(counts, "SummarizedExperiment")) {
        counts <- assay(counts)
    }
    if (downsample) {
        counts <- .downsample_matrix(counts, down.ncells)
    }

    stats <- DataFrame(row.names=colnames(counts))
    idx <- seq_len(ncol(counts))

    if (use.gini) {
        stats[["gini"]] <- vapply(idx, FUN=function(j) .compute_gini(counts[,j]), FUN.VALUE=0)
    }

    for (n in use.top) {
        stats[[paste0("top", n)]] <- vapply(idx, FUN=function(j) .compute_top_prop(counts[,j], n), FUN.VALUE=0)
    }

    if (length(use.hill)) {
        hill <- lapply(idx, function(j) calcDiversity(counts[,j], q=use.hill))
        hill <- do.call(rbind, hill)
        colnames(hill) <- sprintf("hill%i", use.hill)
        stats <- cbind(stats, hill)
    }

    stats
}
