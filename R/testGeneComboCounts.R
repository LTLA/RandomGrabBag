#' Test for DE gene combinations
#'
#' Test for differential expression of gene combinations between pairs of groups using Fisher's exact test.
#'
#' @inheritParams summarizeGeneComboCounts
#' @param ... Further arguments to pass to \code{\link[scran]{combineMarkers}}.
#'
#' @author Aaron Lun
#'
#' @return A \linkS4class{List} containing one \linkS4class{DataFrame} per group.
#' Each row of a DataFrame for a particular group corresponds to a gene combination, 
#' ordered by the degree of upregulation of that combination in the current group compared to all other groups.
#' Columns contain the combined p-value across all pairwise comparisons involving the current group
#' as well as the log-odds ratio against each other group.
#'
#' @details
#' This performs pairwise tests for differential expression of each gene combination in each group compared to every other group.
#' For each gene combination and for each pair of groups,
#' we construct a 2-by-2 contigency matrix using the number of cells expressing that combination in the two groups.
#' We condition on the total number of cells - see \code{?\link{countCellsPerGeneCombo}} for comments on cell count normalization
#' - and then apply a one-sided Fisher's exact test.
#' 
#' Once we obtain p-values and log-odds ratios for all pairwise comparisons, 
#' we consolidate them into a ranked list for each group using the \code{\link[scran]{combineMarkers}} function.
#' This provides a convenient per-group summary of the genes upregulated in each group compared to all others.
#' See \code{?\link[scran]{combineMarkers}} for more details.
#'
#' @examples
#' df <- data.frame(
#'     cell.id=sample(LETTERS, 30, replace=TRUE),
#'     v_gene=sample(c("TRAV1", "TRAV2", "TRAV3"), 30, replace=TRUE),
#'     j_gene=sample(c("TRAJ4", "TRAJ5", "TRAV6"), 30, replace=TRUE)
#' )
#'
#' y <- splitToCells(df, field="cell.id")
#' out <- countCellsPerGeneCombo(y, c("v_gene", "j_gene"), 
#'    group=sample(3, length(y), replace=TRUE))
#'
#' de <- testGeneComboCountsPairwise(out)
#' de[[1]]
#'
#' @export
#' @importFrom stats fisher.test
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors DataFrame 
testGeneComboCountsPairwise <- function(counts, ...) {
    idx <- seq_len(ncol(counts))
    nm <- colnames(counts)
    if (is.null(nm)) nm <- idx

    if (is(counts, "SummarizedExperiment")) {
        counts <- assay(counts)
    }
    totals <- colSums(counts)

    collected.stats <- collected.pairs <- list()
    counter <- 1L

    for (i in idx) {
        for (j in idx) {
            if (i==j) next

            p <- LOR <- numeric(nrow(counts))
            for (k in seq_along(p)) {
                y <- counts[k,c(i,j)]
                test.out <- fisher.test(rbind(y, totals[c(i,j)]-y), alternative="greater")
                p[k] <- test.out$p.value
                LOR[k] <- log2(test.out$estimate)
            }

            collected.pairs[[counter]] <- DataFrame(first=nm[i], second=nm[j])
            collected.stats[[counter]] <- DataFrame(p.value=p, LOR=LOR, row.names=rownames(counts))
            counter <- counter + 1L
        }
    }

    all.pairs <- do.call(rbind, collected.pairs)
    o <- order(all.pairs$first, all.pairs$second)
    scran::combineMarkers(collected.stats[o], all.pairs[o,,drop=FALSE], effect.field="LOR", ...)
}
