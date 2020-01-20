#' Count gene combinations
#'
#' Count the number of cells that express each unique combination of genes.
#'
#' @param x A \linkS4class{SplitDataFrameList} where each \linkS4class{DataFrame} is a cell and each row is a sequence.
#' @param gene.field Character vector of names of columns of \code{x} containing the genes of interest (e.g., VDJ components).
#' @param group Factor of length equal to \code{x} indicating the group to which each cell belongs.
#' @param cov.field String specifying the column of \code{x} containing the read/UMI coverage.
#'
#' @return A \linkS4class{SummarizedExperiment} where each row corresponds to a unique gene combination
#' and each column corresponds to a level of \code{group} (or all cells, if \code{group=NULL}).
#' The \code{assays} contain a single matrix containing the number of cells for each gene combination and grouping level,
#' while the \code{rowData} contains information about the gene combination.
#'
#' @author Aaron Lun
#'
#' @details
#' The aim of this function is to generate a count matrix for use in differential \dQuote{expression} analyses,
#' i.e., does one particular group of cells express a particular gene combination more frequently than another group? 
#' This can be useful to examine the effect of particular experimental conditions or the behavior of different cell states,
#' especially if the specific biological function (e.g., antigen) of each gene combination is known in advance.
#'
#' We quote \dQuote{expression} as this is defined in terms of number of cells expressing the gene,
#' rather than the more typical quantity of the number of reads or UMIs assigned to that gene.
#' Note that this raises some difficult questions about normalization when the sequencing coverage varies between groups;
#' we assume that such changes have the same scaling effect on the probability of detecting each gene combination,
#' such that it cancels out after normalizing by the total number of cells.
#'
#' If \code{cov.field} is set, only the most high-abundance sequence is used from each cell.
#' It is probably safest to set this field to avoid potential complications from dependencies between counts,
#' though any problems are also probably minor.
#'
#' @examples
#' df <- data.frame(
#'     cell.id=sample(LETTERS, 30, replace=TRUE),
#'     v_gene=sample(c("TRAV1", "TRAV2", "TRAV3"), 30, replace=TRUE),
#'     j_gene=sample(c("TRAJ4", "TRAJ5", "TRAV6"), 30, replace=TRUE)
#' )
#'
#' y <- splitToCells(df, field="cell.id")
#' out <- countCellsPerGeneCombo(y, c("v_gene", "j_gene"))
#' rowData(out)
#' assay(out)
#'
#' out2 <- countCellsPerGeneCombo(y, c("v_gene", "j_gene"), 
#'    group=sample(10, length(y), replace=TRUE))
#' rowData(out2)
#' assay(out2)
#'
#' @export
#' @importFrom S4Vectors selfmatch
#' @importFrom IRanges IntegerList
countCellsPerGeneCombo <- function(x, gene.field, group=NULL, cov.field=NULL) {
    if (!is.null(cov.field)) {
        cov <- x[,cov.field]
        x <- x[cov==max(cov)]
    }

    y <- unlist(x)[,gene.field]
    ids <- selfmatch(y)
    keep <- !duplicated(ids)
    y <- y[keep,]

    if (is.null(group)) {
        mat <- table(ids)
    } else {
        mat <- table(ids, rep(group, lengths(x)))
    }
    mat <- as.matrix(mat)

    rownames(mat) <- NULL
    rownames(y) <- NULL
    SummarizedExperiment(rowData=y, assays=list(counts=mat))
}
