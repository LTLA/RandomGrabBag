#' Count gene combinations
#'
#' Count the number of cells that express each unique combination of genes.
#'
#' @param x 
#' Any data.frame-like object where each row corresponds to a single cell and contains its representative sequence.
#' Rows with any \code{NA} values in the specified \code{gene.field} columns are ignored.
#' 
#' Alternatively, a \linkS4class{SplitDataFrameList} where each \linkS4class{DataFrame} corresponds to a cell 
#' and each row in that DataFrame is a sequence in that cell.
#' @param gene.field Character vector of names of columns of \code{x} containing the genes of interest (e.g., VDJ components).
#' @param group Factor of length equal to \code{x} indicating the group to which each cell belongs.
#' @param cov.field String specifying the column of \code{x} containing the read/UMI coverage.
#' @param downsample Logical scalar indicating whether downsampling should be performed.
#' @param down.ncells Integer scalar indicating the number of cells to downsample each group to.
#' Defaults to the number of cells in the smallest group in \code{group}.
#' @param row.names Logical scalar indicating whether row names should be added by concatenating all gene names per combination.
#' @param ... For the generic, further arguments to pass to individual methods.
#'
#' For the \code{CompressedSplitDataFrameList} method, further arguments to pass to the ANY method.
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
#' If \code{cov.field} is set, only the most high-abundance sequence is used from each cell.
#' In contrast, setting \code{cov.field=NULL} will count each sequence separately,
#' such that one cell may contribute multiple times.
#' It is probably safest to set this to some non-\code{NULL} value to avoid complications from dependencies between counts,
#' though any problems are also probably minor.
#'
#' @section Normalization for cell number:
#' Here, expression is defined in terms of number of cells expressing the gene,
#' rather than the more typical quantity of the number of reads or UMIs assigned to that gene.
#' If the sequencing coverage varies between groups, 
#' we assume that such changes have the same scaling effect on the probability of detecting each gene combination,
#' which cancels out after normalizing by the total number of cells.
#'
#' However, the above assumption only works for differential expression analyses between groups.
#' When comparing other metrics such as diversity values (see \code{\link{summarizeGeneComboCounts}}),
#' scaling normalization is not sufficient and we instead resort to downsampling all groups to the same total cell number.
#' This is achieved with \code{downsample=TRUE} with the automatically determined \code{down.ncells},
#' which eliminates uninteresting technical differences between groups from cell capture efficiency or sample size.
#'
#' @examples
#' df <- data.frame(
#'     cell.id=sample(LETTERS, 30, replace=TRUE),
#'     v_gene=sample(c("TRAV1", "TRAV2", "TRAV3"), 30, replace=TRUE),
#'     j_gene=sample(c("TRAJ4", "TRAJ5", "TRAV6"), 30, replace=TRUE),
#'     umi=pmax(1, rpois(30, 1))
#' )
#'
#' y <- splitDataFrameByCell(df, field="cell.id")
#' out <- countCellsPerGeneCombo(y, c("v_gene", "j_gene"), cov.field="umi")
#' rowData(out)
#' assay(out)
#'
#' out2 <- countCellsPerGeneCombo(y, c("v_gene", "j_gene"), cov.field="umi",
#'    group=sample(10, length(y), replace=TRUE))
#' rowData(out2)
#' assay(out2)
#'
#' @export
#' @name countCellsPerGeneCombo
setGeneric("countCellsPerGeneCombo", function(x, ...) standardGeneric("countCellsPerGeneCombo"))

#' @importFrom S4Vectors selfmatch
#' @importFrom IRanges IntegerList
#' @importFrom SummarizedExperiment SummarizedExperiment
.countCellsPerGeneCombo <- function(x, gene.field, group=NULL,
    downsample=FALSE, down.ncells=NULL, row.names=TRUE) 
{
    y <- x[,gene.field]
    discard <- Reduce("|", lapply(y, is.na))
    if (any(discard)) {
        y <- y[discard,]
        group <- group[discard]
    }

    ids <- selfmatch(y)
    keep <- !duplicated(ids)
    y <- y[keep,]

    if (is.null(group)) {
        mat <- table(ids)
    } else {
        mat <- table(ids, group)
    }
    mat <- as.matrix(mat)

    if (downsample) {
        mat <- .downsample_matrix(mat, down.ncells)
    }

    if (row.names) {
        rn <- do.call(paste, c(as.list(y), list(sep="|")))
    } else {
        rn <- NULL
    }
    rownames(mat) <- rownames(y) <- rn

    SummarizedExperiment(rowData=y, assays=list(counts=mat))
}

.downsample_matrix <- function(mat, down.ncells) {
    if (is.null(down.ncells)) {
        down.ncells <- min(colSums(mat))
    }
    for (j in seq_len(ncol(mat))) {
        mat[,j] <- .downsample_counts(mat[,j], down.ncells)
    }
    mat
}

#' @export
#' @rdname countCellsPerGeneCombo
setMethod("countCellsPerGeneCombo", "ANY", .countCellsPerGeneCombo)

#' @export
#' @rdname countCellsPerGeneCombo
#' @importClassesFrom IRanges CompressedSplitDataFrameList
setMethod("countCellsPerGeneCombo", "CompressedSplitDataFrameList", function(x, gene.field, cov.field, group=NULL, ...) {
    if (!is.null(cov.field)) {
        cov <- x[,cov.field]
        x <- x[cov==max(cov)]
    }

    if (!is.null(group)) {
        group <- group[rep(seq_along(x), lengths(x))]
    }

    .countCellsPerGeneCombo(unlist(x), gene.field=gene.field, group=group, ...)
})
