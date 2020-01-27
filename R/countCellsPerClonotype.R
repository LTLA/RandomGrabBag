#' Count cells per clonotype
#' 
#' Count the number of cells that exhibit each clonotype.
#'
#' @param x 
#' Any data.frame-like object where each row corresponds to a single cell and contains its representative sequence.
#' Rows with any \code{NA} values in the specified \code{clone.field} columns are ignored.
#' 
#' Alternatively, a \linkS4class{SplitDataFrameList} where each \linkS4class{DataFrame} corresponds to a cell 
#' and each row in that DataFrame is a sequence in that cell.
#' @param clone.field String specifying the columns of \code{x} containing the clonotype identity.
#' @param group Factor of length equal to \code{x} indicating the group to which each cell belongs.
#' @param cov.field String specifying the column of \code{x} containing the read/UMI coverage.
#' @param downsample Logical scalar indicating whether downsampling should be performed.
#' @param down.ncells Integer scalar indicating the number of cells to downsample each group to.
#' Defaults to the smallest number of sequence-containing cells across all levels in \code{group}.
#' @param ... For the generic, further arguments to pass to individual methods.
#'
#' For the \code{CompressedSplitDataFrameList} method, further arguments to pass to the ANY method.
#'
#' @return An \linkS4class{IntegerList} containing one integer vector per level of \code{group} 
#' (or all cells, if \code{group=NULL}).
#' Each entry of the vector corresponds to a clonotype and contains the number of cells with that clonotype.
#' Each vector is also sorted in decreasing order.
#'
#' @author Aaron Lun
#'
#' @details
#' The aim of this function is to quantify clonal expansion based on the number of cells of a particular clonotype.
#' Clonal expansion is of interest as it serves as a proxy for the strength of the immune response to antigens;
#' we can then compare the degree of expansion between experimental conditions or cell states to gain some biological insights.
#'
#' Greater expansion manifests in the form of (i) more clonotypes with multiple cells
#' and (ii) clonotypes with a greater number of cells.
#' The exact effect probably depends on the nature of the antigen, e.g., the number of exposed epitopes,
#' that determine whether the expansion is spread across a larger number of clones.
#' See the \code{\link{summarizeClonalExpansion}} function for more details.
#' 
#' When \code{cov.field} is specified, only the most high-abundance sequence is used from each cell.
#' In contrast, setting \code{cov.field=NULL} will count each sequence separately,
#' such that one cell may contribute multiple times.
#' It is probably safest to set this to some non-\code{NULL} value to avoid complications from dependencies between counts,
#' though any problems are also probably minor.
#'
#' Cells without any clonotype are completely ignored within this function,
#' as they do not contribute to any of the clonotype counts.
#'
#' @section Normalization for cell number:
#' One difficulty with quantification is that the average cells per clonotype and number of multiple-cell clonotypes
#' is not a linear function with respect to the number of cells.
#' Increasing the number of cells may result in more new clonotypes or more cells assigned to previously observed clonotypes, 
#' depending on the (unknown) clonotype composition of the population.
#' This complicates comparisons between groups that contain different numbers of cells,
#' e.g., diversity metrics in \code{\link{summarizeClonotypeCounts}} cannot be directly compared between groups.
#'
#' We solve this problem by simply downsampling so that all levels of \code{group} have the same number of cells.
#' This eliminates uninteresting technical differences in, e.g., cell capture rates when comparing between groups,
#' without making any assumptions about the clonotype composition of the vast majority of unobserved cells.
#' Of course, this also eliminates the biological effect of the increase in the number of cells upon expansion,
#' but any such expansion should hopefully still be detectable via changes in clonotype composition among the remaining cells.
#' 
#' As discussed in \code{\link{countCellsPerGeneCombo}}, we do not adjust for differences in sequencing depth across groups.
#' Our assumption is that any change in coverage manifests as a scaling of the probability of detecting a clonotype
#' where the magnitude of scaling is the same across all clonotypes.
#' If so, differences in coverage translate to differences in the number of cells with a clonotype,
#' allowing us to use the same downsampling solution above.
#'
#' @examples
#' df <- data.frame(
#'     cell.id=sample(LETTERS, 30, replace=TRUE),
#'     clonotype=sample(paste0("clonotype_", 1:5), 30, replace=TRUE),
#'     umi=pmax(1, rpois(30, 5))
#' )
#'
#' y <- splitToCells(df, field="cell.id")
#' out <- countCellsPerClonotype(y, "clonotype", cov.field="umi")
#' out
#'
#' out2 <- countCellsPerClonotype(y, "clonotype", cov.field="umi",
#'    group=sample(3, length(y), replace=TRUE))
#' out2
#' 
#' @export
setGeneric("countCellsPerClonotype", function(x, ...) standardGeneric("countCellsPerClonotype"))

#' @importFrom IRanges IntegerList
.countCellsPerClonotype <- function(x, clone.field, group=NULL, cov.field=NULL, downsample=FALSE, down.ncells=NULL) {
    ids <- x[,clone.field]

    # Avoid generating zeros from levels absent in a particular group.
    if (is.factor(ids)) { ids <- as.character(ids) } 

    if (is.null(group)) {
        out <- list(table(ids))
    } else {
        out <- split(ids, group)
        out <- lapply(out, table)
    }

    for (i in seq_along(out)) {
        current <- out[[i]]
        current <- sort(current, decreasing=TRUE)
        out[[i]] <- as.integer(current)
        names(out[[i]]) <- names(current)
    }

    if (downsample) {
        out <- .downsample_list(out, down.ncells)
    }

    IntegerList(out)
}

#' @importFrom S4Vectors Rle runLength runValue
.downsample_counts <- function(x, n) {
    stuff <- Rle(seq_along(x), x)
    down <- stuff[sort(sample(length(stuff), n))]

    output <- x
    output[] <- 0L
    output[runValue(down)] <- runLength(down)

    output
}

.downsample_list <- function(counts, down.ncells) {
    if (is.null(down.ncells)) {
        down.ncells <- min(vapply(counts, sum, 0L))
    }

    for (i in seq_along(counts)) {
        current <- counts[[i]]
        current <- .downsample_counts(current, n=down.ncells)

        # Important to strip out zeroes, as we don't report
        # zero-frequency clonotypes anywhere else.
        current <- current[current!=0L]
        counts[[i]] <- current
    }

    counts
}


#' @export
#' @rdname countCellsPerClonotype
setMethod("countCellsPerClonotype", "ANY", .countCellsPerClonotype)

#' @export
#' @rdname countCellsPerClonotype
#' @importClassesFrom IRanges CompressedSplitDataFrameList
setMethod("countCellsPerClonotype", "CompressedSplitDataFrameList", function(x, clone.field, cov.field, group=NULL, ...) {
    if (!is.null(cov.field)) {
        cov <- x[,cov.field]
        x <- x[cov==max(cov)]
    }

    if (!is.null(group)) {
        group <- group[rep(seq_along(x), lengths(x))]
    }

    .countCellsPerClonotype(unlist(x), clone.field=clone.field, group=group, ...)
})
