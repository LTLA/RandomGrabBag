#' Count shared clonotypes across groups
#'
#' Count the number of clonotypes that are shared across groups, usually different cell types.
#'
#' @inheritParams countCellsPerClonotype
#' @param metric String specifying the type of sharing metric to return.
#' @param collapse.cells Logical scalar indicating whether each clonotype should be counted only once.
#'
#' @return 
#' A numeric matrix where the off-diagonal entries contain metrics of clonotype sharing between each pair of groups.
#'
#' @details
#' This function quantifies the sharing of clonotypes are shared across different groups.
#' The most obvious application is to identify shared clonotypes across different cell types (as represented by cluster identity),
#' allowing us to infer that those cell types shared a common ancestor.
#' Examples include clonotypes that are shared across various B or T cell states (e.g., activation, memory),
#' indicating that there is an active transition between states. 
#'
#' If \code{metric="none"}, we return the number of clonotypes that are shared between each pair of groups.
#' If \code{collapse.cells=FALSE}, we instead return the total number of cells across both groups that exhibit a shared clonotype.
#'
#' If \code{metric="jaccard"}, the Jaccard index is computed between every pair of groups.
#' On a practical level, this adjusts for differences in the size of the groups so that large groups do not dominate the output.
#' On a theoretical level, we interpret this value by considering the progenitor population that gives rise to the two groups;
#' the Jaccard index represents the proportion of cells in the progenitor population that can develop into either group.
#' If \code{collapse.cells=FALSE}, a weighted version of the Jaccard index is used, i.e., the Ruzicka similarity.
#' This involves summing the minimum and maximum frequencies of each clonotype in the numerator and denominator, respectively.
#' 
#' If \code{metric="maximum"}, we compute the larger of the proportions of shared clonotypes across the two groups.
#' For example, if 10 clonotypes are shared and one group has 20 clonotypes and the other group has 40 clonotypes,
#' the output value will be 0.5.
#' This is designed to detect high shared proportions in smaller groups;
#' for example, if a rare subpopulation is derived from the same progenitors as a much larger population,
#' the Jaccard index would be small even if the sharing in the small subpopulation is very high.
#' If \code{collapse.cells=FALSE}, we instead compute the proportion of cells in each group that exhibit a shared clonotype.
#' 
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{countCellsPerClonotype}}, to compare between groups.
#'
#' @examples
#' df <- data.frame(
#'     cell.id=sample(LETTERS, 30, replace=TRUE),
#'     clonotype=sample(paste0("clonotype_", 1:5), 30, replace=TRUE),
#'     umi=pmax(1, rpois(30, 5))
#' )
#'
#' y <- splitDataFrameByCell(df, field="cell.id")
#' out <- countSharedClonotypes(y, "clonotype", cov.field="umi",
#'    group=sample(3, length(y), replace=TRUE))
#' out
#'
#' @export
setGeneric("countSharedClonotypes", function(x, ...) standardGeneric("countSharedClonotypes"))

.countSharedClonotypes <- function(x, clone.field, group, metric=c("none", "jaccard", "maximum"), collapse.cells=FALSE) {
    val <- x[,clone.field]
    by.group <- split(val, group)
    if (collapse.cells) {
        by.group <- lapply(by.group, unique)
    }

    metric <- match.arg(metric)
    output <- matrix(NA_real_, length(by.group), length(by.group))
    dimnames(output) <- list(names(by.group), names(by.group))

    for (i in seq_along(by.group)) {
        g1 <- by.group[[i]]
        for (j in seq_len(i-1L)) {
            g2 <- by.group[[j]]

            if (metric=="maximum") {
                common <- intersect(g1, g2)
                current <- max(mean(g1 %in% common), mean(g2 %in% common))
            } else if (metric=="jaccard") {
                universe <- union(g1, g2)
                freq1 <- as.integer(table(factor(g1, levels=universe)))
                freq2 <- as.integer(table(factor(g2, levels=universe)))
                current <- sum(pmin(freq1, freq2))/sum(pmax(freq1, freq2))
            } else {
                common <- intersect(g1, g2)
                current <- sum(g1 %in% common) + sum(g2 %in% common)
            }

            output[i,j] <- current
        }
    }

    output
}

#' @export
#' @rdname countSharedClonotypes
setMethod("countSharedClonotypes", "ANY", .countSharedClonotypes)

#' @export
#' @rdname countSharedClonotypes
#' @importClassesFrom IRanges CompressedSplitDataFrameList
setMethod("countSharedClonotypes", "CompressedSplitDataFrameList", function(x, clone.field, cov.field, group, ...) {
    if (!is.null(cov.field)) {
        cov <- x[,cov.field]
        x <- x[cov==max(cov)]
    }

    if (!is.null(group)) {
        group <- group[rep(seq_along(x), lengths(x))]
    }

    .countSharedClonotypes(unlist(x), clone.field=clone.field, group=group, ...)
})
