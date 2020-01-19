#' Create a SplitDataFrameList
#'
#' Create a SplitDataFrameList to store repertoire data.
#' Each cell corresponds to one DataFrame, while each row of that DataFrame corresponds to a sequence detected in the cell.
#' 
#' @param df A data frame of repertoire data where each row corresponds to a detected sequence.
#' @param field A string specifying the field of \code{df} that contains the identity of the cell.
#' @param ids A character vector of length equal to \code{nrow(df)} specifying the identity of the cell for each sequence.
#' @param universe Character vector of all possible cell identities.
#'
#' @details
#' Setting \code{universe} is useful to ensure that the output object is of the same length as,
#' say, the number of columns in a \linkS4class{SingleCellExperiment} object containing expression data.
#'
#' @return
#' A \linkS4class{SplitDataFrameList} containing a per-cell perspective of repertoire data.
#'
#' @author Aaron Lun
#'
#' @examples
#' df <- data.frame(
#'     cell.id=sample(LETTERS, 30, replace=TRUE),
#'     v_gene=sample(c("TRAV1", "TRAV2", "TRAV3"), 30, replace=TRUE),
#'     j_gene=sample(c("TRAJ4", "TRAJ5", "TRAV6"), 30, replace=TRUE)
#' )
#'
#' splitToCells(df, field="cell.id")
#'
#' splitToCells(df, field="cell.id", universe=c(LETTERS, "AA"))
#'     
#' @export
#' @importFrom S4Vectors split DataFrame
splitToCells <- function(df, field, ids=df[[field]], universe=sort(unique(ids))) {
    split(DataFrame(df), factor(ids, universe))
}


#' Utility counting functions 
#'
#' Pretty much as it says.
#' These are utility functions for counting diagnostics in the repertoire chapter of the OSCA book.
#' Should really find a better home somewhere else.
#'
#' @param x A \linkS4class{SplitDataFrameList} where each \linkS4class{DataFrame} is a cell and each row is a sequence.
#'
#' @return
#' For \code{countSequences}, an integer scalar 
#'
#' @author Aaron Lun
#' 
NULL

#' @export
countSequencesPerCell <- function(x, filter.field=NULL, filter.value=NULL) {
    if (length(field)) {
        keep <- TRUE
        for (i in seq_along(filter.field)) {
            keep <- keep & x[,filter.field[[i]]]==filter.value[[i]]
        }
        sum(keep)
    } else {
        lengths(x) 
    }
}

#' @export
topCoveragePropPerCell <- function(x, cov.field) {
    counts <- x[,cov.field]
    totals <- sum(counts)
    max <- max(counts)

    out <- max/totals
    out[totals==0] <- NA_real_
    out
}

#' @export
umiToReadPerCell <- function(x, read.field, umi.field) {
    sum(x[,umi.field])/sum(x[,read.field])
}

#' @export
#' @importFrom S4Vectors selfmatch
#' @importFrom IRanges IntegerList
countGeneCombo <- function(x, gene.field, weight=FALSE, cov.field=NULL) {
    if (!is.null(cov.field)) {
        cov <- x[,cov.field]
        x <- x[cov==max(cov)]
    } else if (weight) {
        w <- rep(1/lengths(x), lengths(x))        
    }

    y <- unlist(x)[,gene.field]
    ids <- selfmatch(y)
    keep <- !duplicated(ids)
    y <- y[keep,]

    if (weight) { 
        y$count <- sum(NumericList(split(w, ids)))
    } else {
        y$count <- as.integer(table(ids))
    }

    y[order(y$count, decreasing=TRUE),]
}
