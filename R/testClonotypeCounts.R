#' Test for differences in clonotype diversity
#'
#' Test for significant differences in the diversity of clonotypes between groups.
#'
#' @inheritParams summarizeClonotypeCounts
#' @param iterations Positive integer scalar indicating the number of permutation iterations to use for testing.
#' @param adj.method String specifying the multiple testing correction method to use across pairwise comparisons.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how parallelization should be performed.
#'
#' @return A \linkS4class{List} of numeric matrices containing p-values for pairwise comparisons of diversity between groups.
#' Each matrix is lower-triangular as the tests do not consider directionality.
#'
#' @details
#' This function computes permutation p-values to test for significant differences in the diversity values of different groups,
#' as computed using \code{\link{summarizeClonotypeCounts}}.
#' The aim is to help to whether one group is significantly more or less diverse,
#' providing evidence for differences in the rate of clonal expansion between clusters or conditions.
#'
#' Under the null hypothesis, two groups are derived from a pool of cells with the same clonotype composition (see below).
#' We sample without replacement to obtain two permuted groups that match the size of the original groups,
#' recompute the diversity indices for each permuted group and calculate the absolute difference of the indices.
#' Our permutation p-value is computed by comparing the observed absolute difference with the null distribution,
#' using the Phipson and Smyth (2010) approach to avoid p-values of zero.
#'
#' We repeat this process for each diversity index, e.g., Gini index, Hill numbers.
#' This yields a matrix of p-values per index where each row and column represents a group.
#' Within each index, we apply a multiple testing correction over all pairwise comparisons between groups.
#' By default, we use the Holm-Bonferroni correction to control the FWER across all comparisons.
#'
#' @section Comments on testing:
#' The null distribution depends on the composition of the common pool of cells.
#' This is not known so we approximate it by from the composition of the two groups being compared.
#' We rank all clonotypes within each group and sum the frequencies of clonotypes with the same rank between groups,
#' yielding a common population from which sampling without replacement is performed.
#'
#' We do not perform the more obvious strategy of creating a pool of clonotypes from both groups,
#' e.g., by literally concatenating the respective integer vectors from \code{counts}.
#' This strategy effectively doubles the number of available clonotypes used to compute the diversity indices,
#' making it difficult to justify using the null distribution to compute a p-value upon comparison to the observed difference.
#' 
#' Again, it is a good idea to downsample to ensure that all groups are of the same size.
#' Otherwise, the permutation test will not be symmetric;
#' it will only ever be significant if the larger group has the larger index.
#'
#' @author Aaron Lun
#'
#' @examples
#' df <- data.frame(
#'     cell.id=sample(LETTERS, 30, replace=TRUE),
#'     clonotype=sample(paste0("clonotype_", 1:5), 30, replace=TRUE),
#'     umi=pmax(1, rpois(30, 2))
#' )
#' 
#' y <- splitDataFrameByCell(df, field="cell.id")
#' out <- countCellsPerClonotype(y, "clonotype", cov.field="umi",
#'    group=sample(3, length(y), replace=TRUE))
#'
#' test.out <- testClonotypeCountsPairwise(out)
#' test.out$gini
#'
#' @seealso
#' \code{\link{summarizeClonotypeCounts}}, to compute diversity indices.
#'
#' @export
#' @importFrom stats p.adjust
#' @importFrom BiocParallel SerialParam bpmapply
#' @importFrom S4Vectors List
#' @importFrom dqrng generateSeedVectors
testClonotypeCountsPairwise <- function(counts, 
    use.gini=TRUE, use.hill=0:2,
    downsample=TRUE, down.ncells=NULL, 
    iterations=2000, adj.method="holm", BPPARAM=SerialParam()) 
{
    if (downsample) {
        counts <- .downsample_list(counts, down.ncells)
    }

    all.pairs <- combn(length(counts), 2)
    all.seeds <- generateSeedVectors(ncol(all.pairs))
    all.streams <- seq_len(ncol(all.pairs))

    res <- bpmapply(x=all.pairs[1,], y=all.pairs[2,], rand.seed=all.seeds, rand.stream=all.streams,
        FUN=function(x, y, counts, ...) .generate_shuffled_diversity_p(counts[[x]], counts[[y]], ...),
        MoreArgs=list(counts=counts, use.gini=use.gini, use.hill=use.hill, iterations=iterations),
        SIMPLIFY=FALSE, BPPARAM=BPPARAM)

    output <- list()
    if (length(counts) <= 1L) {
        stop("'counts' should contain at least two groups")
    }
    all.stat.names <- names(res[[1]])

    for (stat in all.stat.names) {
        current <- matrix(NA_real_, length(counts), length(counts))
        dimnames(current) <- list(names(counts), names(counts))

        for (i in seq_len(ncol(all.pairs))) {
            y <- all.pairs[1,i]
            x <- all.pairs[2,i]
            current[x,y] <- res[[i]][stat]
        }

        current[] <- p.adjust(current, method=adj.method)
        output[[stat]] <- current
    }

    List(output)
}

#' @importFrom S4Vectors Rle runValue runLength
.generate_shuffled_diversity_p <- function(x, y, iterations, use.gini, use.hill, rand.stream, rand.seed) {
    if (use.gini) {
        ref.gini <- abs(.compute_gini(x) - .compute_gini(y))
        out.gini <- 0L
    }
    if (length(use.hill)) {
        ref.hill <- abs(calcDiversity(x, use.hill) - calcDiversity(x, use.hill))
        out.hill <- integer(length(ref.hill))
    }

    permuted <- permute_diversity(x, y, iterations=iterations, 
        use_gini=use.gini, use_hill=use.hill, 
        seed=rand.seed, stream=rand.stream)
    
    # Using Phipson & Smyth's approach:
    output <- numeric(0)

    if (use.gini) {
        out.gini <- sum(abs(permuted[[1]]) >= ref.gini)
        output["gini"] <- (out.gini + 1)/(iterations+1)
    }
    if (length(use.hill)) {
        out.hill <- rowSums(abs(permuted[[2]]) >= out.hill)
        output[sprintf("hill%i", use.hill)] <- (out.hill + 1)/(iterations+1)
    }

    output
}
