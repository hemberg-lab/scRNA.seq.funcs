#' @export
calc_fpkm <- function(counts, gene2length) {
        order = match(rownames(counts),gene2length[,1]);
        norm = counts[!is.na(order),]
        order = order[!is.na(order)]
        weighted_sum_lengths = colSums(norm*as.numeric(gene2length[order,2]))
        fpkms = t(t(norm)/weighted_sum_lengths)*10^9
        return(fpkms)
}

#' @export
calc_cpm <- function(expr_mat, spikes = NULL) {
    norm_factor <- colSums(expr_mat[-spikes,])
    return(t(t(expr_mat)/norm_factor))*10^6
}

#' @export
calc_sf <- function(expr_mat, spikes = NULL) {
    geomeans <- exp(rowMeans(log(expr_mat[-spikes,])))
    SF<-function(cnts){median((cnts/geomeans)[(is.finite(geomeans) & geomeans > 0)])}
    norm_factor <- apply(expr_mat[-spikes,], 2, SF)
    return(t(t(expr_mat)/norm_factor))
}

#' @export
calc_uq <- function(expr_mat, spikes = NULL) {
    UQ <- function(x){quantile(x[x > 0], 0.75)}
    uq <- unlist(apply(expr_mat[-spikes, ], 2, UQ))
    norm_factor <- uq/median(uq)
    return(t(t(expr_mat)/norm_factor))
}

#' @export
calc_cell_RLE <- function(expr_mat, spikes = NULL) {
    RLE_gene <- function(x) {
        if (median(unlist(x)) > 0) {
            log( (x + 1) / ( median(unlist(x)) + 1 ) ) / log(2)
        } else {
            rep(NA, times = length(x))
        }
    }
    if(!is.null(spikes)) {
        RLE_matrix <- t(apply(expr_mat[-spikes, ], 1, RLE_gene))
    } else {
        RLE_matrix <- t(apply(expr_mat, 1, RLE_gene))
    }
    cell_RLE <- apply(RLE_matrix, 2, median, na.rm = T)
    return(cell_RLE)
}

#' @export
Down_Sample_Matrix <- function(expr_mat) {
    min_lib_size <- min(colSums(expr_mat))
    down_sample <- function(x) {
        prob <- min_lib_size/sum(x)
        return(
            unlist(
                lapply(
                    x,
                    function(y) {
                        rbinom(1, y, prob)
                    }
                )
            )
        )
    }
    down_sampled_mat <- apply(expr_mat, 2, down_sample)
    return(down_sampled_mat)
}
