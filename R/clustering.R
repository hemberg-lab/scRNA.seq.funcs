#' @export
gene_filter <- function(data, fraction) {
    filter.params <- filter_params(data, fraction)
    min.cells <- filter.params$min.cells
    max.cells <- filter.params$max.cells
    min.reads <- filter.params$min.reads
    d <- data[rowSums(data > min.reads) >= min.cells &
                  rowSums(data > 0) <= dim(data)[2] - max.cells, ]
    d <- unique(d)
    return(d)
}

filter_params <- function(dataset, fraction) {
    n.cells <- dim(dataset)[2]

    min.cells <- ceiling(fraction*n.cells)
    max.cells <- ceiling(fraction*n.cells)
    min.reads <- 2

    return(list(min.cells = min.cells, max.cells = max.cells,
                min.reads = min.reads))
}

#' @export
merge_pcaReduce_results <- function(dat, k) {
    res <- data.frame()
    for(i in 1:length(dat)) {
        tmp <- dat[[i]][,32 - k]
        res <- rbind(res, tmp)
    }
    colnames(res) <- colnames(dat)
    return(res)
}

#' @export
tsne_mult <- function(dataset, ks, n) {
    res <- list()
    for(i in 1:n) {
        tsne_out <- Rtsne::Rtsne(t(dataset), perplexity = 10) # Run TSNE
        tmp <- data.frame()
        for(k in ks) {
            t <- kmeans(tsne_out$Y, k, iter.max = 1e9, nstart = 1000)$clust
            tmp <- rbind(tmp, t)
        }
        colnames(tmp) <- colnames(dataset)
        rownames(tmp) <- ks
        res[[i]] <- tmp
    }
    return(res)
}


#########################################################
# This program is part of the SNN-Cliq method           #
# Contact Chen Xu at UNC-Charlotte for more information.#
#########################################################
#----- example of use------#
#data<-read.table(infile, header=TRUE, sep="\t", row.names=1);
#data<-log2(data+1)
#source('SNN.R')
#SNN(data, edge_file, k=3, distance='euclidean')
#--------------------------#
#' @export
SNN<-function(data, outfile, k, distance){

    if(missing(data)){
        stop(paste("Input data missing.",help,sep="\n"))
    }
    if(missing(outfile)){
        stop(paste("Output file name missing.",help,sep="\n"))
    }
    if(missing(k)){
        k=3
    }
    if(missing(distance)){
        distance<-"euclidean"  # other distance options refer to dist() in R
    }
    m<-as.data.frame(data)
    numSpl<-dim(data)[1]
    m<-dist(data, distance, diag=TRUE, upper=TRUE)
    x<-as.matrix(m)
    IDX<-t(apply(x,1,order)[1:k,]) # knn list

    edges<-list()              # SNN graph
    for (i in 1:numSpl){
        j<-i
        while (j<numSpl){
            j<-j+1
            shared<-intersect(IDX[i,], IDX[j,])
            if(length(shared)>0){
                s<-k-0.5*(match(shared, IDX[i,])+match(shared, IDX[j,]))
                strength<-max(s)
                if (strength>0)
                    edges<-rbind(edges, c(i,j,strength))
            }
        }
    }
    write.table(edges, outfile, quote=FALSE, sep='\t',col.names=FALSE,row.names=FALSE)
}

#' @export
z.transform.helper <- function(x) {
    x <- as.numeric(x)
    x.mu <- mean(x)
    x.sd <- sd(x)
    if (x.sd == 0) {
        x <- rep(0, length(x))
    } else {
        x <- (x-x.mu)/x.sd
    }
    return(x)
}
