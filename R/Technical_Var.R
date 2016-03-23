#' @export
plot_expression_heatmap <- function (genes, data, cell_labels=NA, gene_labels=NA, key_genes=genes, key_cells=NA) {
	require("RColorBrewer")
	require("gplots")
	# Check input
	cell_labels = as.character(cell_labels);
	# Get gene rows
	if(!is.numeric(genes)) {
		genes = match(genes, rownames(data));
		nomatch = sum(is.na(genes));
		if (nomatch > 0) {warning(paste(nomatch, " genes could not be matched to data, they will not be included in the heatmap."));}
		genes = genes[!is.na(genes)];
	}
	if (length(genes) < 1) {warning("No genes for heatmap.");return();}
	# Plot heatmap of expression
	heatcolours <- rev(brewer.pal(11,"RdBu"))
	col_breaks = c(-100,seq(-2,2,length=10),100)
	heat_data = as.matrix(data[genes,])
	heat_data = log(heat_data+1)/log(2);
	ColColors = rep("white", times=length(heat_data[1,]))
	RowColors = rep("white", times=length(heat_data[,1]))
	# remove row & column labels
	rownames(heat_data) = rep("", length(heat_data[,1]));
	if (!is.na(key_genes)) {
		rownames(heat_data)[rownames(data[genes,]) %in% key_genes] = rownames(data[genes,])[rownames(data[genes,]) %in% key_genes];
	}
	colnames(heat_data) = rep("", length(heat_data[1,]));
	if (!is.na(key_cells)) {
		colnames(heat_data)[colnames(data[genes,]) %in% key_cells] = colnames(data[genes,])[colnames(data[genes,]) %in% key_cells];
	}
	if (!is.na(cell_labels[1])) {
		colours = as.factor(cell_labels)
		palette = brewer.pal(max(3,length(unique(cell_labels))), "Set3");
		ColColors = palette[colours];
		mylegend<- list(names = unique(cell_labels), fill = unique(ColColors));
	}
	if (!is.na(gene_labels[1])) {
		# lowest factor level = grey (so 0-1 is striking)
		if (!is.numeric(gene_labels)) {
			colours = as.factor(gene_labels)
		} else {
			colours = gene_labels
		}
		palette = c("grey75",brewer.pal(max(3,length(unique(gene_labels))), "Set1"));
		RowColors = palette[colours];
	}
	# Custom Shit
	lwid=c(1,0.2,4)
	lhei=c(1,0.2,4)
	lmat=rbind(c(6,0,5),c(0,0,2),c(4,1,3))


	heatmap_output = heatmap.2(heat_data, ColSideColors = ColColors, RowSideColors = RowColors, col=heatcolours, breaks=col_breaks, scale="row",symbreaks=T, trace="none", dendrogram="column", key=FALSE, Rowv=TRUE, Colv=TRUE,lwid=lwid, lhei=lhei,lmat=lmat, hclustfun=function(x){hclust(x,method="ward.D2")})
	# Custom key
	par(fig = c(0, 1/(5.2),4/(5.2), 1), mar=c(4,1,1,1), new=TRUE)
	scale01 <- function(x, low = min(x), high = max(x)) {
        	x <- (x - low)/(high - low)
        	x
    	}
	par(mar=c(5,1,1,1))
	par(cex=0.75)
	par(mgp=c(2,1,0))
	key_breaks = seq(-2,2,length=10)
	key_col = heatcolours[2:(length(heatcolours)-1)]
	z = seq(min(key_breaks),max(key_breaks), by=min(diff(key_breaks)/4))
	image(z=matrix(z,ncol=1),col=key_col,breaks=key_breaks,xaxt="n",yaxt="n")
	par(usr = c(0, 1, 0, 1))
	lv <- pretty(key_breaks)
        xv <- scale01(as.numeric(lv), min(key_breaks),max(key_breaks))
        xargs <- list(at = xv, labels = lv)
	xargs$side <- 1
	do.call(axis, xargs)
	mtext(side = 1, "Expression Z-Score", line = par("mgp")[1], padj = 0.5,
                cex = par("cex") * par("cex.lab"))

	# Legend
	par(fig = c(0/5.2, 1/(5.2),0/(5.2), 4/5.2), mar=c(0,0,0,0), new=TRUE)
	par(mar=c(0,0,0,0))
	if (!is.na(cell_labels[1])) {
		legend("left", mylegend$names, pt.bg = mylegend$fill,bg="white",col="black", pch=22, pt.cex=2.5, cex=1.25, bty="n",y.intersp = 2);
	}
	invisible(heatmap_output);
}

#' @export
#from : http://www.nature.com/nmeth/journal/v10/n11/full/nmeth.2645.html#supplementary-information
Brennecke_getVariableGenes <- function(expr_mat, spikes=NA, suppress.plot=FALSE, fdr=0.1, minBiolDisp=0.5) {
        require(statmod)

        rowVars <- function(x) { unlist(apply(x,1,var))}

        colGenes = "black"
        colSp = "grey35"


        fullCountTable <- expr_mat;

        if (is.character(spikes)) {
                sp = rownames(fullCountTable) %in% spikes;
                countsSp <- fullCountTable[sp,];
                countsGenes <- fullCountTable[!sp,];
        } else if (is.numeric(spikes)) {
                countsSp <- fullCountTable[spikes,];
                countsGenes <- fullCountTable[-spikes,];
        } else {
                countsSp = fullCountTable;
                countsGenes = fullCountTable;
        }

        meansSp = rowMeans(countsSp)
        varsSp = rowVars(countsSp)
        cv2Sp = varsSp/meansSp^2
        meansGenes = rowMeans(countsGenes)
        varsGenes = rowVars(countsGenes)
        cv2Genes = varsGenes/meansGenes^2
        # Fit Model
        minMeanForFit <- unname( quantile( meansSp[ which( cv2Sp > 0.3 ) ], 0.80))
        useForFit <- meansSp >= minMeanForFit
#        if (sum(useForFit) < 50) {
#                warning("Too few spike-ins exceed minMeanForFit, recomputing using all genes.")
#                meansAll = c(meansGenes, meansSp)
#                cv2All = c(cv2Genes,cv2Sp)
#                minMeanForFit <- unname( quantile( meansAll[ which( cv2All > 0.3 ) ], 0.80))
#                useForFit <- meansSp >= minMeanForFit
#        }
        if (sum(useForFit) < 30) {warning(paste("Only", sum(useForFit), "spike-ins to be used in fitting, may result in poor fit."))}
        fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansSp[useForFit] ), cv2Sp[useForFit] )
        a0 <- unname( fit$coefficients["a0"] )
        a1 <- unname( fit$coefficients["a1tilde"])

        # Test
        psia1theta <- a1
        minBiolDisp <- minBiolDisp^2
        m = ncol(countsSp);
        cv2th <- a0 + minBiolDisp + a0 * minBiolDisp
        testDenom <- (meansGenes*psia1theta + meansGenes^2*cv2th)/(1+cv2th/m)
        p <- 1-pchisq(varsGenes * (m-1)/testDenom,m-1)
        padj <- p.adjust(p,"BH")
        sig <- padj < fdr
        sig[is.na(sig)] <- FALSE
        if (!suppress.plot) {
                plot( meansGenes,cv2Genes, xaxt="n", yaxt="n", log="xy",
                        xlab = "average normalized read count",
                        ylab = "squared coefficient of variation (CV^2)", col="white")
                axis( 1, 10^(-2:5), c( "0.01", "0.1", "1", "10", "100", "1000",
                        expression(10^4), expression(10^5) ) )
                axis( 2, 10^(-2:3), c( "0.01", "0.1", "1", "10", "100","1000" ), las=2 )
                abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )
                # Plot the genes, use a different color if they are highly variable
                points( meansGenes, cv2Genes, pch=20, cex=.2,
                        col = ifelse( padj < .1, "#C0007090", colGenes ) )
		points( meansSp, cv2Sp, pch=20, cex=.5, col="blue1")
                # Add the technical noise fit
                xg <- 10^seq( -2, 6, length.out=1000 )
                lines( xg, (a1)/xg + a0, col="#FF000080", lwd=3 )
                # Add a curve showing the expectation for the chosen biological CV^2 thershold
                lines( xg, psia1theta/xg + a0 + minBiolDisp, lty="dashed", col="#C0007090", lwd=3)
        }
        return(names(meansGenes)[sig])
}

