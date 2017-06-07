globalVariables(names=c("conversor","conversor.mouse","conversor","conversor.mouse",
"org.Hs.egSYMBOL","org.Mm.egSYMBOL","microCosm_v5_18","genes_human_h37",
"mirnas_human_17_h37","grup2","grup1","version.mirnas"),package="miRComb")




checkmiRNAs <- function ( v.miRNAs, to.dataframe=FALSE ) {

	if (!exists("versions.mirnas")) { data(versions.mirnas) }
	versions.mirnas<-versions.mirnas[,-1]

	perc <- function (x,target) {
		return(length(which((target %in% x) == TRUE)) / length(target) *100)
	}

	a<-apply(versions.mirnas,2,perc,v.miRNAs)

	dat<-data.frame(x=names(a),y=a)

	dat$x<-factor(gsub("miRBase_","",dat$x))
	dat$x<-relevel(dat$x,"9.2")
	dat$x<-relevel(dat$x,"9.1")
	dat$x<-relevel(dat$x,"9.0")
	dat$x<-relevel(dat$x,"8.2")
	dat$x<-relevel(dat$x,"8.1")
	dat$x<-relevel(dat$x,"8.0")
	dat$x<-relevel(dat$x,"7.1")
	dat$x<-relevel(dat$x,"7.0")
	dat$x<-relevel(dat$x,"6.0")

	if (to.dataframe==TRUE) {
		return(dat)
	} else {

		qplot(x=dat$x, y=dat$y, geom="bar", stat="identity", fill=dat$x) +
  guides(fill=FALSE) + xlab("miRBase version") + ylab("Coincidence (%)")

	}

}


translatemiRNAs <- function ( x , from = NULL, to = "21", force.translation = FALSE, species = NULL) {

	if (!exists("versions.mirnas")) { data(versions.mirnas) }
	
	if (from=="unknown") {
		options<-checkmiRNAs(x, to.dataframe=TRUE)
		sel<-which(options[,2]==max(options[,2]))
		from<-as.character(options[sel[length(sel)],1])
		cat(paste("Your miRNAs have been assigned to miRBase version ",from,", with a coincidence of ", round(max(options[,2],1)),"%. Check checkmiRNAs function for more details.\n",sep=""))
	}
  
  from<-paste("miRBase_",as.character(from),sep="")
  to<-paste("miRBase_",as.character(to),sep="")
  
	sel<-which((versions.mirnas[,from] %in% x)==TRUE)
	subs<-versions.mirnas[sel,]
	rownames(subs)<-subs[,from]
	
	translated <- rep(NA, length(x))
	names(translated)<-x

	sel<-intersect(names(translated),rownames(subs))
	
	for (i in 1:length(sel)) {
		translated[which(x %in% sel[i])]<-as.character(subs[sel[i],to])
	}

	return(translated)
}



plotCordist <- function (obj, subset, type="cor", method.cor="pearson", method.dist="euclidean", hierarchical=FALSE, ...) {
		
	nf<-layout(mat=matrix(c(1),ncol=1,nrow=1))
	par(mar=c(5.1, 4.1, 4.1, 2.1))
#	plot.new()


myImagePlot <- function(x, ...){
     min <- min(x)
     max <- max(x)
     yLabels <- rownames(x)
     xLabels <- colnames(x)
     title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
       min <- Lst$zlim[1]
       max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
       yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
       xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
       title <- Lst$title
    }
  }
# check for null values
if( is.null(xLabels) ){
   xLabels <- c(1:ncol(x))
}
if( is.null(yLabels) ){
   yLabels <- c(1:nrow(x))
}

layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))

 # Red and green range from 0 to 1 while Blue ranges from 1 to 0
 ColorRamp <- rgb( seq(0,1,length=256),  # Red
                   seq(0,1,length=256),  # Green
                   seq(1,0,length=256))  # Blue
 ColorLevels <- seq(min, max, length=length(ColorRamp))

 # Reverse Y axis
 reverse <- nrow(x) : 1
 yLabels <- yLabels[reverse]
 x <- x[reverse,]

 # Data Map
 par(mar = c(3,5,2.5,2))
 image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
 ylab="", axes=FALSE, zlim=c(min,max))
 if( !is.null(title) ){
    title(main=title)
 }
axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7, las=2)
 axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
 cex.axis=0.7)

 # Color Scale
 par(mar = c(3,2.5,2.5,2))
 image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n")

 layout(1)
}


	if (subset=="miRNA") {
		if (type == "cor") {
			toplot<-cor(obj@dat.miRNA,method=method.cor)
			if (hierarchical==TRUE) {
				ord<-hclust(dist(toplot))$order
				toplot<-toplot[ord,ord]
			}
			myImagePlot(toplot)			
		}
		if (type == "dist") {
			toplot<-as.matrix(dist(t(obj@dat.miRNA),upper=TRUE,diag=TRUE,method=method.dist))
			if (hierarchical==TRUE) {
				ord<-hclust(dist(toplot))$order
				toplot<-toplot[ord,ord]
			}
			myImagePlot(toplot)
		}
	}


	if (subset=="mRNA") {
		if (type == "cor") {
			toplot<-cor(obj@dat.mRNA,method=method.cor)
			if (hierarchical==TRUE) {
				ord<-hclust(dist(toplot))$order
				toplot<-toplot[ord,ord]
			}
			myImagePlot(toplot)			
		}
		if (type == "dist") {
			toplot<-as.matrix(dist(t(obj@dat.mRNA),upper=TRUE,diag=TRUE,method=method.dist))
			if (hierarchical==TRUE) {
				ord<-hclust(dist(toplot))$order
				toplot<-toplot[ord,ord]
			}
			myImagePlot(toplot)
		}
	}
}


plotHclust <- function (obj, subset) {

	if (subset=="miRNA") {
		plot(hclust(dist(t(obj@dat.miRNA))),xlab="miRNA samples",col=c(1,2))

	}
	if (subset=="mRNA") {
		plot(hclust(dist(t(obj@dat.mRNA))),xlab="mRNA samples")

	}

}




plotCorrelation <- function (obj, miRNA, mRNA, type="cor",samples="all",col.color=1,i.legend=!is.na(col.color), pos.legend="topright", sample.names=FALSE, pos.sample.names=1, cex.main=1.35, alternative="two.sided", colors=c("turquoise","violet")) {

	if (samples=="all") {
		comm<-intersect(colnames(obj@dat.mRNA),colnames(obj@dat.miRNA))
	} else {comm<-samples}

	linmod<-lm(obj@dat.mRNA[mRNA,comm]~obj@dat.miRNA[miRNA,comm])
	
	if (type=="cor") {

	#	nf<-layout(mat=matrix(c(1),ncol=1,nrow=1))
		par(mar=c(5.1, 4.1, 4.1, 2.1))
	#	plot.new()


		if (!is.na(col.color)) {
			if (i.legend) {noms<-levels(as.factor(obj@pheno.miRNA[comm,col.color]))}
		#	col.color<-as.numeric(as.factor(obj@pheno.miRNA[comm,col.color]))

	if (length(levels(as.factor(obj@pheno.miRNA[comm,col.color])))==2) {
		col.color<-as.character(factor(as.numeric(as.factor(obj@pheno.miRNA[comm,col.color])),labels=colors))
	} else {
		col.color<-as.character(factor(as.numeric(as.factor(obj@pheno.miRNA[comm,col.color]))+1))
	}


		} else {col.color<-1;legend<-FALSE}
		corr<-round(cor(obj@dat.miRNA[miRNA,comm],obj@dat.mRNA[mRNA,comm]),3)
		pvall<-round(cor.test(obj@dat.miRNA[miRNA,comm],obj@dat.mRNA[mRNA,comm],alternative=alternative)$p.value,3)

		
		if (is.null(obj@net[paste(miRNA,mRNA,sep=":"),"adj.pval"])) {
			pvall<-round(cor.test(obj@dat.miRNA[miRNA,comm],obj@dat.mRNA[mRNA,comm],alternative=alternative)$p.value,3)
			main<-paste("Pearson Cor.: ",corr,"; p.val: ",pvall,sep="")
		} else {
			pvall<-round(obj@net[paste(miRNA,mRNA,sep=":"),"adj.pval"],4)
			main<-paste("Pearson Cor.: ",corr,"; adj.pval: ",pvall,sep="")
		}


		plot(obj@dat.miRNA[miRNA,comm], obj@dat.mRNA[mRNA,comm], pch=19, xlab=miRNA, ylab=mRNA, col=col.color, main=main, cex.main=cex.main, cex.lab=1.25)
		abline(linmod)
		if (i.legend) {
			legend(pos.legend,noms,col=levels(as.factor(col.color)),pch=19)
		}
		if (sample.names) {
			text(obj@dat.miRNA[miRNA,comm],obj@dat.mRNA[mRNA,comm],comm,pos=pos.sample.names)
		}
		
	}		
	if (type=="residuals") {
		nf<-layout(mat=matrix(c(1:4),ncol=2,nrow=2,byrow=TRUE))
		par(mar=c(5.1, 4.1, 4.1, 2.1))
		plot(linmod)
	}		
}

#get.info<-function (obj, miRNA=NULL, mRNA=NULL, list.pairs=NULL) {
#	if (is.null(list.pairs)) {
#		sel<-merge(miRNA,mRNA)
#		list.pairs<-paste(sel[,1],sel[,2],sep=":")
#	}
#	return(obj@net[list.pairs,])
#}


# ----- Define a function for plotting a matrix ----- #

# ----- END plot function ----- #


plot3d <- function (obj, subset, col.color=1, angle=45, colors=c("violet","turquoise"), lty=0, cex.points=1, ...) {

	nf<-layout(mat=matrix(c(1),ncol=1,nrow=1))
	par(mar=c(5.1, 4.1, 4.1, 2.1))
#	plot.new()

	#library(scatterplot3d)

	if (subset=="miRNA") {
		color<-as.character(factor(as.numeric(as.factor(obj@pheno.miRNA[,col.color]))+1,labels=colors))
		groups<-levels(as.factor(obj@pheno.miRNA[,col.color]))
		#remove with 0 variance
		sel<-which(apply(obj@dat.miRNA,1,sd)==0)
		if (length(sel)>0) {pca <- prcomp(t(obj@dat.miRNA[-sel,]), scale=T)}
		else {pca <- prcomp(t(obj@dat.miRNA), scale=T)}
		labels.dat<-colnames(obj@dat.miRNA)
	}

	if (subset=="mRNA") {
		color<-as.character(factor(as.numeric(as.factor(obj@pheno.mRNA[,col.color]))+1,labels=colors))
		groups<-levels(as.factor(obj@pheno.mRNA[,col.color]))
		#remove with 0 variance
		sel<-which(apply(obj@dat.mRNA,1,sd)==0)
		if (length(sel)>0) {pca <- prcomp(t(obj@dat.mRNA[-sel,]), scale=T)}
		else {pca <- prcomp(t(obj@dat.mRNA), scale=T)}
		labels.dat<-colnames(obj@dat.mRNA)
	}

	# create column indicating point color
	pca$pcolor<- color
	with(pca, {
  	  s3d <- scatterplot3d(pca$x[,1], pca$x[,2], pca$x[,3],        # x y and z axis
                  color=pcolor, pch=19,  
                  type="h", lty.hplot=lty,      
                  scale.y=.75,                
        #          main="3-D Scatterplot",
                  xlab=paste("Comp 1: ",round(pca$sdev[1]^2/sum(pca$sdev^2)*100,1),"%",sep=""),
                  ylab=paste("Comp 2: ",round(pca$sdev[2]^2/sum(pca$sdev^2)*100,1),"%",sep=""),
                  zlab=paste("Comp 3: ",round(pca$sdev[3]^2/sum(pca$sdev^2)*100,1),"%",sep=""),
		  angle=angle, cex.symbols=cex.points, ...)
   	  s3d.coords <- s3d$xyz.convert(pca$x[,1], pca$x[,2], pca$x[,3])
   #	  text(s3d.coords$x, s3d.coords$y,     # x and y coordinates
   #	       labels=labels.dat,       # text to plot
   #	       pos=4, cex=.5)                  # shrink text 50% and place to right of points)
# add the legend
	legend("topleft", inset=.05,      # location and inset
 	   bty="n", cex=1,              # suppress legend box, shrink text 50%
 	   title="Group",
 	   groups, col=c(levels(as.factor(color))),pch=19)
})
### ojuuuu, problemes amb el label!!!

}




selOutliers <- function (obj, subset, method="aq.plot", delete=FALSE, add.pheno=TRUE, n.dim=2) {

	if (subset=="miRNA") {
		data<-t(obj@dat.miRNA)
	}
	if (subset=="mRNA") {
		data<-t(obj@dat.mRNA)
	}

	if (method=="aq.plot") {
		#library("mvoutlier")
		pca <- prcomp(data, scale=TRUE, center=TRUE)
		outs <- aq.plot(pca$x[,1:n.dim],quan=1)
		outs.true <- which(outs$outliers==TRUE)
		num.out<-as.numeric(outs$outliers)
	}

	if (delete==TRUE) {
		if (subset=="miRNA") {
			if (length(outs.true > 0)) {
				obj@dat.miRNA<-as.matrix(t(data[-c(outs.true),]))
				obj@pheno.miRNA<-obj@pheno.miRNA[-c(outs.true),]
			} 

		}
		if (subset=="mRNA") {
			if (length(outs.true > 0)) {
				obj@dat.mRNA<-as.matrix(t(data[-c(outs.true),]))
				obj@pheno.mRNA<-obj@pheno.mRNA[-c(outs.true),]
			}

		}
		return(obj)
	} else {

		if (add.pheno==TRUE) { 
			if (subset=="miRNA") {
				obj@pheno.miRNA$is.outlier <- num.out
				plotPca(obj, "miRNA", col.color="is.outlier")
			}
			if (subset=="mRNA") {
				obj@pheno.mRNA$is.outlier <- num.out
				plotPca(obj, "mRNA", col.color="is.outlier")
			}
		return(obj)
		} else {
			return(rownames(data)[outs.true])
		}
	}


}

#br<-selOutliers(data.tcga,"mRNA", delete=FALSE ,n.dim=2)
#plotPca(br, "mRNA", col.color="is.outlier")



plotPca <- function (obj, subset, col.color=1, colors=c("turquoise","violet"), pos.leg="topleft", names=FALSE , ...) {

	#library(scatterplot3d)
	nf<-layout(mat=matrix(c(1),ncol=1,nrow=1))
	par(mar=c(5.1, 4.1, 4.1, 2.1))
#	plot.new()

	if (subset=="miRNA") {

		if (is.numeric(obj@pheno.miRNA[,col.color])==TRUE) {

		substna <- function (x) {
				x[which(is.na(x)==TRUE)] <- min(x,na.rm=TRUE)
				return(x)
			}			

			remap <- function(x) { (( x ) / max( abs(x) ) / 2)+0.5 }  # map x onto [0, 1]
			fun.col <- function(x) {rgb(colorRamp(c("violet","gray", "black"))(remap(substna(x))),
	                        maxColorValue = 255)
			}
			#nodeinfo2<-with(scorescomp, fun.col(as.numeric(scorescomp$scores)))		


			color <- with(obj@pheno.miRNA, fun.col(obj@pheno.miRNA[,col.color]))
			groups<-paste("min:", min(obj@pheno.miRNA[,col.color],na.rm=TRUE), ", max:",max(obj@pheno.miRNA[,col.color],na.rm=TRUE),sep="")


		} else {
			if (length(levels(as.factor(obj@pheno.miRNA[,col.color])))==2) {
				color<-as.character(factor(as.numeric(as.factor(obj@pheno.miRNA[,col.color])),labels=colors))
			} else {
				color<-as.character(factor(as.numeric(as.factor(obj@pheno.miRNA[,col.color]))+1))
			}
			groups<-levels(as.factor(obj@pheno.miRNA[,col.color]))
		}
		#remove with 0 variance
		sel<-which(apply(obj@dat.miRNA,1,sd)==0)
		if (length(sel)>0) {pca <- prcomp(t(obj@dat.miRNA[-sel,]), scale=T)}
		else {pca <- prcomp(t(obj@dat.miRNA), scale=T)}
		labels.dat<-colnames(obj@dat.miRNA)
	}

	if (subset=="mRNA") {
		if (length(levels(as.factor(obj@pheno.mRNA[,col.color])))==2) {
			color<-as.character(factor(as.numeric(as.factor(obj@pheno.mRNA[,col.color])),labels=colors))
		} else {
			color<-as.character(factor(as.numeric(as.factor(obj@pheno.mRNA[,col.color]))+1))
		}
		groups<-levels(as.factor(obj@pheno.mRNA[,col.color]))
		#remove with 0 variance
		sel<-which(apply(obj@dat.mRNA,1,sd)==0)
		if (length(sel)>0) {pca <- prcomp(t(obj@dat.mRNA[-sel,]), scale=T)}
		else {pca <- prcomp(t(obj@dat.mRNA), scale=T)}
		labels.dat<-colnames(obj@dat.mRNA)
	}

	# create column indicating point color
	pca$pcolor<- color
	plot(pca$x[,1], pca$x[,2],        # x y and z axis
                  col=pca$pcolor, pch=19,  
        #          main="3-D Scatterplot",
                  xlab=paste("Comp 1: ",round(pca$sdev[1]^2/sum(pca$sdev^2)*100,1),"%",sep=""),
                  ylab=paste("Comp 2: ",round(pca$sdev[2]^2/sum(pca$sdev^2)*100,1),"%",sep=""),
		main=paste(subset,"dataset"))

	if (names==TRUE) {
		text(pca$x[,1], pca$x[,2], labels.dat )

	}


	legend(pos.leg, groups, col=levels(as.factor(color)),pch=19)
### ojuuuu, problemes amb el label!!!

}

#plotPca(data.obj,"miRNA",col.color="age")



boxplotSamples <- function (obj, subset, col.color=1, las=1, colors=c("turquoise","violet")) {
	nf<-layout(mat=matrix(c(1),ncol=1,nrow=1))
	par(mar=c(5.1, 4.1, 4.1, 2.1))
#	plot.new()
	if (subset=="miRNA") {
		if (length(levels(as.factor(obj@pheno.miRNA[,col.color])))==2) {
			color<-as.character(factor(as.numeric(as.factor(obj@pheno.miRNA[,col.color])),labels=colors))
		} else {
			color<-as.character(factor(as.numeric(as.factor(obj@pheno.miRNA[,col.color]))+1))
		}
		boxplot(obj@dat.miRNA,col=color,las=las)
	}
	if (subset=="mRNA") {
		if (length(levels(as.factor(obj@pheno.mRNA[,col.color])))==2) {
			color<-as.character(factor(as.numeric(as.factor(obj@pheno.mRNA[,col.color])),labels=colors))
		} else {
			color<-as.character(factor(as.numeric(as.factor(obj@pheno.mRNA[,col.color]))+1))
		}
		boxplot(obj@dat.mRNA,col=color,las=las)
	}


}

#comprovar lo dels "..." que no sé si funciona
boxplotCorrelation <- function (obj, miRNA, mRNA, col.color=1, pos.leg="topright", colors=c("turquoise","violet"), ...) {
  
  #make 4-figure layout
	nf<-layout(mat=matrix(c(1:4),ncol=2,nrow=2),heights=c(1,3),widths=c(1,3))
	par(mar=c(5.1, 4.1, 1.1, 2.1))
	plot.new()
	

	if (length(levels(as.factor(obj@pheno.miRNA[,col.color])))==2) {
		color<-as.character(factor(as.numeric(as.factor(obj@pheno.miRNA[,col.color])),labels=colors))
		colors.plot <- as.character(factor(as.numeric(as.factor(obj@pheno.mRNA[,col.color]))+1,labels=colors))
	} else {
		color<-as.character(factor(as.numeric(as.factor(obj@pheno.miRNA[,col.color]))+1))
		colors.plot <- as.character(factor(as.numeric(as.factor(obj@pheno.mRNA[,col.color]))+1))
	}

	print(color)

	boxplot(obj@dat.mRNA[mRNA,]~obj@pheno.mRNA[,col.color],las=2,col=levels(factor(color)), ...)
	boxplot(obj@dat.miRNA[miRNA,]~obj@pheno.miRNA[,col.color],horizontal=TRUE,las=2,col=levels(factor(color)), ...)

	plot(obj@dat.miRNA[miRNA,],obj@dat.mRNA[mRNA,],xlab=paste(miRNA," expression",sep=""),ylab=paste(mRNA," expression",sep=""),col=colors.plot,pch=19, ...)
	legend(pos.leg,levels(as.factor(obj@pheno.mRNA[,col.color])),col=levels(factor(colors.plot)),pch=19, ...)
	abline(lm(obj@dat.mRNA[mRNA,]~obj@dat.miRNA[miRNA,]), ...)
	
	#return to previous layout
	nf<-layout(mat=matrix(c(1),ncol=1,nrow=1),heights=c(1),widths=c(1))
	
}



addCorrelation.R<- function (obj,method="pearson",subset.miRNA=obj@sig.miRNA,subset.mRNA=obj@sig.mRNA,common=NULL, d.influences=FALSE, alternative="less", kfold.cv=FALSE) {

	obj@net<-data.frame()
	obj@info[["pcomb.method"]]<-NULL
	obj@info[["padjust.method"]]<-NULL
	obj@info[["cor.alternative.hypothesis"]]<-NULL


	#if no differential expression has been performed
	if (is.null(subset.miRNA)) {
		subset.miRNA<-rownames(obj@dat.miRNA)
	}
	if (is.null(subset.mRNA)) {
		subset.mRNA<-rownames(obj@dat.mRNA)
	}

	#select always common samples
	cat("Correlating miRNA and mRNA\n")
	if (is.null(common)) {
		common<-intersect(colnames(obj@dat.mRNA),colnames(obj@dat.miRNA))
	}


	#seleccionar subset dels comuns
	if (length(subset.miRNA)>1) {
		subset.miRNA.sel<-intersect(subset.miRNA,rownames(obj@dat.miRNA))		
		miRNA.data<-obj@dat.miRNA[subset.miRNA.sel,common]
		if (!all(subset.mRNA==obj@sig.mRNA)) {obj@info[["miRNA.criteria"]]<-"manual!"} else {if ((obj@info[["miRNA.criteria"]][1]=="manual!") | (obj@info[["miRNA.criteria"]][1]=="All miRNA")) {obj@info[["miRNA.criteria"]]<-"stored for addSig (please check)"}}
	}else {miRNA.data<-obj@dat.miRNA[,common]
		if (class(miRNA.data)=="numeric") {
			miRNA.data<-t(as.matrix(miRNA.data))
			rownames(miRNA.data)<-rownames(obj@dat.miRNA)
			}
		subset.miRNA<-rownames(miRNA.data)
		obj@info[["miRNA.criteria"]]<-"All miRNA"}

	if (length(subset.mRNA)>1) {
		subset.mRNA.sel<-intersect(subset.mRNA,rownames(obj@dat.mRNA))		
		mRNA.data<-obj@dat.mRNA[subset.mRNA.sel,common]
		if (!all(subset.mRNA==obj@sig.mRNA)) {obj@info[["mRNA.criteria"]]<-"manual!"} else {if ((obj@info[["mRNA.criteria"]][1]=="manual!") | (obj@info[["mRNA.criteria"]][1]=="All mRNA")) {obj@info[["mRNA.criteria"]]<-"stored for addSig (please check)"}}
	}else {mRNA.data<-obj@dat.mRNA[,common]	
		if (class(mRNA.data)=="numeric") {
			mRNA.data<-t(as.matrix(mRNA.data))
			rownames(mRNA.data)<-rownames(obj@dat.mRNA)
			}
	
		subset.mRNA<-rownames(mRNA.data)
		obj@info[["mRNA.criteria"]]<-"All mRNA"}
	
	
	
	

	### fer les correlacions
	correlation.matrix<-matrix(NA,nrow=nrow(miRNA.data),ncol=nrow(mRNA.data))
	colnames(correlation.matrix)<-rownames(mRNA.data)
	rownames(correlation.matrix)<-rownames(miRNA.data)
	pval.matrix<-matrix(NA,nrow=nrow(miRNA.data),ncol=nrow(mRNA.data))
	colnames(pval.matrix)<-rownames(mRNA.data)
	rownames(pval.matrix)<-rownames(miRNA.data)

	if ((d.influences)) {
	someData <- rep(0, nrow(miRNA.data)*nrow(mRNA.data)*length(common))

	inf.matr <- array(someData, c(nrow(miRNA.data), nrow(mRNA.data), length(common)))

	for (i in 1:nrow(miRNA.data)) {
		for (j in 1:nrow(mRNA.data)) {
			#x<-as.numeric(miRNA.data[i,])
			#y<-as.numeric(mRNA.data[j,])
			x<-miRNA.data[i,]
			y<-mRNA.data[j,]
	#		print(x)
	#		print(y)
	#		print(paste("Comparativa",i,"versus",j))
			cor.t<-cor.test(x,y,method=method,alternative=alternative)
			correlation.matrix[i,j]<-cor.t$estimate
			pval.matrix[i,j]<-cor.t$p.value
			inf.matr[i,j,]<-cooks.distance(lm(mRNA.data[j,]~miRNA.data[i,]))	

		}
	#	print(i)
	}
	obj@info[["influenting.sample"]]<-inf.matr

	}


	if (!(d.influences)) {

	for (i in 1:nrow(miRNA.data)) {
		for (j in 1:nrow(mRNA.data)) {
			#x<-as.numeric(miRNA.data[i,])
			#y<-as.numeric(mRNA.data[j,])
			x<-miRNA.data[i,]
			y<-mRNA.data[j,]
	#		print(x)
	#		print(y)
	#		print(paste("Comparativa",i,"versus",j))
			cor.t<-cor.test(x,y,method=method,alternative=alternative)
			correlation.matrix[i,j]<-cor.t$estimate
			pval.matrix[i,j]<-cor.t$p.value
		}
	#	print(i)
	}

	}



	if (kfold.cv) {


	for (i in 1:nrow(miRNA.data)) {
		for (j in 1:nrow(mRNA.data)) {

			x<-miRNA.data[i,]
			y<-mRNA.data[j,]

			my.cors <- vector()

			print("hey")
			for (z in 1:10) {

				keep <- sample(1:length(x),length(x)/10*9-1)
				print(keep)
				my.cors[z]<-cor(x[keep],y[keep],method=method)


			}

			print(my.cors)
			correlation.matrix[i,j]<-mean(my.cors)
			pval.matrix[i,j]<-sd(my.cors)
		}
	#	print(i)
	}






	}





	obj@cor<-correlation.matrix
	obj@pval<-pval.matrix
	#obj@common<-common
	obj@info[["correlation.type"]]<-method
	obj@info[["correlation.function.used"]]<-"correlation"
	obj@info[["correlation.samples.used"]]<-common
	obj@info[["cor.alternative.hypothesis"]]<-alternative
	return(obj)
}




 ## addcorrelation before removing d.influences
addCorrelation<- function (obj,method="pearson",subset.miRNA=obj@sig.miRNA,subset.mRNA=obj@sig.mRNA,common=NULL, alternative="less", voom=FALSE) {

	obj@net<-data.frame()
	obj@info[["pcomb.method"]]<-NULL
	obj@info[["padjust.method"]]<-NULL
	obj@info[["cor.alternative.hypothesis"]]<-NULL


	cat("Correlating miRNA and mRNA\n")
	if (is.null(common)) {
		common<-intersect(colnames(obj@dat.mRNA),colnames(obj@dat.miRNA))
	}
	
	#seleccionar subset dels comuns
	if (length(subset.miRNA)>1) {
				
		subset.miRNA.sel<-intersect(subset.miRNA,rownames(obj@dat.miRNA))		
		miRNA.data<-obj@dat.miRNA[subset.miRNA.sel,common]
				
		if (!all(subset.mRNA==obj@sig.mRNA)) {
			obj@info[["miRNA.criteria"]]<-"manual!"
		} else {
			if ((obj@info[["miRNA.criteria"]][1]=="manual!") | (obj@info[["miRNA.criteria"]][1]=="All miRNA")) {
				obj@info[["miRNA.criteria"]]<-"stored for addSig (please check)"
			}
		}
	} else {
		miRNA.data<-obj@dat.miRNA[,common]
		if (class(miRNA.data)=="numeric") {
			miRNA.data<-t(as.matrix(miRNA.data))
			rownames(miRNA.data)<-rownames(obj@dat.miRNA)
		}
		subset.miRNA<-rownames(miRNA.data)
		obj@info[["miRNA.criteria"]]<-"All miRNA"
	}
	

	if (length(subset.mRNA)>1) {
		subset.mRNA.sel<-intersect(subset.mRNA,rownames(obj@dat.mRNA))		
		mRNA.data<-obj@dat.mRNA[subset.mRNA.sel,common]
		if (!all(subset.mRNA==obj@sig.mRNA)) {obj@info[["mRNA.criteria"]]<-"manual!"} else {if ((obj@info[["mRNA.criteria"]][1]=="manual!") | (obj@info[["mRNA.criteria"]][1]=="All mRNA")) {obj@info[["mRNA.criteria"]]<-"stored for addSig (please check)"}}
	}else {mRNA.data<-obj@dat.mRNA[,common]	
		if (class(mRNA.data)=="numeric") {
			mRNA.data<-t(as.matrix(mRNA.data))
			rownames(mRNA.data)<-rownames(obj@dat.mRNA)
			}
	
		subset.mRNA<-rownames(mRNA.data)
		obj@info[["mRNA.criteria"]]<-"All mRNA"}
	
	
	
	#### en cas de DESeq o edgeR, dividir els pesos
	if (obj@info[["diffexp.miRNA.method"]][1] %in% c("DESeq","edgeR") | voom=TRUE) {
	  if (!voom) {
	    miRNA.data<-t ( t(miRNA.data) / obj@info[["weights.miRNA"]][colnames(miRNA.data)] )
	  } else {
	    voom.d<-voom(obj@dat.miRNA)$E
	    rownames(voom.d)<-rownames(obj@dat.miRNA)
	    colnames(voom.d)<-colnames(obj@dat.miRNA)
	    miRNA.data<-voom.d[rownames(miRNA.data),colnames(miRNA.data)]
	  }
	}
	
	if (obj@info[["diffexp.mRNA.method"]][1] %in% c("DESeq","edgeR")  | voom=TRUE) {
	  if (!voom) {
	    mRNA.data<-t ( t(mRNA.data) / obj@info[["weights.mRNA"]][colnames(mRNA.data)] )
	  } else {
	    voom.d<-voom(obj@dat.mRNA)$E
	    rownames(voom.d)<-rownames(obj@dat.mRNA)
	    colnames(voom.d)<-colnames(obj@dat.mRNA)
	    mRNA.data<-voom.d[rownames(mRNA.data),colnames(mRNA.data)]
	  }
	}
	
	
	##### if the method is BaySeq, warning that the method is still in test
	if (obj@info[["diffexp.miRNA.method"]][1] %in% c("baySeq")) {
	  print("BaySeq still in test, try other method for NGS")
	}
	if (obj@info[["diffexp.mRNA.method"]][1] %in% c("baySeq")) {
	  print("BaySeq still in test, try other method for NGS")
	}
	
	
	
	

	### fer les correlacions
	correlation.matrix<-matrix(NA,nrow=nrow(miRNA.data),ncol=nrow(mRNA.data))
	colnames(correlation.matrix)<-rownames(mRNA.data)
	rownames(correlation.matrix)<-rownames(miRNA.data)
	pval.matrix<-matrix(NA,nrow=nrow(miRNA.data),ncol=nrow(mRNA.data))
	colnames(pval.matrix)<-rownames(mRNA.data)
	rownames(pval.matrix)<-rownames(miRNA.data)

	d.influences<-FALSE

	if ((d.influences)) {
	someData <- rep(0, nrow(miRNA.data)*nrow(mRNA.data)*length(common))

	inf.matr <- array(someData, c(nrow(miRNA.data), nrow(mRNA.data), length(common)))

	for (i in 1:nrow(miRNA.data)) {
		for (j in 1:nrow(mRNA.data)) {
			#x<-as.numeric(miRNA.data[i,])
			#y<-as.numeric(mRNA.data[j,])
			x<-miRNA.data[i,]
			y<-mRNA.data[j,]
	#		print(x)
	#		print(y)
	#		print(paste("Comparativa",i,"versus",j))
			cor.t<-cor.test(x,y,method=method,alternative=alternative)
			correlation.matrix[i,j]<-cor.t$estimate
			pval.matrix[i,j]<-cor.t$p.value
			inf.matr[i,j,]<-cooks.distance(lm(mRNA.data[j,]~miRNA.data[i,]))	

		}
	#	print(i)
	}
	obj@info[["influenting.sample"]]<-inf.matr

	}


	if (!(d.influences)) {


	correlation.matrix<-matrix(0.0,nrow=nrow(miRNA.data),ncol=nrow(mRNA.data))
	colnames(correlation.matrix)<-rownames(mRNA.data)
	rownames(correlation.matrix)<-rownames(miRNA.data)
	pval.matrix<-matrix(0.0,nrow=nrow(miRNA.data),ncol=nrow(mRNA.data))
	colnames(pval.matrix)<-rownames(mRNA.data)
	rownames(pval.matrix)<-rownames(miRNA.data)



out<-.C("pearson",
   X = as.vector(as.numeric(t(mRNA.data))),
   Y = as.vector(as.numeric(t(miRNA.data))),
   pnrx = as.integer(nrow(mRNA.data)),
   pnry = as.integer(nrow(miRNA.data)),   
   pnc = as.integer(ncol(miRNA.data)),   
   corr = as.vector(as.numeric(correlation.matrix)),
   pval = as.vector(as.numeric(pval.matrix))
)

#### GET RESULTS
# matrix of correlations (r)
correlation.matrix <- matrix(out$corr,nrow=nrow(miRNA.data),ncol=nrow(mRNA.data))
rownames(correlation.matrix)<-rownames(miRNA.data)
colnames(correlation.matrix)<-rownames(mRNA.data)

# matrix of p-values (P)
pval.matrix <- matrix(out$pval,nrow=nrow(miRNA.data),ncol=nrow(mRNA.data))

if (alternative=="less") {
  pval.matrix <- pt(pval.matrix, ncol(mRNA.data)-2)
}
if (alternative=="greater") {
  pval.matrix <- 1-pt(pval.matrix, ncol(mRNA.data)-2)
}
if (alternative=="two.sided") {
  pval.matrix.two <- matrix(abs(out$pval),nrow=nrow(miRNA.data),ncol=nrow(mRNA.data))
  pval.matrix <- (1-pt(pval.matrix.two, ncol(mRNA.data)-2))*2
  pval.matrix <- apply(pval.matrix,c(1,2),min,1)
}

rownames(pval.matrix)<-rownames(miRNA.data)
colnames(pval.matrix)<-rownames(mRNA.data)

}


	obj@cor<-correlation.matrix
	obj@pval<-pval.matrix
	#obj@common<-common
	obj@info[["correlation.type"]]<-method
	obj@info[["correlation.function.used"]]<-"correlation"
	obj@info[["correlation.samples.used"]]<-common
	obj@info[["cor.alternative.hypothesis"]]<-alternative
	return(obj)
}


addGlmnet <- function (obj , response = "mRNAs", alpha=0.5, upper.limit=0, cluster="manual", plot=TRUE) {

	#library(glmnet)	

#obj<-data.txell
#response = "mRNAs"
#alpha=0.5
#upper.limit=0


	x<-t(obj@dat.miRNA[obj@sig.miRNA,])
	y<-t(obj@dat.mRNA[obj@sig.mRNA,])
	fullcoefs<-matrix(0,ncol=ncol(y), nrow=ncol(x))
	rownames(fullcoefs)<-colnames(x)
	colnames(fullcoefs)<-colnames(y)


	if (response=="mRNAs" | response=="mrnas") {
		res<-cv.glmnet(x, y,alpha=alpha,family = "mgaussian", upper.limit=upper.limit)

		opt <- coef(res, s="lambda.min")
		for (i in 1:ncol(fullcoefs)) {
			fullcoefs[opt[[i]]@i[-1],i] <-as.vector(opt[[i]]@x)[-1]
 			print(i)
		}
		obj@cor<-fullcoefs
	}

	if (response=="mRNA" | response=="mrna") {
	
		for (i in 1:ncol(y)) {
			res<-cv.glmnet(x, y[,i],alpha=alpha,family = "gaussian", upper.limit=upper.limit)
			fullcoefs[,i] <-coef(res, s="lambda.min")[-1]
		}
		obj@cor<-fullcoefs


	}


	if (response=="mRNA-clust" | response=="mrna-clust") {

		#clusterization
		print(" Making clusters")


		## manual identification:
		if (cluster=="manual") {
		
			print("Select miRNA clusters")

			scaled.mirnas<-t(x)
			for (i in 1:nrow(x)) {
				scaled.mirnas[i,]<-(t(x)[i,]-mean(t(x)[i,]))/sd(t(x)[i,])
	
			}
	
			d_dist <- dist(as.matrix(scaled.mirnas))   # find distance matrix 
			plot(hclust(d_dist)) 
			mirna.clusters <- identify(hclust(d_dist),method="ward.D")
		#	print(mirna.clusters)
		#	mirna.class<-unlist(mirna.clusters)[colnames(x)]
			mirna.clusters.o<-mirna.clusters
			for (i in 1:length(mirna.clusters)) {
				mirna.clusters[[i]]<-names(mirna.clusters[[i]])
			}

		#	print(mirna.class)

			print("Select mRNA clusters")
			d_dist <- dist(as.matrix(t(y)))   # find distance matrix 
			plot(hclust(d_dist)) 
			mrna.clusters <- identify(hclust(d_dist))
		#	print(mrna.clusters)
		#	mrna.class<-unlist(mrna.clusters)[colnames(x)]
			mrna.clusters.o<-mrna.clusters
			for (i in 1:length(mrna.clusters)) {
				mrna.clusters[[i]]<-names(mrna.clusters[[i]])
			}
	

			if (plot==TRUE) {
				mirna.class.o<-list()
				print(mirna.clusters.o)
				for (i in 1:length(mirna.clusters.o)) {
					mirna.class.o[[i]]<-rep(i,length(mirna.clusters.o[[i]]))
					names(mirna.class.o[[i]])<-names(mirna.clusters.o[[i]])
				}
				mirna.class<-unlist(mirna.class.o)[colnames(x)]

				mrna.class.o<-list()
				print(mrna.clusters.o)
				for (i in 1:length(mrna.clusters.o)) {
					mrna.class.o[[i]]<-rep(i,length(mrna.clusters.o[[i]]))
					names(mrna.class.o[[i]])<-names(mrna.clusters.o[[i]])
				}
				mrna.class<-unlist(mrna.class.o)[colnames(x)]



				pca <- prcomp(t(x), scale=T)
				plot(pca$x[,1],pca$x[,2],pch=19,col=as.numeric(mirna.class),main="Clustered miRNAs")
				x11()


				plotHeatmap(obj,"miRNA",grouping.row=mirna.class)
				pca <- prcomp(t(y), scale=T)
				x11()		
				plot(pca$x[,1],pca$x[,2],pch=19,col=as.numeric(mrna.class),main="Clustered mRNAs")
				x11()
				plotHeatmap(obj,"mRNA",grouping.row=mrna.class)
			}





		}



		if (cluster=="automatic") {
			#library(mclust)
			d_clust <- Mclust(as.matrix(t(x)), G=1:ceiling(sqrt(ncol(x))*1.5))
			m.best <- dim(d_clust$z)[2]
			#cat("model-based optimal number of clusters:", m.best, "\n")
			mirna.class<-d_clust$classification
			print(paste("Number of miRNA clusters:",m.best))
			print(table(mirna.class))

			mirna.clusters<-list()
			for (i in 1:length(table(mirna.class))) {
				mirna.clusters[[i]]<-names(mirna.class[which(mirna.class==i)])

			}

			d_clust <- Mclust(as.matrix(t(y)), G=1:ceiling(sqrt(ncol(y))*0.75),initialization=list(subset=sample(1:nrow(t(y)), size=ncol(y)/1)) )
			m.best <- dim(d_clust$z)[2]
			print(paste("Number of mRNA clusters:",m.best))
			mrna.class<-d_clust$classification
			print(table(mrna.class))

			#cat("model-based optimal number of clusters:", m.best, "\n")
			mrna.class<-d_clust$classification

			mrna.clusters<-list()
			for (i in 1:length(table(mrna.class))) {
				mrna.clusters[[i]]<-names(mrna.class[which(mrna.class==i)])

			}


			if (plot==TRUE) {
				pca <- prcomp(t(x), scale=T)
				plot(pca$x[,1],pca$x[,2],pch=19,col=as.numeric(mirna.class),main="Clustered miRNAs")
				x11()	
				plotHeatmap(obj,"miRNA",grouping.row=mirna.class)
				pca <- prcomp(t(y), scale=T)
				x11()		
				plot(pca$x[,1],pca$x[,2],pch=19,col=as.numeric(mrna.class),main="Clustered mRNAs")
				x11()
				plotHeatmap(obj,"mRNA",grouping.row=mrna.class)
			}



		}


		print("Computing relations")


		for (ii in 1:length(mirna.clusters)) {
			for (jj in 1:length(mrna.clusters)) {

				print(paste("MiRNA cluster",ii,", mRNA cluster",jj))

				res<-cv.glmnet(x[,mirna.clusters[[ii]]], y[,mrna.clusters[[jj]]],alpha=alpha,family = "mgaussian", upper.limit=upper.limit)
				opt <- coef(res, s="lambda.min")

				for (zzz in 1:length(mrna.clusters[[jj]])) {
					fullcoefs[ mirna.clusters[[ii]], mrna.clusters[[jj]][zzz]  ] <- as.vector(opt[[mrna.clusters[[jj]][zzz]]][-1])
					#fullcoefs[ names(mirna.clusters[[ii]]), names(mrna.clusters[[jj]][zzz])  ] <- as.vector(opt[[mrna.clusters[[jj]][zzz]]][-1])			
				}
			}
		}


		obj@cor<-fullcoefs

#source("http://www.r-statistics.com/wp-content/uploads/2012/01/source_https.r.txt") # Making sure we can source code from github
#source_https("http://raw.github.com/talgalili/R-code-snippets/master/clustergram.r")
 
#data(iris)
#set.seed(250)
#par(cex.lab = 1.5, cex.main = 1.2)
#Data <- scale(iris[,-5]) # notice I am scaling the vectors)
#clustergram(Data, k.range = 2:8, line.width = 0.004) # notice how I am using line.width.  Play with it on your problem, according to the scale of Y.

	}


	if (response=="miRNA" | response=="mirna") {

	
		for (i in 1:ncol(x)) {
			res<-cv.glmnet(y, x[,i],alpha=alpha,family = "gaussian", upper.limit=upper.limit)
			fullcoefs[i,] <-coef(res, s="lambda.min")[-1]
		}
		obj@cor<-fullcoefs

	}

	return(obj)

}

#a<-addGlmnet(data.obj,response="mRNAs")
#b<-addGlmnet(data.obj,response="mRNA")
#c<-addGlmnet(data.obj,response="miRNA")

#	x<-t(data.obj@dat.miRNA[data.obj@sig.miRNA,])
#	y<-t(data.obj@dat.mRNA[data.obj@sig.mRNA,])


#res<-cv.glmnet(  as.matrix(x),  as.vector(y[,1]),alpha=alpha,family = "gaussian", upper.limit=upper.limit)



addSurv <- function (obj, dataset, time, event, adjusting=NULL) {

	#require(survival)
	if (dataset=="miRNA") {
		time.v <- obj@pheno.miRNA[,time]
		event.v <- as.numeric(obj@pheno.miRNA[,event])
		data<-obj@dat.miRNA
		if (!is.null(adjusting)) {
			adjust.v <- obj@pheno.miRNA[,adjusting]
			print("c")
		}
		else {
			adjust.v <- NA
			print("d")
		}
	}
	if (dataset=="mRNA") {
		time.v <- obj@pheno.mRNA[,time]
		event.v <- as.numeric(obj@pheno.mRNA[,event])
		data<-obj@dat.mRNA
		if (!is.null(adjusting)) {
			adjust.v <- obj@pheno.mRNA[,adjusting]
			print("a")
		}
		else {
			adjust.v <- NA
			print("b")
		}
	}

	comp.surv <- function (x, time.f, event.f, adjust.f) {

		if (!is.na(adjust.f)) {
			reg<-as.matrix(data.frame(x=x,adjust.f))
		}
		if (is.na(adjust.f)) {
			reg<-as.matrix(data.frame(x=x))
		}

		mod.uni<-coxph(Surv(time.f,event.f)~reg ) 
		res<-as.vector(data.frame(summary(mod.uni)$coefficients)[1,c(1,5)])
	}

	res.surv<-data.frame(matrix(unlist(apply(data,1,comp.surv,time.v,event.v,adjust.v)), ncol=2, byrow=T),stringsAsFactors=FALSE)
	colnames(res.surv)<-c("coef","pval")
	rownames(res.surv)<-rownames(data)
	res.surv$adj.pval<-p.adjust(res.surv$pval,method="BH")

	if (dataset=="miRNA") {
		obj@diffexp.miRNA<-res.surv
		obj@info[["miRNA.diffexp.method"]][1]<-"survival"
		obj@info[["miRNA.diffexp.method"]][2]<-time
		obj@info[["miRNA.diffexp.method"]][3]<-event
	}
	if (dataset=="mRNA") {
		obj@diffexp.mRNA<-res.surv
		obj@info[["mRNA.diffexp.method"]][1]<-"survival"
		obj@info[["mRNA.diffexp.method"]][2]<-time
		obj@info[["mRNA.diffexp.method"]][3]<-event
	}

	return(obj)

}



plotSurv <- function (obj, subset, item, time, event) {

	#require(survival)
	if (subset=="miRNA") {
		time.v <- obj@pheno.miRNA[,time]
		event.v <- as.numeric(obj@pheno.miRNA[,event])
		item.v<-obj@dat.miRNA[item,]
		item.cat<-as.character(item.v)
		item.cat[which(item.v>=median(item.v))]<-"High"
		item.cat[which(item.v<median(item.v))]<-"Low"
		item.cat<-factor(item.cat)
		item.cat<-relevel(item.cat,"Low")

	}
	if (subset=="mRNA") {
		time.v <- obj@pheno.mRNA[,time]
		event.v <- as.numeric(obj@pheno.mRNA[,event])
		item.v<-obj@dat.mRNA[item,]
		item.cat<-as.character(item.v)
		item.cat[which(item.v>=median(item.v))]<-"High"
		item.cat[which(item.v<median(item.v))]<-"Low"
		item.cat<-factor(item.cat)
		item.cat<-relevel(item.cat,"Low")
	}

	plot(survfit(Surv(time.v,event.v)~item.cat ), main=item ,col=c(1,2), xlab="Time", ylab="Survival", bty="l")


	# draw an axis on the left

	# draw an axis on the right, with smaller text and ticks

	legend("topright",c("Low","High"),col=c(1,2),lty=1)

}










addNet <- function (obj) {

	obj@info[["pcomb.method"]]<-NULL
	obj@info[["padjust.method"]]<-NULL

	cat("Converting to net\n")
	# convertir a fitxer per xarxa

	if (!is.null(rownames(obj@cor)) & !is.null(rownames(obj@pval))) {

		corr<-t(obj@cor)
		pvals<-t(obj@pval)

		obj@net<-data.frame(miRNA = rep(colnames(corr), each = nrow(corr)), 
        	   mRNA = rep(rownames(corr), ncol(corr)), 
        	   cor = as.vector(corr),
		   pval = as.vector(pvals))

	} else {


		if (!is.null(rownames(obj@cor)) ) {

			corr<-t(obj@cor)

			obj@net<-data.frame(miRNA = rep(colnames(corr), each = nrow(corr)), 
        		   mRNA = rep(rownames(corr), ncol(corr)), 
        		   cor = as.vector(corr))

		} else {

		
			if (!is.null(obj@sig.miRNA)) {
				mirnas <- obj@sig.miRNA
			} else {
				mirnas <- rownames(obj@dat.miRNA)
			}
			if (!is.null(obj@sig.mRNA)) {
				mrnas <- obj@sig.mRNA
			} else {
				mrnas <- rownames(obj@dat.mRNA)
			}
	
			obj@net<-data.frame(miRNA = rep(mirnas, each = length(mrnas)), 
        		   mRNA = rep(mrnas, length(mirnas)))

		}

	}



	rownames(obj@net)<-paste(obj@net$miRNA,obj@net$mRNA,sep=":")

	return(obj)	

#system.time(a2<-melt(m))
#system.time(a2.t<-melt(t(m)))
#a3<-m
#dimnames(a3) <- list( One=colnames(a3), Two=rownames(a3) )
#a3.3<-as.data.frame( as.table(a3) )

}



# This way is faster then reshape! :D
#d <- 
#as.matrix(read.table(text = "
#    month    Q1    Q2   Q3    Q4    Q5    Q6    Q7    Q8    Q9   Q10   Q11   Q12   Q13         
#X10    10  7.04  8.07  9.4  8.17  9.39  8.13  9.43  9.06  8.59  9.37  9.79  8.47  8.86  
#X11    11 12.10 11.50 12.6 13.70 11.90 11.50 13.10 17.20 19.00 14.60 13.70 13.20 16.10  
#X12    12 24.00 22.00 22.2 20.50 21.60 22.50 23.10 23.30 30.50 34.10 36.10 37.40 28.90 
#X1      1 18.30 16.30 16.2 14.80 16.60 15.40 15.20 14.80 16.70 14.90 15.00 13.80 15.90  
#X2      2 16.70 14.40 15.3 14.10 15.50 16.70 15.20 16.10 18.00 26.30 28.00 31.10 34.20",
#header=TRUE))


#> system.time(res <- reshape(as.data.frame(d), idvar="month", timevar="day", 
#+                varying = -1, direction = "long", sep = "")
#+ )
#   user  system elapsed 
#  0.004   0.000   0.004 
#> 
#> system.time(res2 <- data.frame(mirna=rep(colnames(d),each=nrow(d)),
#+ mrna = rep(rownames(d),ncol(d)),
#+ iuju = as.vector(d))
#+ )
#   user  system elapsed 
#  0.000   0.000   0.001 


addFoldchanges <- function (obj, add.pvals=FALSE) {

	if (all(dim(obj@diffexp.miRNA)!=0)) {

		if (obj@info[["miRNA.diffexp.method"]][1]=="time.point" | obj@info[["miRNA.diffexp.method"]][1]=="linear.regression" | obj@info[["miRNA.diffexp.method"]][1]=="anova") {
			if (obj@info[["mRNA.diffexp.method"]][1]=="anova") {
				obj@net$ss.prop.miRNA<-obj@diffexp.miRNA[as.character(obj@net$miRNA),1]
			} else {
				obj@net$slope.miRNA<-obj@diffexp.miRNA[as.character(obj@net$miRNA),1]
			}
			obj@net$meanExp.miRNA<-obj@diffexp.miRNA[as.character(obj@net$miRNA),"meanExp"]
		} else {
			obj@net$logratio.miRNA<-obj@diffexp.miRNA[as.character(obj@net$miRNA),"logratio"]
			obj@net$meanExp.miRNA<-obj@diffexp.miRNA[as.character(obj@net$miRNA),"meanExp"]
		}
		if (add.pvals) {
			obj@net$adj.pval.miRNA<-obj@diffexp.miRNA[as.character(obj@net$miRNA),"adj.pval"]
			obj@net$pval.miRNA<-obj@diffexp.miRNA[as.character(obj@net$miRNA),"pval"]
		}

	}

	if (all(dim(obj@diffexp.mRNA)!=0)) {

		if (obj@info[["mRNA.diffexp.method"]][1]=="time.point" | obj@info[["mRNA.diffexp.method"]][1]=="linear.regression"| obj@info[["mRNA.diffexp.method"]][1]=="anova") {
			if (obj@info[["mRNA.diffexp.method"]][1]=="anova") {
				obj@net$ss.prop.mRNA<-obj@diffexp.mRNA[as.character(obj@net$mRNA),1]
			} else {
				obj@net$slope.mRNA<-obj@diffexp.mRNA[as.character(obj@net$mRNA),1]
			}
			obj@net$meanExp.mRNA<-obj@diffexp.mRNA[as.character(obj@net$mRNA),"meanExp"]
		} else {
			obj@net$logratio.mRNA<-obj@diffexp.mRNA[as.character(obj@net$mRNA),"logratio"]
			obj@net$meanExp.mRNA<-obj@diffexp.mRNA[as.character(obj@net$mRNA),"meanExp"]
		}
		if (add.pvals) {
			obj@net$adj.pval.mRNA<-obj@diffexp.mRNA[as.character(obj@net$mRNA),"adj.pval"]
			obj@net$pval.mRNA<-obj@diffexp.mRNA[as.character(obj@net$mRNA),"pval"]
		}
	}
	return(obj)
}






removeSamp <- function (obj, dataset, samples=NA, genes=NA, keep=FALSE) {

## add error messages if the matrix is too low

	if (dataset == "miRNA") {
		if (!is.na(genes[1])) {
			if (keep==FALSE) {obj@dat.miRNA<-obj@dat.miRNA[-which(rownames(obj@dat.miRNA) %in% genes),]}
			if (keep==TRUE) {
				if (length(genes)==1) {
					obj@dat.miRNA<-t(as.matrix(obj@dat.miRNA[which(rownames(obj@dat.miRNA) %in% genes),]))
					rownames(obj@dat.miRNA)<-genes
				} else {
					obj@dat.miRNA<-obj@dat.miRNA[which(rownames(obj@dat.miRNA) %in% genes),]

				}
			}
		}
		if (!is.na(samples[1])) {
			if (keep==FALSE) {
				if (nrow(obj@dat.miRNA)>1) {
					obj@dat.miRNA<-obj@dat.miRNA[,-which(colnames(obj@dat.miRNA) %in% samples)]
				} else {
					gene.names<-rownames(obj@dat.miRNA)
					obj@dat.miRNA<-t(as.matrix(obj@dat.miRNA[,-which(colnames(obj@dat.miRNA) %in% samples)]))
					rownames(obj@dat.miRNA)<-gene.names
				}
				obj@pheno.miRNA<-obj@pheno.miRNA[-which(rownames(obj@pheno.miRNA) %in% samples),]
			}

			if (keep==TRUE) {
				if (nrow(obj@dat.miRNA)>1) {
					obj@dat.miRNA<-obj@dat.miRNA[,which(colnames(obj@dat.miRNA) %in% samples)]
				} else {
					gene.names<-rownames(obj@dat.miRNA)
					obj@dat.miRNA<-t(as.matrix(obj@dat.miRNA[,which(colnames(obj@dat.miRNA) %in% samples)]))
					rownames(obj@dat.miRNA)<-gene.names
				}
				obj@pheno.miRNA<-obj@pheno.miRNA[which(rownames(obj@pheno.miRNA) %in% samples),]
			}
		}
	}

	if (dataset == "mRNA") {
		if (!is.na(genes[1])) {
			if (keep==FALSE) {obj@dat.mRNA<-obj@dat.mRNA[-which(rownames(obj@dat.mRNA) %in% genes),]}
			if (keep==TRUE) {
				if (length(genes)==1) {
					obj@dat.mRNA<-t(as.matrix(obj@dat.mRNA[which(rownames(obj@dat.mRNA) %in% genes),]))
					rownames(obj@dat.mRNA)<-genes
				} else {
					obj@dat.mRNA<-obj@dat.mRNA[which(rownames(obj@dat.mRNA) %in% genes),]
				}

			}
		}
		if (!is.na(samples[1])) {
			if (keep==FALSE) {
				if (nrow(obj@dat.mRNA)>1) {
					obj@dat.mRNA<-obj@dat.mRNA[,-which(colnames(obj@dat.mRNA) %in% samples)]
				} else {
					gene.names<-rownames(obj@dat.mRNA)
					obj@dat.mRNA<-t(as.matrix(obj@dat.mRNA[,-which(colnames(obj@dat.mRNA) %in% samples)]))
					rownames(obj@dat.mRNA)<-gene.names
				}
				obj@pheno.mRNA<-obj@pheno.mRNA[-which(rownames(obj@pheno.mRNA) %in% samples),]
			}

			if (keep==TRUE) {
				if (nrow(obj@dat.mRNA)>1) {
					obj@dat.mRNA<-obj@dat.mRNA[,which(colnames(obj@dat.mRNA) %in% samples)]
				} else {
					gene.names<-rownames(obj@dat.mRNA)
					obj@dat.mRNA<-t(as.matrix(obj@dat.mRNA[,which(colnames(obj@dat.mRNA) %in% samples)]))
					rownames(obj@dat.mRNA)<-gene.names
				}
				obj@pheno.mRNA<-obj@pheno.mRNA[which(rownames(obj@pheno.mRNA) %in% samples),]
			}
		}
	}
	return(obj)

}



addDiffexp <- function (obj, dataset, classes, method.dif="t.test", method.adj="BH", var.t.test=FALSE, trend = FALSE ) {

#	require(gtools)
	#per poder determinar els grups a comparar
#	require(RankProd)
	if (is.character(classes)) {
		if (dataset == "miRNA") {
			lev<-obj@pheno.miRNA[,classes]
			levi<-levels(as.factor(obj@pheno.miRNA[,classes]))
			lev1<-which(lev==1)
			lev2<-which(lev==0)

		}
		if (dataset == "mRNA") {
			lev<-obj@pheno.mRNA[,classes]
			levi<-levels(as.factor(obj@pheno.mRNA[,classes]))
			lev1<-which(lev==1)
			lev2<-which(lev==0)
		}

	} else {
		lev<-classes
		levi<-levels(as.factor(classes))
		lev1<-which(lev==levi[1])
		lev2<-which(lev==levi[2])
	}
	
	calc.lograt <- function (dats,lev1,lev2) {
		lograt<-mean(dats[lev1])-mean(dats[lev2])
		return(lograt)		
	}

	ttestfd <- function (dats,lev1,lev2) {
		t.test(dats[lev1],dats[lev2],var.equal=var.t.test)$p.value
	}

	wilcoxfd <- function (dats,lev1,lev2) {
		wilcox.test(dats[lev1],dats[lev2])$p.value
	}

	aov_pval <- function (dats, comp) {
		summary(aov(dats~comp))[[1]][1,5]
	}

	aov_ss <- function (dats, comp) {
		summary(aov(dats~comp))[[1]][1,2] / ( summary(aov(dats~comp))[[1]][2,2] + summary(aov(dats~comp))[[1]][1,2] )
	}


	if (dataset == "miRNA") {

		mirdat<-obj@dat.miRNA[,c(lev1,lev2)]

		if (method.dif=="anova" | method.dif=="max.var") {
			mirdat<-obj@dat.miRNA
		}

		rlev1<-1:length(lev1)
		rlev2<-((length(lev1)+1):(length(lev1)+length(lev2)))

		if (method.dif=="t.test" | method.dif=="wilcoxon") {
			lograt<-apply(mirdat,1,calc.lograt,rlev1,rlev2)	
			FC<-logratio2foldchange(lograt)
			if (method.dif=="t.test") {
				pval.uncor<-apply(mirdat,1,ttestfd,rlev1,rlev2)}
			if (method.dif=="wilcoxon") {
				pval.uncor<-apply(mirdat,1,wilcoxfd,rlev1,rlev2)}
			adj.pval<-p.adjust(pval.uncor,method=method.adj)
			obj@diffexp.miRNA<-data.frame(FC,logratio=lograt,pval=pval.uncor,adj.pval=adj.pval)
			rownames(obj@diffexp.miRNA)<-rownames(mirdat)
		}

###### !

		if (method.dif=="only.fc") {
			lograt<-apply(mirdat,1,calc.lograt,rlev1,rlev2)	
			FC<-logratio2foldchange(lograt)
			obj@diffexp.miRNA<-data.frame(FC,logratio=lograt)
			rownames(obj@diffexp.miRNA)<-rownames(mirdat)
		}




		if (method.dif=="limma") {
			em<-as.matrix(mirdat)
			grup<-as.factor(obj@pheno.miRNA[,classes][c(lev1,lev2)]+1)
			edesign<-model.matrix(~0+grup)
			efit<-lmFit(em,edesign)
			contr<-makeContrasts(dif=grup2-grup1,levels=edesign)
			efit2<-contrasts.fit(efit,contr)
			efit2<-eBayes(efit2,trend=trend)
			obj@diffexp.miRNA<-data.frame(FC=logratio2foldchange(as.vector(efit2$coefficients)),
				logratio=as.vector(efit2$coefficients),
				pval=as.vector(efit2$p.value),
				adj.pval=p.adjust(efit2$p.value,method=method.adj)
				)
			rownames(obj@diffexp.miRNA)<-rownames(mirdat)
		}


		if (method.dif=="rankprod") {
			RP.out<-RP(mirdat,obj@pheno.miRNA[,classes][c(lev1,lev2)])
			#invertir el FC class1/class2
			obj@diffexp.miRNA<-data.frame(FC=logratio2foldchange(as.vector(-RP.out$AveFC)),
				logratio=as.vector(-RP.out$AveFC),
				pval=apply(RP.out$pval,1,min),
				adj.pval=apply(RP.out$pfp,1,min))
			rownames(obj@diffexp.miRNA)<-rownames(mirdat)
		}


		if (method.dif=="DESeq") {
			conds <- obj@pheno.miRNA[,classes][c(lev1,lev2)]
			cds<-newCountDataSet(round(mirdat,0),conds)

			cds<-estimateSizeFactors(cds)
			cds<-estimateDispersions(cds,fitType="local")
			comparative <- nbinomTest (cds, "0", "1")

			obj@diffexp.miRNA<-data.frame(FC=comparative$foldChange,
				logratio=comparative$log2FoldChange,
				pval=comparative$pval,
				adj.pval=comparative$padj)
			rownames(obj@diffexp.miRNA)<-comparative$id
			
			obj@info[["miRNA.weights"]]<-cds@phenoData@data$sizeFactor
		}

#library(DESeq2)
#			cds<-newCountDataSet(round(mirdat,0),conds)
#ddsFull <- DESeqDataSet( cds, design = ~ 0 + conds)

		
		if (method.dif=="edgeR") {
			conds <- obj@pheno.miRNA[,classes][c(lev1,lev2)]
			d <- DGEList(counts=mirdat,group=factor(conds))	
			design.mat <- model.matrix(~ 0 + d$samples$group)

			d <- calcNormFactors(d)
			d1 <- estimateCommonDisp(d, verbose=T)
			d1 <- estimateTagwiseDisp(d1)

			colnames(design.mat) <- levels(d$samples$group)
			d2 <- estimateGLMCommonDisp(d,design.mat)
			d2 <- estimateGLMTrendedDisp(d2,design.mat, method="auto")
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
			d2 <- estimateGLMTagwiseDisp(d2,design.mat)

			et12 <- exactTest(d1, pair=c("0","1")) # compare groups 1 and 2

			obj@diffexp.miRNA<-data.frame(FC=logratio2foldchange(et12$table$logFC),
				logratio=et12$table$logFC,
				pval=et12$table$PValue,
				adj.pval=p.adjust(et12$table$PValue,method="BH"))
			rownames(obj@diffexp.miRNA)<-rownames(et12$table)
			
			obj@info[["miRNA.weights"]]<-d@.Data[[2]]$norm.factors
		}



		if (method.dif=="baySeq") {
			conds <- obj@pheno.miRNA[,classes][c(lev1,lev2)]
			cds<-newCountDataSet(round(mirdat,0),conds)

			cds<-estimateSizeFactors(cds)
			cds<-estimateDispersions(cds,fitType="local")
			comparative <- nbinomTest (cds, "0", "1")

			obj@diffexp.miRNA<-data.frame(FC=comparative$foldChange,
				logratio=comparative$log2FoldChange,
				pval=comparative$pval,
				adj.pval=comparative$padj)
			rownames(obj@diffexp.miRNA)<-comparative$id
		}

#		CD <- new("countData", data = mirdat, replicates = conds, libsizes = as.integer(apply(mirdat,2,sum)), groups = list(c1=conds))

#cl <- NULL
#CDP.NBML <- getPriors.NB(CD, samplesize = 1000, estimation = "QL", cl = cl)
#CDPost.NBML <- getLikelihoods(CDP.NBML, prs=pr/sum(pr),  pET = 'BIC', cl = cl)

#pr<-apply(mirdat,1,sum)

#prs=pr/sum(pr)

#prs=rep(1/nrow(mirdat),nrow(mirdat))

#rep(1/ncol(mirdat),ncol(mirdat))

#bayseq_de = topCounts(CDPost.NBML, group=1, number=20)

#CD@annotation <- as.data.frame(cname)


		if (method.dif=="anova") {
			pval.aov<-apply(mirdat,1,aov_pval,as.factor(obj@pheno.miRNA[,classes]))
			ss.aov<-apply(mirdat,1,aov_ss,as.factor(obj@pheno.miRNA[,classes]))
			obj@diffexp.miRNA<-data.frame(ss.prop=ss.aov,
				pval=pval.aov,
				adj.pval=p.adjust(pval.aov,method=method.adj))
			rownames(obj@diffexp.miRNA)<-rownames(mirdat)
		}


		if (method.dif=="max.var") {
			myvar <- apply(mirdat,1,var)
			conv.pvals<- abs((myvar-max(myvar))/max(myvar))

			obj@diffexp.miRNA<-data.frame(variance=myvar,
				pval=conv.pvals,
				adj.pval=p.adjust(conv.pvals,method=method.adj))
			rownames(obj@diffexp.miRNA)<-rownames(mirdat)

		}

	#	if (method.dif=="time-course") {
	#	linear <- 

	#	}

		obj@diffexp.miRNA$meanExp<-apply(mirdat,1,mean)
		if (method.dif != "anova" & method.dif != "max.var" & method.dif!="only.fc") {
			obj@diffexp.miRNA<-obj@diffexp.miRNA[,c		("FC","logratio","meanExp","pval","adj.pval")]
		}
		if (method.dif=="anova") {
			obj@diffexp.miRNA<-obj@diffexp.miRNA[,c		("ss.prop","meanExp","pval","adj.pval")]
		}
		if (method.dif=="max.var") {
			obj@diffexp.miRNA<-obj@diffexp.miRNA[,c		("variance","meanExp","pval","adj.pval")]
		}
		if (method.dif=="only.fc") {
			obj@diffexp.miRNA<-obj@diffexp.miRNA[,c("FC","logratio","meanExp")]
		}
		obj@info[["miRNA.diffexp.method"]][1]<-method.dif
		obj@info[["miRNA.diffexp.method"]][2]<-classes
	}

	if (dataset == "mRNA") {

		mrdat<-obj@dat.mRNA[,c(lev1,lev2)]

		if (method.dif=="anova" | method.dif=="max.var") {
			mrdat<-obj@dat.mRNA
		}


		rlev1<-1:length(lev1)
		rlev2<-((length(lev1)+1):(length(lev1)+length(lev2)))


		if (method.dif=="t.test" | method.dif=="wilcoxon") {
			lograt<-apply(mrdat,1,calc.lograt,rlev1,rlev2)	
			FC<-logratio2foldchange(lograt)
			if (method.dif=="t.test") {
				pval.uncor<-apply(mrdat,1,ttestfd,rlev1,rlev2)}
			if (method.dif=="wilcoxon") {
				pval.uncor<-apply(mrdat,1,wilcoxfd,rlev1,rlev2)}
			adj.pval<-p.adjust(pval.uncor,method=method.adj)
			obj@diffexp.mRNA<-data.frame(FC,logratio=lograt,pval=pval.uncor,adj.pval=adj.pval)
			rownames(obj@diffexp.mRNA)<-rownames(mrdat)
		}

		if (method.dif=="only.fc") {
			lograt<-apply(mrdat,1,calc.lograt,rlev1,rlev2)	
			FC<-logratio2foldchange(lograt)
			obj@diffexp.mRNA<-data.frame(FC,logratio=lograt)
			rownames(obj@diffexp.mRNA)<-rownames(mrdat)
		}


		if (method.dif=="limma") {
			em<-as.matrix(mrdat)
			grup<-as.factor(obj@pheno.mRNA[,classes][c(lev1,lev2)]+1)
			edesign<-model.matrix(~0+grup)
			efit<-lmFit(em,edesign)
			contr<-makeContrasts(dif=grup2-grup1,levels=edesign)
			efit2<-contrasts.fit(efit,contr)
			efit2<-eBayes(efit2,trend=trend)
			obj@diffexp.mRNA<-data.frame(FC=logratio2foldchange(as.vector(efit2$coefficients)),
				logratio=as.vector(efit2$coefficients),
				pval=as.vector(efit2$p.value),
				adj.pval=p.adjust(efit2$p.value,method=method.adj)
				)
			rownames(obj@diffexp.mRNA)<-rownames(mrdat)
		}

		if (method.dif=="rankprod") {
			RP.out<-RP(mrdat,obj@pheno.mRNA[,classes][c(lev1,lev2)])
			#invertir el FC class1/class2
			obj@diffexp.mRNA<-data.frame(FC=logratio2foldchange(as.vector(-RP.out$AveFC)),
				logratio=as.vector(-RP.out$AveFC),
				pval=apply(RP.out$pval,1,min),
				adj.pval=apply(RP.out$pfp,1,min))
			rownames(obj@diffexp.mRNA)<-rownames(mrdat)
		}

		if (method.dif=="DESeq") {
			conds <- obj@pheno.mRNA[,classes][c(lev1,lev2)]
			cds<-newCountDataSet(round(mrdat,0),conds)

			cds<-estimateSizeFactors(cds)
			cds<-estimateDispersions(cds,fitType="local")
			comparative <- nbinomTest (cds, "0", "1")

			obj@diffexp.mRNA<-data.frame(FC=comparative$foldChange,
				logratio=comparative$log2FoldChange,
				pval=comparative$pval,
				adj.pval=comparative$padj)
			rownames(obj@diffexp.mRNA)<-comparative$id
			
			obj@info[["mRNA.weights"]]<-cds@phenoData@data$sizeFactor
		}
		
		
		if (method.dif=="edgeR") {
		  conds <- obj@pheno.mRNA[,classes][c(lev1,lev2)]
		  d <- DGEList(counts=mirdat,group=factor(conds))	
		  design.mat <- model.matrix(~ 0 + d$samples$group)
		  
		  d <- calcNormFactors(d)
		  d1 <- estimateCommonDisp(d, verbose=T)
		  d1 <- estimateTagwiseDisp(d1)
		  
		  colnames(design.mat) <- levels(d$samples$group)
		  d2 <- estimateGLMCommonDisp(d,design.mat)
		  d2 <- estimateGLMTrendedDisp(d2,design.mat, method="auto")
		  # You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
		  # The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
		  d2 <- estimateGLMTagwiseDisp(d2,design.mat)
		  
		  et12 <- exactTest(d1, pair=c("0","1")) # compare groups 1 and 2
		  
		  obj@diffexp.mRNA<-data.frame(FC=logratio2foldchange(et12$table$logFC),
		                                logratio=et12$table$logFC,
		                                pval=et12$table$PValue,
		                                adj.pval=p.adjust(et12$table$PValue,method="BH"))
		  rownames(obj@diffexp.mRNA)<-rownames(et12$table)
		  
		  obj@info[["mRNA.weights"]]<-d@.Data[[2]]$norm.factors
		}
		

		if (method.dif=="anova") {
			pval.aov<-apply(mrdat,1,aov_pval,as.factor(obj@pheno.mRNA[,classes]))
			ss.aov<-apply(mrdat,1,aov_ss,as.factor(obj@pheno.mRNA[,classes]))
			obj@diffexp.mRNA<-data.frame(ss.prop=ss.aov,
				pval=pval.aov,
				adj.pval=p.adjust(pval.aov,method=method.adj))
			rownames(obj@diffexp.mRNA)<-rownames(mrdat)
		}

		if (method.dif=="max.var") {
			myvar <- apply(mrdat,1,var)
			conv.pvals<- abs((myvar-max(myvar))/max(myvar))

			obj@diffexp.mRNA<-data.frame(variance=myvar,
				pval=conv.pvals,
				adj.pval=p.adjust(conv.pvals,method=method.adj))
			rownames(obj@diffexp.mRNA)<-rownames(mrdat)

		}


		obj@diffexp.mRNA$meanExp<-apply(mrdat,1,mean)
		if (method.dif != "anova" & method.dif != "max.var" & method.dif!="only.fc") {
			obj@diffexp.mRNA<-obj@diffexp.mRNA[,c		("FC","logratio","meanExp","pval","adj.pval")]
		}
		if (method.dif=="anova") {
			obj@diffexp.mRNA<-obj@diffexp.mRNA[,c		("ss.prop","meanExp","pval","adj.pval")]
		}
		if (method.dif=="max.var") {
			obj@diffexp.mRNA<-obj@diffexp.mRNA[,c		("variance","meanExp","pval","adj.pval")]
		}
		if (method.dif=="only.fc") {
			obj@diffexp.mRNA<-obj@diffexp.mRNA[,c("FC","logratio","meanExp")]
		}
		obj@info[["mRNA.diffexp.method"]][1]<-method.dif
		obj@info[["mRNA.diffexp.method"]][2]<-classes
	}

	return(obj)
}






addLong <- function (obj, dataset, classes, method.dif="time.point", method.adj="BH", var.t.test=FALSE) {

#function similar to addDiffexp, but for

	if (method.dif=="time.point") {
		obj<-addDiffexp(obj, dataset=dataset, classes=classes, method.dif="t.test",method.adj=method.adj, var.t.test=var.t.test)

		if (dataset == "miRNA") {
			obj@info[["miRNA.diffexp.method"]][1]<-method.dif
		}
		if (dataset == "mRNA") {
			obj@info[["mRNA.diffexp.method"]][1]<-method.dif
		}

		return(obj)
	}


	if (method.dif=="linear.regression") {

		if (dataset == "miRNA") {
			dat<-obj@dat.miRNA
			time<-obj@pheno.miRNA[,classes]	
		}

		if (dataset == "mRNA") {
			dat<-obj@dat.mRNA
			time<-obj@pheno.mRNA[,classes]	
		}

		comp.long<- function (x, t) {
			mod<-summary(lm(x~t))
			return(list=c(mod$coefficients[2,1],mod$coefficients[2,4]))
		}

		slopes<-t(apply(dat,1,comp.long,time))

		final<-data.frame(slope=slopes[,1], meanExp=apply(dat,1,mean), pval=slopes[,2],adj.pval=p.adjust(slopes[,2],method=method.adj))

		if (dataset == "miRNA") {
			obj@diffexp.miRNA<-final
			obj@info[["miRNA.diffexp.method"]][1]<-method.dif
			return(obj)
		}

		if (dataset == "mRNA") {
			obj@diffexp.mRNA<-final
			obj@info[["miRNA.diffexp.method"]][1]<-method.dif
			return(obj)
		}

	}

}










selSubsetExprs <- function (obj, dataset, FC=NA, logratio=foldchange2logratio(FC), slope=NA, pval=NA, adj.pval=NA, min.meanExp=NA, up=FALSE, dw=FALSE) {

	if (dataset == "miRNA") {
		subset<-obj@diffexp.miRNA
	}

	if (dataset == "mRNA") {
		subset<-obj@diffexp.mRNA
	}

	if (!is.na(logratio)) {
		sel.log<-which( abs(subset$logratio) >= abs(logratio) & is.finite(subset$logratio)==TRUE )	
		subset<-subset[sel.log,]
	}

	if (!is.na(slope)) {
		sel.log<-which( abs(subset[,1]) >= abs(slope) & is.finite(subset[,1])==TRUE )	
		subset<-subset[sel.log,]
		}

	if (!is.na(pval)) {
		sel.pval<-which( subset$pval <= pval)	
		subset<-subset[sel.pval,]		
	}
	if (!is.na(adj.pval)) {
		sel.adj.pval<-which( subset$adj.pval <= adj.pval)	
		subset<-subset[sel.adj.pval,]		
	}

	if (!is.na(adj.pval)) {
		sel.adj.pval<-which( subset$adj.pval <= adj.pval)	
		subset<-subset[sel.adj.pval,]		
	}

	if (!is.na(min.meanExp)) {
		sel.min.meanExp<-which(subset$meanExp >= min.meanExp)	
		subset<-subset[sel.min.meanExp,]		
	}
	
	if (up==TRUE & dw==FALSE) {
		if (is.na(slope)) {
			sel.up<-which(subset$logratio>0)
			subset<-subset[sel.up,]
		} else {
			sel.up<-which(subset[,1]>0)
			subset<-subset[sel.up,]
		}
	}
	if (up==FALSE & dw==TRUE) {
		if (is.na(slope)) {
			sel.dw<-which(subset$logratio<0)
			subset<-subset[sel.dw,]
		} else {
			sel.dw<-which(subset[,1]<0)
			subset<-subset[sel.dw,]
		}
	}

	return(subset)

}


selSubsetCor <- function (obj, pval.cutoff=1, dat.sum=0, sub.miRNA=NULL, sub.mRNA=NULL) {
	sel<-1:nrow(obj@net)
	if (!is.null(obj@net$adj.pval)==TRUE) {
		sel.p<-which(obj@net$adj.pval<=pval.cutoff)
		sel<-intersect(sel,sel.p)
	}
	if (!is.null(obj@net$dat.sum)==TRUE) {
		sel.ds<-which(obj@net$dat.sum>=dat.sum)
		sel<-intersect(sel,sel.ds)
	}
	if (!is.null(sub.miRNA)==TRUE) {
		sel.miRNA<-which(obj@net$miRNA %in% sub.miRNA)
		sel<-intersect(sel,sel.miRNA)
	}	
	if (!is.null(sub.mRNA)==TRUE) {
		sel.mRNA<-which(obj@net$mRNA %in% sub.mRNA)
		sel<-intersect(sel,sel.mRNA)
	}	
	return(obj@net[sel,])
}


#net <- function (obj) {
#	return(obj@net)
#}

addDatabase <- function (obj, database, pval.ref=1, dat.sum=1) {
	obj@info[["pcomb.method"]]<-NULL
	obj@info[["padjust.method"]]<-NULL
	obj@info[["dat.sum"]]<-dat.sum

	cat("Intersecting with database\n")

	if (database[1]=="microCosm_v5_18_numeric") {

		cat(" microCosm_v5_18 database chosen\n")
		a<-intersect(rownames(obj@net),rownames(microCosm_v5_18))
		#subsubsetmicro<-microCosm_unic[a,]
		obj@net$pval.database<-1	
		obj@net[a,"pval.database"]<-microCosm_v5_18[a,"pval"]
		obj@info[["database"]]<-paste(database,"_numeric",sep="")
	}

	else {

	for (i in database) {
		cat(paste(" ",i," database chosen\n",sep=""))
		sel<-intersect(rownames(obj@net),rownames(get(i)))
		#subsubsetmicro<-microCosm_unic[a,]
		name<-paste("dat.",i,sep="")
		obj@net[,name]<-0
		
		obj@net[sel,name]<-1

	}
	obj@info[["database"]]<-database
	if (length(database)>1) {
		obj@net[,"dat.sum"]<-apply(obj@net[,grep("^dat.",colnames(obj@net))],1,sum)
		} else {
		if (database!="microCosm_v5_18_numeric") {
			obj@net[,"dat.sum"]<-obj@net[,grep("^dat.",colnames(obj@net))]
		}

}
	
	
	}

	return(obj)
}



addSig <- function (obj, dataset, FC=NA, logratio=foldchange2logratio(FC), slope=NA, pval=NA, adj.pval=NA, min.meanExp=NA, up=FALSE, dw=FALSE, manual=NULL) {

	#cleaning of previous information
	obj@cor<-matrix()
	obj@pval<-matrix()
	obj@net<-data.frame()
	obj@info[["pcomb.method"]]<-NULL
	obj@info[["padjust.method"]]<-NULL
	obj@info[["correlation.samples.used"]]<-NULL
	obj@info[["correlation.function.used"]]<-NULL
	obj@info[["correlation.type"]]<-NULL

	if (is.null(manual)) {
		llista<-rownames(selSubsetExprs(obj,dataset,FC=FC,logratio=logratio, slope=slope, pval=pval, adj.pval=adj.pval, min.meanExp=min.meanExp, up=up, dw=dw))

		criteria<-vector()
		i<-1
		if (!is.na(logratio)) {criteria[i]<-paste("abs(logratio) >",round(logratio,2));i<-i+1}
		if (!is.na(slope)) {criteria[i]<-paste("abs(slope)",slope);i<-i+1}
		if (!is.na(FC)) {criteria[i]<-paste("abs(FC) >",FC);i<-i+1}
		if (!is.na(pval)) {criteria[i]<-paste("pval <",pval);i<-i+1}
		if (!is.na(adj.pval)) {criteria[i]<-paste("adj.pval <",adj.pval);i<-i+1}
		if (!is.na(min.meanExp)) {criteria[i]<-paste("min.meanExp >",min.meanExp);i<-i+1}
		if (up) {criteria[i]<-c("filter for UP");i<-i+1}
		if (dw) {criteria[i]<-c("filter for DOWN");i<-i+1}

		if (dataset == "miRNA") {
			obj@sig.miRNA<-llista
			obj@info[["miRNA.criteria"]]<-criteria
		}
	
		if (dataset == "mRNA") {
			obj@sig.mRNA<-llista
			obj@info[["mRNA.criteria"]]<-criteria
		}
	}
	else {
		llista<-manual
		if (dataset == "miRNA") {
			obj@sig.miRNA<-llista
			obj@info[["miRNA.criteria"]]<-"manual!"
		}

		if (dataset == "mRNA") {
			obj@sig.mRNA<-llista
			obj@info[["mRNA.criteria"]]<-"manual!"
		}
	}
	
	return(obj)
}






addScore <- function (obj) {
	compute.score<-function(x,y) {
		res<-(( x*cos(pi/4)-y*sin(pi/4) ) ^2 -
	      ( x*sin(pi/4)+y*cos(pi/4) ) ^2 )
		return(res)
	}

	if (obj@info[["miRNA.diffexp.method"]][1]=="time.point" | obj@info[["miRNA.diffexp.method"]][1]=="linear.regression" | obj@info[["miRNA.diffexp.method"]][1]=="anova") {
		if (obj@info[["miRNA.diffexp.method"]][1]=="anova") {
			obj@net$score <- obj@net$ss.prop.miRNA * obj@net$ss.prop.mRNA
		} else {
			obj@net$score <- compute.score (obj@net$slope.miRNA,obj@net$slope.mRNA)            
		}                    
	} else {
		obj@net$score <- compute.score (obj@net$logratio.miRNA,obj@net$logratio.mRNA)
	}

	return(obj)
}



combinePval <- function (obj, pval.1 = "pval", pval.2 = "pval.database", method="stouffer", w=c(1,1)) {
	obj@info[["padjust.method"]]<-NULL

	cat("Combining p.values\n")
	pval.1<-obj@net[,pval.1]
	pval.2<-obj@net[,pval.2]

	if (method=="stouffer") {
		p.comb<-1-pnorm(1/sqrt(w[1]^2+w[2]^2)*(w[1]*qnorm(1-pval.1)+w[2]*qnorm(1-pval.2)))
		}

	if (method=="fisher") {
		fishert<-(-2*log(apply(data.frame(pval.1,pval.2),1,prod)))
		p.comb<-1-pchisq(fishert,4)
		}
	obj@net$p.comb<-p.comb
	obj@info[["pcomb.method"]]<-method
	return(obj)
}


#setMethod("correctPval", "corObject", correctPval) #això no xuta


correctPval <- function (obj, method.adj="BH", pval="pval") {
	cat("Correcting p.values\n")
	obj@net$adj.pval<-p.adjust(obj@net[,pval],method=method.adj)
	obj@info[["padjust.method"]]<-method.adj
	return(obj)
}



##### evaluate databases

evaluate <- function ( obj, method=c("hypergeometric", "logistic", "GSEA"), databases="all", adj.pval=0.05 , plot=TRUE, miRNAs="all" , mRNAs= "all" , nperm=nperm ) {
    
  if (databases=="all") {
    test.dat<-colnames(obj@net)[grep("dat.",colnames(obj@net))]
  } else {
    if (!all(paste("dat.",databases,sep="") %in% colnames(obj@net)[grep("dat.",colnames(obj@net))])) stop("Selected databases not in the corObject")
    test.dat<-c(paste("dat.",databases,sep=""),"dat.sum")
  }
  
  #subnet <- obj@net[which(obj@net$miRNA %in% miRNAs & obj@net$mRNA %in% mRNAs) , ]
  subnet <- obj@net
  
  result<-data.frame(miRNA=as.character(unique(subnet$miRNA)))
  
  result$tot.targets <- table(subnet$miRNA[which(subnet$dat.sum>0 & subnet$adj.pval < adj.pval)])[result$miRNA]
  
  
  if ("hypergeometric" %in% method | "hyper" %in% method) {
    
    cat("Hypergeometric testing\n")
    
    for (i in 1:length(test.dat)) {
      cat(paste("  Testing database",test.dat[i],"\n"))
      k<-dim(result)[2]
      
      targets<-table(subnet$miRNA[which(subnet$adj.pval<adj.pval & subnet[,test.dat[i]]>0)])[result$miRNA]
      cor<-table(subnet$miRNA[which(subnet$adj.pval<adj.pval)])[result$miRNA]
      pred.targets<-table(subnet$miRNA[which(subnet[,test.dat[i]]>0)])[result$miRNA]
      size<-table(subnet$miRNA)[result$miRNA]
      
      #saving the results
      result[,k+1]<-targets#[result$miRNA]
      result[,k+2]<-phyper(targets-1, pred.targets, size-pred.targets, cor, lower.tail=FALSE )
      result[,k+3]<-p.adjust(result[,k+2],method="BH")
      result[,k+4]<-(size-cor) * targets / ((pred.targets-targets) * cor) 
       
      #### do fisher test
      allt<-table(subnet[,test.dat[i]]>0,subnet$adj.pval<0.05,subnet$miRNA)
      result[,k+5]<-NA
      for (ii in 1:length(result$miRNA)) {
        result[ii,k+5]<-chisq.test(allt[,,result$miRNA[ii]])$p.value
      }
      result[,k+6]<-p.adjust(result[,k+5],method="BH")
          
      #changing column names
      colnames(result)[(k+1):(k+6)] <- paste("h",test.dat[i],c("targets","pval","FDR","OR","fish.p","fish.FDR"),sep=".")
      
    }
    
    if (plot==TRUE) {
      x11()
      boxplot(as.matrix(log2(result[,grep("h*OR",colnames(result))])))
    } 
    
  }
   
  
  if ("logistic" %in% method | "log" %in% method | "regression" %in% method) {
    
    cat("Logistic regression\n")   
    #library(pROC)
    #library(verification)
    
    for (i in 1:length(test.dat)) {
      k<-dim(result)[2]
      
      result[,(k+1):(k+6)]<-NA
      
      colnames(result)[(k+1):(k+6)]<-paste("l",test.dat[i],c("beta1","pval","FDR","auc","pval.a","FDR.a"),sep=".")
      
      for (j in 1:dim(result)[1]) {
        cat(paste("  Testing database",test.dat[i],"miRNA",result[j,"miRNA"],"\n"))
        
        resp <- as.numeric(subnet[which(subnet$miRNA == result[j,"miRNA"]),test.dat[i]]>0)
        
        if (length(table(resp))>1) {
          cor <- subnet$cor[which(subnet$miRNA == result[j,"miRNA"])] 
          
          #logistic regression
          mod <- glm(resp ~ cor, family=binomial)
          prob<-predict(mod,type=c("response"))
          g1 <- roc(resp ~ prob)
          
          result[j,k+1]<-summary(mod)$coefficients[2,1]
          result[j,k+2]<-summary(mod)$coefficients[2,4]
          
          result[j,k+4]<-g1$auc
          result[j,k+5]<-roc.area(resp,prob)$p.value
          
        }
        
        result[,k+3]<-p.adjust(result[,k+2],method="BH")
        result[,k+6]<-p.adjust(result[,k+5],method="BH")
        
      }
    } 
  }
  
  
  if ("GSEA" %in% method | "gsea" %in% method) {
    cat("GSEA\n")
    
    #library(fgsea)
    micro<-result$miRNA
    k<-dim(result)[2]
    result[,(k+1):(k+4*length(test.dat))]<-NA
    colnames(result)[(k+1):(k+4*length(test.dat))]<-paste("g",rep(test.dat,each=4),c("NES","ES","pval","FDR"),sep=".")
      
    for (j in 1:length(micro)) {  
	    
      cat(paste("  Testing miRNA",micro[j],"\n"))
      
      rank<-subnet$cor[which(subnet$miRNA==micro[j])]
      names(rank)<-subnet$mRNA[which(subnet$miRNA==micro[j])]		
      rank<-sort(rank)
      
      pathways<-list()
      for (i in 1:length(test.dat)) {
        pathways[[i]] <- as.character(subnet$mRNA[which(subnet$miRNA==micro[j] & subnet[,test.dat[i]]>0)])
      }
      names(pathways)<-test.dat
      
      fgseaRes <- fgsea(pathways = pathways, 
                        stats = rank,
                        minSize=1,
                        maxSize=5000,
                        nperm=nperm)
      
      res<-data.frame(fgseaRes)
      
      for (i in 1:length(test.dat)) {
	      
        if (test.dat[i] %in% res[,1]) {
          sel<-which(res[,1] %in% test.dat[i])
          result[j,k+4*(i-1)+1]<-res[sel,"NES"]
          result[j,k+4*(i-1)+2]<-res[sel,"ES"]
          result[j,k+4*(i-1)+3]<-res[sel,"pval"]
          result[j,k+4*(i-1)+4]<-p.adjust(res[sel,"pval"],method="BH") 
        } 
      } 
    } 
  }    
  return(result)
}






plotHeatmap <- function(obj, class, n=50, col.color=1, min.exp=NULL, main=NULL, pval.cutoff=NULL, grouping.col=NULL, grouping.row=NULL, order="pval", cex.lab=0.75) {
	#library(pheatmap)
#	plot.new()


	triming<-function(X){
		num<-X
		for(i in 1:dim(num)[1]) {
			x<-num[i,]
        		trim = 0.05
        		lo = quantile(x, trim)
        		hi = quantile(x, 1 - trim)
        		x[x < lo] = lo
        		x[x > hi] = hi
			num[i,]<-x
		}
		return(num)
	}

if (length(col.color)==1) {

	nf<-layout(mat=matrix(c(1),ncol=1,nrow=1))
	par(mar=c(5.1, 4.1, 4.1, 2.1))

	colors.cat<-rep(palette()[-1],length.out=nrow(obj@pheno.miRNA))

	#colors.cat<-c("blue","red","green","violet","orange","pink","yellow")

	col<-c("turquoise","violet")

	if (class=="miRNA") {
		#cc<-as.factor(obj@pheno.miRNA[,col.color])


		if (length(levels(as.factor(obj@pheno.miRNA[,col.color])))==2) {
			cc.color<-as.character(factor(as.numeric(as.factor(obj@pheno.miRNA[,col.color])),labels=col))
		} else {
			cc.color<-as.character(factor(as.numeric(as.factor(obj@pheno.miRNA[,col.color]))+1))
		}



		#cc.color<-colors.cat[as.numeric(cc)]
		
		data<-obj@diffexp.miRNA
		if (!is.null(min.exp)) {
			data<-data[which(data$meanExp>=min.exp),]
		}

		if (order=="FDR" | order=="pval") {
			sel.sort<-rownames(data[with(data,order(pval)),])
		} else {
			if (order=="fc" | order=="logratio") {
				sel.sort<-rownames(data[order(abs(data$logratio),decreasing=TRUE),])
			} else {
				stop("no matching order")
			}
		}
		
		
		

		xmat<-obj@dat.miRNA[sel.sort[1:n],]
		rownames(xmat)<-gsub("hsa-","",rownames(xmat))
		xmat<-xmat[which(apply(xmat,1,sd)>0),]

		cc.row<-hclust(dist(xmat,method="euclidean"))
		cc.col<-hclust(dist(t(xmat),method="euclidean"))

		if (is.null(main)) {main <- paste("Top",n,class)}
		

		if (!is.null(grouping.row)) {

			data<-obj@dat.miRNA

			lev.order<-names(table(grouping.row))
			for (i in 1:length(lev.order)) {

				distdat<-hclust(dist(data[which(grouping.row==lev.order[i]),]))
				if (i==1) {
					new.order<-data[which(grouping.row==lev.order[i])[distdat$order],]
					rownames(new.order)<-rownames(data[which(grouping.row==lev.order[i]),])[distdat$order]
				} else {
					lines<-nrow(new.order)
					new.order<-rbind(new.order,data[which(grouping.row==lev.order[i])[distdat$order],])
					print(i)
					print(new.order)

					rownames(new.order)[(lines+1):nrow(new.order)]<-rownames(data[which(grouping.row==lev.order[i]),])[distdat$order]

				}		
				

			}

			#heatmap.2(triming(new.order), col=greenred(75), scale="row", ColSideColors=cc.color, RowSideColors=as.character(rep(1:length(lev.order),table(grouping.row))), key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.75, main=main,labCol=NA,dendrogram="col",distfun=function(x) as.dist(1-cor(t(x), method="pearson")),hclustfun=function(x) hclust(x,method="average"))
			print("entered")
			heatmap.2(triming(new.order), col=greenred(75), scale="row", ColSideColors=cc.color, RowSideColors=as.character(rep(1:length(lev.order),table(grouping.row))), Rowv=FALSE, Colv=TRUE, key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cex.lab, main=main,labCol=NA,dendrogram="col")
			leg<-levels(as.factor(obj@pheno.miRNA[,col.color]))
			legend("topright",leg,col=levels(as.factor(cc.color)),lwd=7)




		} else {

		#	heatmap.2(triming(xmat)[cc.row$order,cc.col$order], col=redgreen(75), scale="row", ColSideColors=cc.color[cc.col], Colv=FALSE, Rowv=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.75, main=paste("Top",n,class),labCol=NA,dendrogram="none")
			heatmap.2(triming(xmat), col=greenred(75), scale="row", ColSideColors=cc.color, key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cex.lab, main=main,labCol=NA,dendrogram="both",distfun=function(x) as.dist(1-cor(t(x), method="pearson")),hclustfun=function(x) hclust(x,method="average"), margins=c(10,10))
			leg<-levels(as.factor(obj@pheno.miRNA[,col.color]))
			legend("topright",leg,col=levels(as.factor(cc.color)),lwd=7)

		}

	}


	if (class=="mRNA") {
		#cc<-as.factor(obj@pheno.mRNA[,col.color])
		#cc.color<-colors.cat[as.numeric(cc)]
		

		if (length(levels(as.factor(obj@pheno.mRNA[,col.color])))==2) {
			cc.color<-as.character(factor(as.numeric(as.factor(obj@pheno.mRNA[,col.color])),labels=col))
		} else {
			cc.color<-as.character(factor(as.numeric(as.factor(obj@pheno.mRNA[,col.color]))+1))
		}


		data<-obj@diffexp.mRNA
		if (!is.null(min.exp)) {
			data<-data[which(data$meanExp>=min.exp),]
		}

		if (order=="FDR" | order=="pval") {
			sel.sort<-rownames(data[with(data,order(pval)),])
		} else {
			if (order=="fc" | order=="logratio") {
				sel.sort<-rownames(data[order(abs(data$logratio),decreasing=TRUE),])
			} else {
				stop("no matching order")
			}
		}


		xmat<-obj@dat.mRNA[sel.sort[1:n],]
		xmat<-xmat[which(apply(xmat,1,sd)>0),]

		cc.row<-hclust(dist(xmat,method="euclidean"))$order
		cc.col<-hclust(dist(t(xmat),method="euclidean"))$order

		if (is.null(main)) {main <- paste("Top",n,class)}





		if (!is.null(grouping.row)) {

			data<-obj@dat.mRNA

			lev.order<-names(table(grouping.row))
			for (i in 1:length(lev.order)) {

				distdat<-hclust(dist(data[which(grouping.row==lev.order[i]),]))
				print(distdat$order)
				if (i==1) {
					new.order<-data[which(grouping.row==lev.order[i])[distdat$order],]
					rownames(new.order)<-rownames(data[which(grouping.row==lev.order[i]),])[distdat$order]
				} else {
					lines<-nrow(new.order)
					new.order<-rbind(new.order,data[which(grouping.row==lev.order[i])[distdat$order],])
					rownames(new.order)[(lines+1):nrow(new.order)]<-rownames(data[which(grouping.row==lev.order[i]),])[distdat$order]

				}		
				

			}

			#heatmap.2(triming(new.order), col=greenred(75), scale="row", ColSideColors=cc.color, RowSideColors=as.character(rep(1:length(lev.order),table(grouping.row))), key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.75, main=main,labCol=NA,dendrogram="col",distfun=function(x) as.dist(1-cor(t(x), method="pearson")),hclustfun=function(x) hclust(x,method="average"))
			print("entered")
			heatmap.2(triming(new.order), col=greenred(75), scale="row", ColSideColors=cc.color, RowSideColors=as.character(rep(1:length(lev.order),table(grouping.row))), Rowv=FALSE, Colv=TRUE, key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cex.lab, main=main,labCol=NA,dendrogram="col")
			leg<-levels(as.factor(obj@pheno.mRNA[,col.color]))
			legend("topright",leg,col=levels(as.factor(cc.color)),lwd=7)


		} else {

		#heatmap.2(triming(xmat)[cc.row,cc.col], col=redgreen(75), scale="row", ColSideColors=cc.color[cc.col], Colv=FALSE, Rowv=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.75, main=paste("Top",n,class),labCol=NA)
			heatmap.2(triming(xmat), col=greenred(75), scale="row", ColSideColors=cc.color, key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cex.lab, main=main,labCol=NA,dendrogram="both",distfun=function(x) as.dist(1-cor(t(x), method="pearson")),hclustfun=function(x) hclust(x,method="average") , margins=c(10,10))

			leg<-levels(as.factor(obj@pheno.mRNA[,col.color]))

			legend("topright",leg,col=levels(as.factor(cc.color)),lwd=7)
		}

	}


	if (class=="both") {
		cyto<-obj@net[with(obj@net,order(adj.pval)),]
		names.miRNA<-unique(cyto[1:n,"miRNA"])
		names.mRNA<-unique(cyto[1:n,"mRNA"])

		xmat.miRNA<-obj@dat.miRNA[names.miRNA,]
		rownames(xmat.miRNA)<-gsub("hsa-","",rownames(xmat.miRNA))

		xmat.mRNA<-obj@dat.mRNA[names.mRNA,]
		xmat<-rbind(xmat.miRNA,xmat.mRNA)
		xmat<-xmat[which(apply(xmat,1,sd)>0),]

		#cc<-as.factor(obj@pheno.mRNA[,col.color])
		#cc.color<-colors.cat[as.numeric(cc)]



		if (length(levels(as.factor(obj@pheno.miRNA[,col.color])))==2) {
			cc.color<-as.character(factor(as.numeric(as.factor(obj@pheno.miRNA[,col.color])),labels=col))
		} else {
			cc.color<-as.character(factor(as.numeric(as.factor(obj@pheno.miRNA[,col.color]))+1))
		}


		if (is.null(main)) {main <- paste("Top",n,class)}


		cc.row<-hclust(dist(xmat,method="euclidean"))
		cc.col<-hclust(dist(t(xmat),method="euclidean"))

		#heatmap.2(triming(xmat)[cc.row$order,cc.col$order], col=redgreen(75), scale="row", ColSideColors=cc.color[cc.col$order], Colv=FALSE, Rowv=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.75, main=paste("Top",n,class),labCol=NA)
		heatmap.2(triming(xmat), col=greenred(75), scale="row", ColSideColors=cc.color, key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cex.lab, main=main,labCol=NA,dendrogram="both",distfun=function(x) as.dist(1-cor(t(x), method="pearson")),hclustfun=function(x) hclust(x,method="average"))

		leg<-levels(as.factor(obj@pheno.mRNA[,col.color]))
		legend("topright",leg,col=levels(as.factor(cc.color)),lwd=7)


	}




}


if (length(col.color)>1) {


	if (class=="miRNA") {
		data<-obj@diffexp.miRNA
		if (!is.null(min.exp)) {
			data<-data[which(data$meanExp>=min.exp),]
		}

		sel.sort<-rownames(data[with(data,order(pval)),])
		xmat<-obj@dat.miRNA[sel.sort[1:n],]
		rownames(xmat)<-gsub("hsa-","",rownames(xmat))

		xmat<-triming(xmat)
		xmat<-xmat[which(apply(xmat,1,sd)>0),]

		if (is.null(main)) {main <- paste("Top",n,class)}

		pheatmap(xmat, 	annotation_col = data.frame(obj@pheno.miRNA[,col.color]), color=greenred(75), main=main, scale="row")
	}


	if (class=="mRNA") {
		data<-obj@diffexp.mRNA
		if (!is.null(min.exp)) {
			data<-data[which(data$meanExp>=min.exp),]
		}

		sel.sort<-rownames(data[with(data,order(pval)),])
		xmat<-obj@dat.mRNA[sel.sort[1:n],]

		xmat<-triming(xmat)
		xmat<-xmat[which(apply(xmat,1,sd)>0),]

		if (is.null(main)) {main <- paste("Top",n,class)}

		pheatmap(xmat, 	annotation_col = data.frame(obj@pheno.mRNA[,col.color]), color=greenred(75), main=main, scale="row")
	}


	if (class=="both") {
		cyto<-obj@net[with(obj@net,order(adj.pval)),]
		names.miRNA<-unique(cyto[1:n,"miRNA"])
		names.mRNA<-unique(cyto[1:n,"mRNA"])

		xmat.miRNA<-obj@dat.miRNA[names.miRNA,]
		rownames(xmat.miRNA)<-gsub("hsa-","",rownames(xmat.miRNA))

		xmat.mRNA<-obj@dat.mRNA[names.mRNA,]
		xmat<-rbind(xmat.miRNA,xmat.mRNA)

		xmat<-triming(xmat)
		xmat<-xmat[which(apply(xmat,1,sd)>0),]

		if (is.null(main)) {main <- paste("Top",n,class)}

		pheatmap(xmat, 	annotation_col = data.frame(obj@pheno.miRNA[,col.color]), color=greenred(75), main=main, scale="row")

	}

}

}



plotVolcano <- function (obj, subset, FC1=1.5, FC2=2, FDR=0.05, cex=1, cex.lab=1, cex.axis=1) {

	nf<-layout(mat=matrix(c(1),ncol=1,nrow=1))
	par(mar=c(5.1, 4.1, 4.1, 2.1))
#	plot.new()


	if (subset=="miRNA") {data<-obj@diffexp.miRNA}
	if (subset=="mRNA") {data<-obj@diffexp.mRNA}

	plot(data$logratio,-log10(data$adj.pval),pch=19,xlab="log2ratio",ylab="-log10(FDR)",cex=cex, cex.lab=cex.lab, cex.axis=cex.axis)
	sel.1<-which(data$adj.pval <= FDR)
	points(data$logratio[sel.1],-log10(data$adj.pval[sel.1]),pch=19,col="yellow",cex=cex)

	sel.2<-which(data$adj.pval <= FDR & abs(data$logratio) > log2(FC1))
	points(data$logratio[sel.2],-log10(data$adj.pval[sel.2]),pch=19,col="orange",cex=cex)

	sel.3<-which(data$adj.pval <= FDR & abs(data$logratio) > log2(FC2))
	points(data$logratio[sel.3],-log10(data$adj.pval[sel.3]),pch=19,col="red",cex=cex)
	
	abline(h=-log10(FDR),lty=3)
	abline(v=c(log2(FC1),-log2(FC1)),lty=2)
	abline(v=c(log2(FC2),-log2(FC2)),lty=4)

	legend("topright", c(
		paste("FDR <",FDR),
		paste("FDR <",FDR,"& abs(FC) >",FC1),
		paste("FDR <",FDR,"& abs(FC) >",FC2)),
		col=c("yellow","orange","red"),pch=19)
}


plotMA <- function (obj, subset, sample1=NULL, sample2=NULL, cex.lab=1, cex.axis=1 ) {

	if (subset=="miRNA") {
		plot(apply(obj@dat.miRNA[,c(sample1,sample2)],1,mean),obj@dat.miRNA[,sample1]-obj@dat.miRNA[,sample2],pch=19,cex=0.5, xlab="Mean log2(intensity)", ylab="log2 ratio", cex.lab=cex.lab, cex.axis=cex.axis)
	}

	if (subset=="mRNA") {
		plot(apply(obj@dat.mRNA[,c(sample1,sample2)],1,mean),obj@dat.mRNA[,sample1]-obj@dat.mRNA[,sample2],pch=19,cex=0.5, xlab="Mean log2(intensity)", ylab="log2 ratio", cex.lab=cex.lab, cex.axis=cex.axis)
	}
}





#plotHeatmap(data.obj,"miRNA")
#Error in cc.color[cc.col] : invalid subscript type 'list'


#setMethod("heatmap", "corObject", heatmap_corObject) #això no xuta
#plotHeatmap(data.obj,"both")


#obj<-data.obj.GSE17498
#pval.cutoff=0.05
#sub.miRNA=NULL
#sub.mRNA=NULL
#names=TRUE
#dat.sum=obj@info[["dat.sum"]]
#add.other=NULL
#vertex.cex=NULL
#n=NULL
#node.size=1.5



plotNetwork <- function (obj, pval.cutoff=0.05, score.cutoff=NULL, sub.miRNA=NULL, sub.mRNA=NULL, names=TRUE, dat.sum=obj@info[["dat.sum"]], add.other=NULL, vertex.cex=NULL, n=NULL, node.size=1.5) {

	nf<-layout(mat=matrix(c(1),ncol=1,nrow=1))
	par(mar=c(5.1, 4.1, 4.1, 2.1))
#	plot.new()

#vertex.cex="interact.table.human"

if (is.null(obj@net$score)==TRUE) {
	if ( (!is.null(obj@net$logratio.miRNA) | !is.null(obj@net$slope.miRNA)) & (!is.null(obj@net$logratio.mRNA) | !is.null(obj@net$slope.mRNA) ) ) {
		obj<-addScore(obj)
	} else {
		obj@net$score<-1
	}
}

if (obj@info[["miRNA.diffexp.method"]][1]=="anova") {
  cat("Colours of the plot won't render correctly if you've used anova method for differential expression in the miRNAs.\n")
}	
if (obj@info[["mRNA.diffexp.method"]][1]=="anova") {
  cat("Colours of the plot won't render correctly if you've used anova method for differential expression in the mRNAs.\n")
}	
	
	
#	library(network)
llistacomp<-obj@net

## if there are scores coming from addLong function, treat them as logratios
if ("slope.miRNA" %in% colnames(llistacomp)) {
  colnames(llistacomp)[which(colnames(llistacomp)=="slope.miRNA")]<-"logratio.miRNA"
}
if ("slope.mRNA" %in% colnames(llistacomp)) {
  colnames(llistacomp)[which(colnames(llistacomp)=="slope.mRNA")]<-"logratio.mRNA"
}



if (!is.null(sub.miRNA)) {
	llistacomp<-llistacomp[llistacomp$miRNA %in% sub.miRNA,]
}
if (!is.null(sub.mRNA)) {
	llistacomp<-llistacomp[llistacomp$mRNA %in% sub.mRNA,]
}

sel<-1:nrow(llistacomp)

if (!is.null(dat.sum)) {
	sel1<-which(llistacomp$dat.sum>=dat.sum)
	sel<-intersect(sel,sel1)
}

if (!is.null(obj@net$adj.pval)) {
	sel2<-which(llistacomp$adj.pval<=pval.cutoff)
	sel<-intersect(sel,sel2)
}

if (!is.null(score.cutoff)) {
	sel3<-which(llistacomp$score>=score.cutoff)
	sel<-intersect(sel,sel3)
}


#if (!is.null(dat.sum)) {
#	sel<-which(llistacomp$dat.sum>=dat.sum & llistacomp$adj.pval<=pval.cutoff)
#} else {
#	sel<-which(llistacomp$adj.pval<=pval.cutoff)
#}


if (length(sel)==0) {
	n<-50
	sel<-1:nrow(llistacomp)
	cat(" No miRNA-mRNA pairs meeting the specified criteria, selecting the top 50 miRNA-mRNA pairs according to its correlation value\n")
}

llistacomp<-llistacomp[sel,]

if (!is.null(n)) {
	if (n<nrow(llistacomp)) {
		if (!is.null(score.cutoff)) {
			llistacomp<-llistacomp[with(llistacomp,order(score, decreasing=TRUE)),]
		}
		if (!is.null(obj@net$adj.pval)) {
			llistacomp<-llistacomp[with(llistacomp,order(adj.pval)),]
		}
		llistacomp<-llistacomp[1:n,]
	}
}


### colorejar\n
remap <- function(x) { (( x ) / max( abs(x) ) / 2)+0.5 }  # map x onto [0, 1]
remap.light <- function(x) { (( x ) / max( abs(x) + quantile(abs(x),0.2,na.rm=TRUE) ) / 2)+0.5 }  # map x onto [0, 1]

fun.col.edges <- function(x) {rgb(colorRamp(c("green","gray95", "red"))(remap(x)),
                        maxColorValue = 255)
}
fun.col <- function(x) {rgb(colorRamp(c("green","white", "red"))(remap(x)),
                        maxColorValue = 255)
}


#triming.list
triming.list<-function(x){
 	trim = 0.05
 	lo = quantile(x, trim)
 	hi = quantile(x, 1 - trim)
       	x[x < lo] = lo
       	x[x > hi] = hi
	return(x)
}

### if there is no score information available
if (!is.null(llistacomp$score)) {
	color <- with(llistacomp, fun.col.edges(as.numeric(llistacomp$score)) )
} else {
	cat("  No score defined, stating arrow color to 'red'\n")
	color <- rep(1,nrow(llistacomp))
}


### if there is no logratio information available
if (is.null(llistacomp$logratio.miRNA)) {
	if (obj@info[["miRNA.diffexp.method"]][1]=="anova") {
		llistacomp$logratio.miRNA<-llistacomp$ss.prop.miRNA	
	}
	if (obj@info[["miRNA.diffexp.method"]][1]=="linear.regression") {
		llistacomp$logratio.miRNA<-llistacomp$slope.miRNA	
	}
	if (obj@info[["miRNA.diffexp.method"]][1]!="linear.regression" & obj@info[["miRNA.diffexp.method"]][1]!="anova") {
		llistacomp$logratio.miRNA<-0	
	}
}
if (is.null(llistacomp$logratio.mRNA)) {
	if (obj@info[["mRNA.diffexp.method"]][1]=="anova") {
		llistacomp$logratio.mRNA<-llistacomp$ss.prop.mRNA	
	}
	if (obj@info[["mRNA.diffexp.method"]][1]=="linear.regression") {
		llistacomp$logratio.mRNA<-llistacomp$slope.mRNA	
	}
	if (obj@info[["mRNA.diffexp.method"]][1]!="linear.regression" & obj@info[["mRNA.diffexp.method"]][1]!="anova") {
		llistacomp$logratio.mRNA<-0	
	}
}

### if there is no p-value
if (!is.null(obj@net$adj.pval)) {
	pvals.transform<-log(llistacomp$adj.pval)
	sel<-which(pvals.transform=="-Inf")
	if (length(sel)>0) {
		pvals.transform[sel]<-min(pvals.transform[-sel])-1
		edge.width <- remap.light(-pvals.transform)
	} else {
		pvals.transform<-log(0.5)
		edge.width <- remap.light(-pvals.transform)
	}
}

if (!is.null(dat.sum)) {
	edge.width<-llistacomp$dat.sum/max(llistacomp$dat.sum)*3
}



#color<-remap(llistacomp$score)
#edge.width <- with (llistacomp, fun.col(as.numeric(-log(llistacomp$adj.pval))))

cadena<-c(as.character(llistacomp[,"miRNA"]),as.character(llistacomp[,"mRNA"]))
cadenau<-unique(cadena)
nodeinfo<-c(rep(100,length(unique(llistacomp[,"miRNA"]))),rep(4,length(unique(llistacomp[,"mRNA"]))))
noderot<-c(rep(0,length(unique(llistacomp[,"miRNA"]))),rep(45,length(unique(llistacomp[,"mRNA"]))))

taula<-matrix(0,ncol=length(cadenau),nrow=length(cadenau))
rownames(taula)<-as.character(cadenau)
colnames(taula)<-as.character(cadenau)

taula.adj<-matrix(0,ncol=length(cadenau),nrow=length(cadenau))
rownames(taula.adj)<-as.character(cadenau)
colnames(taula.adj)<-as.character(cadenau)

taula.edge<-matrix(0,ncol=length(cadenau),nrow=length(cadenau))
rownames(taula.edge)<-as.character(cadenau)
colnames(taula.edge)<-as.character(cadenau)



for (i in 1:nrow(llistacomp))	
	{
	taula[as.character(llistacomp[i,"miRNA"]),as.character(llistacomp[i,"mRNA"])]<-1
	taula.adj[as.character(llistacomp[i,"miRNA"]),as.character(llistacomp[i,"mRNA"])]<-color[i]
	taula.edge[as.character(llistacomp[i,"miRNA"]),as.character(llistacomp[i,"mRNA"])]<-edge.width[i]
	}



if (!is.null(add.other)) {

#data(interact)
interact<-get(add.other)

add<-interact[which((as.character(interact[,1]) %in% as.character(llistacomp$mRNA) & as.character(interact[,3]) %in% as.character(llistacomp$mRNA) )==TRUE) , ] 

add[,1]<-as.character(add[,1])
add[,3]<-as.character(add[,3])


for (i in 1:nrow(add))	
	{
	taula[add[i,1],add[i,3]]<-1
	taula[add[i,3],add[i,1]]<-1

	taula.adj[add[i,1],add[i,3]]<-8
	taula.adj[add[i,3],add[i,1]]<-8

	taula.edge[add[i,1],add[i,3]]<-0.5
	taula.edge[add[i,3],add[i,1]]<-0.5

	}


}

taulanet<-network(taula,directed=TRUE)

#library(gplots)

onlymir<-llistacomp[which(duplicated(as.character(llistacomp[,1]))==FALSE),]
scoresmir<-onlymir[,"logratio.miRNA"]
onlygene<-llistacomp[which(duplicated(as.character(llistacomp[,2]))==FALSE),]
scoresgene<-onlygene[,"logratio.mRNA"]
scorescomp<-c(scoresmir,scoresgene)
scorescomp<-data.frame(names=c(onlymir[,1],onlygene[,2]),scores=scorescomp)

nodeinfo2<-with(scorescomp, fun.col(as.numeric(scorescomp$scores)))

nodetam<-rep(1.5,length(nodeinfo2))
#nodetam<-rep(0.5,length(nodeinfo2))

nodetam<-rep(node.size,length(nodeinfo2))

if (!is.null(vertex.cex)==TRUE) {
	if (!exists(vertex.cex)) {data(list=vertex.cex)}
	interact.table<-get(vertex.cex)
	names(nodetam)<-rownames(taula)
	sel<-intersect(names(nodetam),names(interact.table))
	nodetam[sel]<-interact.table[sel]/max(interact.table[sel])*1.5+1
	mir<-topTable(obj,"miRNA",pval.cutoff=pval.cutoff,dat.sum=dat.sum)
	sel<-intersect(names(nodetam),names(mir))
	nodetam[sel]<-mir[sel]/max(mir[sel])*1.25+1
}

#print(summary(as.vector(taula.edge)))

if (names==TRUE) {
#plot(taulanet,displaylabels=TRUE,vertex.col=nodeinfo2,
#boxed.labels=FALSE,label.cex=0.6,edge.col=taula.adj,
#vertex.sides=nodeinfo,vertex.rot=noderot,arrowhead.cex = 1,edge.lwd = taula.edge*10)

#plot(taulanet,displaylabels=TRUE,vertex.col=nodeinfo2,
#boxed.labels=FALSE,label.cex=0.5,edge.col=taula.edge,
#vertex.sides=nodeinfo,vertex.rot=noderot,arrowhead.cex = 1,edge.lwd = taula.adj^3*5)

#plot(taulanet,displaylabels=TRUE,vertex.col=nodeinfo2,
#boxed.labels=FALSE,label.cex=0.5,edge.col=taula.adj,
#vertex.sides=nodeinfo,vertex.rot=noderot,arrowhead.cex = 1,edge.lwd = taula.edge,vertex.cex=nodetam)
#print("Hello1")

if (!is.null(dat.sum)==TRUE) {
	plot(taulanet,displaylabels=TRUE,vertex.col=nodeinfo2,
boxed.labels=FALSE,label.cex=0.5,edge.col=taula.adj,
vertex.sides=nodeinfo,vertex.rot=noderot,arrowhead.cex = 1,edge.lwd = taula.edge,vertex.cex=nodetam)

} else {
plot(taulanet,displaylabels=TRUE,vertex.col=nodeinfo2,
boxed.labels=FALSE,label.cex=0.5,edge.col=taula.adj,
vertex.sides=nodeinfo,vertex.rot=noderot,arrowhead.cex = 1,edge.lwd = (taula.edge)^6*10,vertex.cex=nodetam)

}

} else {
#plot(taulanet,displaylabels=FALSE,vertex.col=nodeinfo2,
#boxed.labels=FALSE,label.cex=0.5,edge.col=taula.adj,
#vertex.sides=nodeinfo,vertex.rot=noderot)

#plot(taulanet,displaylabels=FALSE,vertex.col=nodeinfo2,
#boxed.labels=FALSE,label.cex=0.5,edge.col=taula.edge,
#vertex.sides=nodeinfo,vertex.rot=noderot,arrowhead.cex = 1,edge.lwd = taula.adj^3*7.5)

#plot(taulanet,displaylabels=FALSE,vertex.col=nodeinfo2,
#boxed.labels=FALSE,label.cex=0.5,edge.col=taula.adj,
#vertex.sides=nodeinfo,vertex.rot=noderot,arrowhead.cex = 1,edge.lwd = sqrt(taula.edge)*1.5)


#plot(taulanet,displaylabels=FALSE,vertex.col=nodeinfo2,
#boxed.labels=FALSE,label.cex=0.5,vertex.cex=0.5,edge.col=taula.adj,
#vertex.sides=nodeinfo,vertex.rot=noderot,arrowhead.cex = 0.5,edge.lwd = taula.edge*1.5)

if (!is.null(dat.sum)==TRUE) {
plot(taulanet,displaylabels=FALSE,vertex.col=nodeinfo2,
boxed.labels=FALSE,label.cex=1,edge.col=taula.adj,
vertex.sides=nodeinfo,vertex.rot=noderot,arrowhead.cex = 1,edge.lwd = taula.edge,vertex.cex=nodetam)

} else {
plot(taulanet,displaylabels=FALSE,vertex.col=nodeinfo2,
boxed.labels=FALSE,label.cex=1,edge.col=taula.adj,
vertex.sides=nodeinfo,vertex.rot=noderot,arrowhead.cex = 1,edge.lwd = taula.edge^3,vertex.cex=nodetam)}

}

}




plotCircos <- function (obj, pval.cutoff=0.05, dat.sum=obj@info[["dat.sum"]], n=NULL, sub.miRNA=NULL, sub.mRNA=NULL) {

	nf<-layout(mat=matrix(c(1),ncol=1,nrow=1))
	par(mar=c(5.1, 4.1, 4.1, 2.1))
#	plot.new()


	llistacirc<-obj@net

	if (!is.null(sub.miRNA)) {
		llistacirc<-llistacirc[as.character(llistacirc$miRNA) %in% sub.miRNA,]
	}
	if (!is.null(sub.mRNA)) {
		llistacirc<-llistacirc[as.character(llistacirc$mRNA) %in% sub.mRNA,]
	}

	sel<-which(llistacirc$adj.pval<=pval.cutoff)
	if (length(sel)==0) {
		cat(" No miRNA-mRNA pairs meet that pval.cutoff criteria\n")
		stop()
	}

	if (!is.null(dat.sum)) {
		seldat<- which(llistacirc[,"dat.sum"] >= dat.sum)
		sel <- intersect(sel,seldat)
	}
	data<-llistacirc[sel,]


	if (is.null(data$logratio.miRNA)) {data$logratio.miRNA<-rep(0,nrow(data))}
	if (is.null(data$logratio.mRNA)) {data$logratio.mRNA<-rep(0,nrow(data))}
	if (is.null(data$meanExp.miRNA)) {data$meanExp.miRNA<-rep(1,nrow(data))}
	if (is.null(data$meanExp.mRNA)) {data$meanExp.mRNA<-rep(1,nrow(data))}


	#library(RCircos)
	#library(circlize)

	if (!is.null(n)) {
		if (n<nrow(data)) {
			data<-data[with(data,order(adj.pval)),]
			data<-data[1:n,]
		}
	}
	
	if (!exists("genes_human_h37")) { data(genes_human_h37) }
	if (!exists("mirnas_human_17_h37")) { data(mirnas_human_17_h37) }
	#if (!exists("UCSC.HG19.Human.CytoBandIdeogram")) { data(UCSC.HG19.Human.CytoBandIdeogram) }

#	cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
#	chr.exclude <- NULL
#	tracks.inside <- 3
#	tracks.outside <- 0


#layout(matrix(data=c(1,2,3), nrow=1, ncol=3), widths=c(8,1,1),
#heights=8);

#	RCircos.Set.Core.Components(cyto.info,chr.exclude,tracks.inside,tracks.outside)



#	RCircos.Set.Plot.Area()
#	RCircos.Chromosome.Ideogram.Plot()

	mir.hairpins<-as.character(data$miRNA)
	mir.hairpins<-gsub("\\*","",mir.hairpins)
	mir.hairpins<-gsub("miR","mir",mir.hairpins)
	mir.hairpins<-gsub("-5p","",mir.hairpins)
	mir.hairpins<-gsub("-3p","",mir.hairpins)

	links<-data.frame(mirnas_human_17_h37[mir.hairpins,1:3] , genes_human_h37[as.character(data$mRNA),1:3]  )

	miR.exp<-data.frame(miRNA=data$miRNA,logratio=data$logratio.miRNA,meian=data$meanExp.miRNA)

	miR.exp$miRNA<-as.character(miR.exp$miRNA)
	miR.exp$miRNA<-gsub("\\*","",miR.exp$miRNA)
	miR.exp$miRNA<-gsub("miR","mir",miR.exp$miRNA)
	miR.exp$miRNA<-gsub("-5p","",miR.exp$miRNA)
	miR.exp$miRNA<-gsub("-3p","",miR.exp$miRNA)


	miR.exp<-miR.exp[with(miR.exp,order(meian,decreasing=TRUE)),]
	miR.exp<-miR.exp[!duplicated(miR.exp$miRNA),]

	track.miRNAs<-data.frame(mirnas_human_17_h37[miR.exp$miRNA,1:3],exp=miR.exp$meian,logratio=miR.exp$logratio)

#	RCircos.Heatmap.Plot(track.miRNAs,track.num=1,side="in",data.col=5)


	mR.exp<-data.frame(mRNA=data$mRNA,logratio=data$logratio.mRNA,meian=data$meanExp.mRNA)
	mR.exp$mRNA<-as.character(mR.exp$mRNA)

	track.mRNAs<-data.frame(genes_human_h37[mR.exp$mRNA,1:3],exp=mR.exp$meian,logratio=mR.exp$logratio)
#	RCircos.Heatmap.Plot(track.mRNAs,track.num=2,side="in",data.col=5)
#	RCircos.Link.Plot(links,3)
	


#ColorRamp <- RCircos.Get.Heatmap.ColorScales("BlueWhiteRed");
#ColorLevels <- seq(min(miR.exp$logratio),
#max(miR.exp$logratio),
#length=length(ColorRamp));
#par(mai=c(1.5, 0.1, 1.5, 0.5));
#image(1, ColorLevels, matrix(data=ColorLevels,
#ncol=length(ColorLevels), nrow=1),
#col=ColorRamp,
#xlab="", ylab="", xaxt="n");





#ColorRamp <- RCircos.Get.Heatmap.ColorScales("BlueWhiteRed");
#ColorLevels <- seq(min(mR.exp$logratio),
#max(mR.exp$logratio),
#length=length(ColorRamp));
#par(mai=c(1.5, 0.1, 1.5, 0.5));
#image(1, ColorLevels, matrix(data=ColorLevels,
#ncol=length(ColorLevels), nrow=1),
#col=ColorRamp,
#xlab="", ylab="", xaxt="n");



#par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.7)
circos.par("track.height" = 0.05)


circos.initializeWithIdeogram(plotType = c("ideogram", "labels"))

links.filt<-links[complete.cases(links),]


#circos.initialize(factors = a$factor, x = a$x)


a<-track.mRNAs[complete.cases(track.mRNAs),-4]
a<-a[!duplicated(a),]

a$col<-NA
a$col[which(a$logratio>=0)]<-"red"
a$col[which(a$logratio<0)]<-"green"

#a<-a[with(a,order(chr,start)),]

circos.genomicTrackPlotRegion(a[,-5], stack = TRUE, panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = "red", border = "orange")
}, bg.border = NA, bg.col = "white", track.height = 0.05)


#a$text<-rownames(a)

#a[,4]<-a$text
#circos.genomicTrackPlotRegion(a[,-5], ylim = c(0, 1),
# panel.fun = function(region, value, ...) {
# circos.genomicText(region, value, y = 0, labels.column = 1,
# facing = "clockwise", adj = c(0, 0.5), posTransform = posTransform.text,
# cex = 0.3, padding = 0.2)
# }, track.height = 0.1, bg.border = NA)
# i_track = get.cell.meta.data("track.index") # previous track
# circos.genomicPosTransformLines(a[-5],posTransform = quote(posTransform.textregion, y = 0, labels = value[[1]],
# cex = 0.3, padding = 0.2, track.index = i_track)), direction = "outside", track.height=0.05)


a<-track.miRNAs[complete.cases(track.miRNAs),-4]
a<-a[!duplicated(a),]
#a[,4]<-rnorm(nrow(a),0,50)
#rownames(a)<-NULL
#colnames(a)<-c("chr","start","start.1","value")

#a$col<-NA
#a$col[which(a$logratio>=0)]<-"red"
#a$col[which(a$logratio<0)]<-"green"


circos.genomicTrackPlotRegion(a[,-5], stack = TRUE, panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = "red", border = "blue")
}, bg.border = NA, bg.col = "white", track.height = 0.05)

#a$text<-rownames(a)

#a[,4]<-a$text
#circos.genomicTrackPlotRegion(a[,-5], ylim = c(0, 1),
# panel.fun = function(region, value, ...) {
# circos.genomicText(region, value, y = 0, labels.column = 1,
# facing = "clockwise", adj = c(0, 0.5), posTransform = posTransform.text,
# cex = 0.3, padding = 0.2)
# }, track.height = 0.05, bg.border = NA,col="violet")

# i_track = get.cell.meta.data("track.index") # previous track
# circos.genomicPosTransformLines(a[-5],posTransform = quote(posTransform.text(region, y = 0, labels = value[[1]],
# cex = 0.3, padding = 0.2, track.index = i_track)), direction = "outside",
# track.height = 0.05,col="violet")



circos.genomicLink(data.frame(links.filt[,1:3],1), data.frame(links.filt[,4:6],1), col = sample(10, nrow(links.filt), replace = TRUE),border="gray")

}










writeExcel <- function (obj, name, pval.cutoff=0.05, cor=NULL, alternative="less", dat.sum=obj@info[["dat.sum"]], slot="net", pval="adj.pval") {
	#libary(WriteXLS)

	if (slot == "net") {
		res <- obj@net[which(obj@net[,pval] <= pval.cutoff),]
		if (!is.null(dat.sum)) {
			res <- res[which(res$dat.sum>=dat.sum),]
		}
		if (!is.null(cor)) {
		  if (alternative=="less") {
		    selcor<-which(obj@net$cor <= cor)
		  }
		  if (alternative=="greater") {
		    selcor<-which(obj@net$cor >= cor)
		  }
		  if (alternative=="two.sided") {
		    selcor<-which(abs(obj@net$cor) >= abs(cor))		    
		  }
		  sel <- intersect(sel,selcor)  
		}
		

		write.datt<-list(
			#miRNAdat<-data.frame(obj@dat.miRNA),
			#mRNAdat<-data.frame(obj@dat.mRNA),
			#common<-data.frame(obj@common),
			results=res
			#diffexp.miRNA<-data.frame(obj@diffexp.miRNA),
			#diffexp.mRNA<-data.frame(obj@diffexp.mRNA)
		)
	}

	if (slot == "diffexp.miRNA") {
		write.datt<-list(diffexp.miRNA=data.frame(obj@diffexp.miRNA[which(obj@diffexp.miRNA[,pval] <= pval.cutoff),]))
	}
	if (slot == "diffexp.mRNA") {
		write.datt<-list(diffexp.mRNA=data.frame(obj@diffexp.mRNA[which(obj@diffexp.mRNA[,pval] <= pval.cutoff),]))
	}
	if (slot == "dat.miRNA") {
		write.datt<-list(dat.mRNA=data.frame(obj@dat.miRNA))
	}
	if (slot == "dat.mRNA") {
		write.datt<-list(dat.mRNA=data.frame(obj@dat.mRNA))
	}
	if (slot == "pheno.miRNA") {
		write.datt<-list(pheno.miRNA=data.frame(obj@pheno.miRNA))
	}
	if (slot == "pheno.mRNA") {
		write.datt<-list(pheno.mRNA=data.frame(obj@pheno.mRNA))
	}

	WriteXLS("write.datt",name)
}


writeCsv <- function (obj, name, pval.cutoff=1, cor=NULL, alternative="less", dat.sum=obj@info[["dat.sum"]], slot="net", pval="adj.pval") {
	#libary(WriteXLS)

	if (slot == "net") {
		sel<-which(obj@net$adj.pval<=pval.cutoff)
		if (!is.null(dat.sum)) {
			seldat<- which(obj@net[,"dat.sum"] >= dat.sum)
			sel <- intersect(sel,seldat)
		}
		if (!is.null(cor)) {
		  if (alternative=="less") {
		    selcor<-which(obj@net$cor <= cor)
		  }
		  if (alternative=="greater") {
		    selcor<-which(obj@net$cor >= cor)
		  }
		  if (alternative=="two.sided") {
		    selcor<-which(abs(obj@net$cor) >= abs(cor))		    
		  }
		  sel <- intersect(sel,selcor)  
		}
		write.datt<-data.frame(obj@net[sel,])
	}
	if (slot == "diffexp.miRNA") {
		sel<-which(obj@diffexp.miRNA$adj.pval<=pval.cutoff)
		write.datt<-data.frame(obj@diffexp.miRNA[sel,])
	}
	if (slot == "diffexp.mRNA") {
		sel<-which(obj@diffexp.mRNA$adj.pval<=pval.cutoff)
		write.datt<-data.frame(obj@diffexp.mRNA[sel,])
	}

	if (slot == "dat.miRNA") {
		write.datt<-data.frame(obj@dat.miRNA)
	}
	if (slot == "dat.mRNA") {
		write.datt<-data.frame(obj@dat.mRNA)
	}
	if (slot == "pheno.miRNA") {
		write.datt<-data.frame(obj@pheno.miRNA)
	}
	if (slot == "pheno.mRNA") {
		write.datt<-data.frame(obj@pheno.mRNA)
	}

	write.table(write.datt,name,sep="\t",quote=FALSE,row.names=FALSE)
}


writeSif<- function( obj, file, pval.cutoff=0.05, cor=NULL, alternative="less", dat.sum=obj@info[["dat.sum"]], add.other=NULL, sub.miRNA=NULL, sub.mRNA=NULL, expand=FALSE, vertex.cex="interact.table") {
	sel<-which(obj@net$adj.pval<=pval.cutoff)
	if (!is.null(dat.sum)) {
		seldat<- which(obj@net[,"dat.sum"] >= dat.sum)
		sel <- intersect(sel,seldat)
	}
	if (!is.null(cor)) {
	  if (alternative=="less") {
	    selcor<-which(obj@net$cor <= cor)
	  }
	  if (alternative=="greater") {
	    selcor<-which(obj@net$cor >= cor)
	  }
	  if (alternative=="two.sided") {
	    selcor<-which(abs(obj@net$cor) >= abs(cor))		    
	  }
	  sel <- intersect(sel,selcor)  
	}
			if (!is.null(sub.miRNA)) {
			selmirna<-which(obj@net$miRNA %in% sub.miRNA)
			sel <- intersect(sel,selmirna)
		}
		if (!is.null(sub.mRNA)) {
			selmrna<-which(obj@net$mRNA %in% sub.mRNA)
			sel <- intersect(sel,selmirna)
		}
	results<-data.frame(miRNA=obj@net$miRNA[sel],type=rep("pd",length(sel)),mRNA=obj@net$mRNA[sel])
	if (!is.null(add.other)) {
	#data(interact)
		interact<-get(add.other)
if (!expand) {
	add<-interact[which((as.character(interact[,1]) %in% as.character(results$mRNA) & as.character(interact[,3]) %in% as.character(results$mRNA) )==TRUE) , ] 
} else {
	add<-interact[which((as.character(interact[,1]) %in% as.character(results$mRNA) | as.character(interact[,3]) %in% as.character(results$mRNA) )==TRUE & 
(as.character(interact[,1]) %in% as.character(obj@sig.mRNA) & as.character(interact[,3]) %in% as.character(obj@sig.mRNA) )==TRUE) , ]
}
		add<-rbind(add,add[,c(3,2,1)])
		colnames(add)<-c("miRNA","type","mRNA")
		results<-rbind(results,add)
	}
	write.table(results, paste(file,".sif",sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

	#write node attributes
	node.attr<-data.frame(names=unique(c(as.character(results[,1]),as.character(results[,3]))))
	rownames(node.attr)<-node.attr$names
	sel<-intersect(rownames(obj@diffexp.miRNA),node.attr$names)
	node.attr$logratio<-rep(NA,nrow(node.attr))
	node.attr$type<-rep(NA,nrow(node.attr))
	node.attr$size<-rep(NA,nrow(node.attr))
	node.attr[sel,"logratio"]<-obj@diffexp.miRNA[sel,"logratio"]
	node.attr[sel,"type"]<-"miRNA"
	tam.miRNA<-table(as.character(results[,1]))
	node.attr[sel,"size"]<-tam.miRNA[sel]

	sel<-intersect(rownames(obj@diffexp.mRNA),node.attr$names)
	node.attr[sel,"logratio"]<-obj@diffexp.mRNA[sel,"logratio"]
	node.attr[sel,"type"]<-"mRNA"
#	data(list=vertex.cex)
	if (!exists(vertex.cex)) {data(list=vertex.cex)}
	node.attr[sel,"size"]<-get(vertex.cex)[sel]

	write.table(node.attr,paste(file,"_node_attributes.csv",sep="") , sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
#	falta posar el tamany els nodes




	#write edge attributes
	edge.attr<-data.frame(node1=as.character(results[,1]),int=rep("pd",nrow(results)),node2=as.character(results[,3]))
	rownames(edge.attr)<-paste(edge.attr[,1],edge.attr[,3],sep=":")
	sel<-intersect(rownames(edge.attr),rownames(obj@net))
	edge.attr$score<-rep(NA,nrow(edge.attr))
	edge.attr[sel,"score"]<-obj@net[sel,"score"]
	edge.attr$dat.sum<-rep(NA,nrow(edge.attr))
	edge.attr[sel,"dat.sum"]<-obj@net[sel,"dat.sum"]

	write.table(edge.attr,paste(file,"_with_attributes.sif",sep="") , sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


}


openCytoscape <- function (obj=NULL, pval.cutoff=0.05, dat.sum=obj@info[["dat.sum"]], file=NULL, cytoscape.folder="/home/mvila/Cytoscape_v2.8.3", sub.miRNA=NULL, sub.mRNA=NULL, add.other=NULL,expand=FALSE) {
	if (is.null(file)) {
		file<-"network_default.sif"
#		if (!is.null(sub.miRNA)) {
#			obj@net<-obj@net[obj@net$miRNA %in% sub.miRNA,]
#		}
#		if (!is.null(sub.mRNA)) {
#			obj@net<-obj@net[obj@net$mRNA %in% sub.mRNA,]
#		}
		writeSif(obj, file, pval.cutoff, dat.sum=dat.sum, add.other=add.other, sub.miRNA=sub.miRNA, sub.mRNA=sub.mRNA, expand=expand)
	}
	jar.location<-paste(cytoscape.folder,"/cytoscape.jar",sep="")
	plugin.location<-paste(cytoscape.folder,"/plugins",sep="")
	system(paste("java -Xmx800M -jar",jar.location,"-N",file,"-p",plugin.location,"&"))
}



GOanalysis <- function (obj, type, ontology, pval.cutoff = 0.05, dat.sum=obj@info[["dat.sum"]], score.cutoff = NULL, sub.miRNA=NULL, exclude.miRNA=NULL, sub.mRNA=NULL, organism="human", FC=NULL, up=FALSE, dw=FALSE, add.miRNA=FALSE) 
{
    	#funció per anotar
  #  	annotentrez <- 
	
	if (organism=="human") {annotation <- "org.Hs.eg.db"}
	if (organism=="mouse") {annotation <- "org.Mm.eg.db"}	

	#library("org.Mm.eg.db")
	#seleccionar
	sel <- 1:nrow(obj@net)

	if (!is.null(obj@net$adj.pval)) {
	    	selpval <- which(obj@net$adj.pval <= pval.cutoff)
		sel <- intersect(sel,selpval)
	}

	if (!is.null(score.cutoff)) {
		selscore<- which(obj@net[,"score"] >= score.cutoff)
		sel <- intersect(sel,selscore)
	}

	if (!is.null(FC)) {
		selFC<- which(abs(obj@net[,"logratio.mRNA"]) >= foldchange2logratio(FC))
		sel <- intersect(sel,selFC)
	}

	if (up) {
		selup<- which(obj@net[,"logratio.mRNA"] > 0)
		sel <- intersect(sel,selup)
	}
	if (dw) {
		seldw<- which(obj@net[,"logratio.mRNA"] < 0)
		sel <- intersect(sel,seldw)
	}
	if (!is.null(dat.sum)) {
		seldat<- which(obj@net[,"dat.sum"] >= dat.sum)
		sel <- intersect(sel,seldat)
	}
	
	subs1 <- obj@net[sel,]

	if (!is.null(sub.miRNA)) {
		subs1 <- subs1[subs1$miRNA %in% sub.miRNA,]

	}
	
    	mRNA.name <- as.character(unique(subs1[,"mRNA"]))
	
	if (!is.null(sub.mRNA)) {
		mRNA.name<-sub.mRNA
	}

	if (!is.null(exclude.miRNA)) {
		subs1 <- obj@net[which(obj@net$adj.pval <= pval.cutoff),]
		subs1 <- subs1[(subs1$miRNA %in% exclude.miRNA),]
		subs <- as.character(unique(subs1[,"mRNA"]))
		mRNA.name <- setdiff(mRNA.name, subs)
	}

	if (organism=="human") {
    		data(conversor)
    		mRNA.id <- as.character(unique(conversor[mRNA.name, "Entrez.Gene.ID"]))
    		#all <- unique(conversor[, "Entrez.Gene.ID"])
    		entrezIds.all <- as.character(unique(conversor[, 2]))
	}

	if (organism=="mouse") {
    		data(conversor.mouse)
    		mRNA.id <- as.character(unique(conversor.mouse[mRNA.name, "Entrez.Gene.ID"]))
    		#all <- unique(conversor[, "Entrez.Gene.ID"])
    		entrezIds.all <- as.character(unique(conversor.mouse[, 2]))
	}


	#print(length(mRNA.id))
	#print("\n")
	#print(length(entrezIds.all))

	if (type=="GO") {
    	params <- new("GOHyperGParams", 
		geneIds = mRNA.id, 
		universeGeneIds = entrezIds.all, 
    		annotation = annotation, 
		ontology = ontology, 
		pvalueCutoff = 1, 
    		conditional = FALSE, 
		testDirection = "over")
	}

	if (type=="KEGG") {
    	params <- new("KEGGHyperGParams", 
		geneIds = mRNA.id, 
		universeGeneIds = entrezIds.all, 
    		annotation = annotation, 
		pvalueCutoff = 1, 
		testDirection = "over")
	}

	#assignInNamespace("lapply", lapply, "base")
	if (type=="GO" | type=="KEGG") {
    	results.test <- hyperGTest(params)
    	df1 <- summary(results.test)
	xxx <- geneIdsByCategory(results.test)
    	genescat <- unlist(lapply(xxx, function(mapped_genes) {
	    if (organism=="human") {
    	    	x <- org.Hs.egSYMBOL
	    }
	    if (organism=="mouse") {
    	    	x <- org.Mm.egSYMBOL
	    }

    	    mapped_genes <- as.character(mapped_genes)
	#cat(class(x))
    	    xx <- as.list(x[mapped_genes])
    	    return(paste(unlist(xx), collapse = ", "))
    	}))[df1[,1]]


    	fdr <- p.adjust(df1$Pvalue, method = "BH")
    	GO.results <- data.frame(Ontology = rep(ontology, length = length(fdr)), 
        	df1, fdr, genescat)

	if (add.miRNA) {
		GO.results$miRNAs.count<-NA
		GO.results$miRNAs<-NA
		GO.results$genescat<-as.character(GO.results$genescat)
		for (i in 1:nrow(GO.results)) {
			genes <- strsplit(GO.results$genescat[i],", ")[[1]]
			sels <- which((obj@net$mRNA %in% genes) & (obj@net$adj.pval <= pval.cutoff))
			miRNAs <- unique(obj@net[sels,"miRNA"])
			GO.results$miRNAs.count[i]<-length(miRNAs)
			GO.results$miRNAs[i]<-paste(miRNAs,collapse=", ")
		}
	}

	}

	if (type=="REACTOME") {

		mRNA.id<-mRNA.id[which(!is.na(mRNA.id)==TRUE)]
		df1<-summary(enrichPathway(gene=mRNA.id, organism=organism, pvalueCutoff=1, readable=TRUE))

		GO.results<-data.frame(df1)

	}

	obj@GO.results[[paste(type,ontology,sep=":")]]<-GO.results
	obj@info[[paste(type,ontology,sep=":")]]<-list(pval.cutoff=pval.cutoff, dat.sum=dat.sum, sub.miRNA=sub.miRNA, exclude.miRNA=exclude.miRNA, sub.mRNA=sub.mRNA, organism=organism, FC=FC, up=up, dw=dw)
    	return(obj)
}


topTable <- function (obj, class, pval.cutoff=0.05, dat.sum=obj@info[["dat.sum"]], score.cutoff=NULL, plot=FALSE, names=FALSE, n=NULL, remove.names=FALSE, table.items = NULL) {

	sel<-1:nrow(obj@net)
	if (!is.null(obj@net$adj.pval)) {
		seldat<-which(obj@net$adj.pval<=pval.cutoff)
		sel <- intersect(sel,seldat)
	}

	if (!is.null(dat.sum)) {
		seldat<- which(obj@net[,"dat.sum"] >= dat.sum)
		sel <- intersect(sel,seldat)
	}

	if (!is.null(score.cutoff)) {
		seldat<- which(obj@net[,"score"] >= score.cutoff)
		sel <- intersect(sel,seldat)
	}


	sub<-obj@net[sel,]
	if (class=="miRNA") {
		res<-table(sub$miRNA)
		res<-sort(res,decreasing=TRUE)
		mirs<-names(res)
		resnum<-as.numeric(res)
		if (names) {
			total.regulated<-vector()
			perc.regulated<-vector()
			freq<-vector()
			for (i in 1:length(mirs)) {
				if (is.null(table.items)) {
				  freq[i]<-paste(sub$mRNA[sub$miRNA==mirs[i]],collapse=", ")
				}
			  else {
				  freq[i]<-paste(sub$mRNA[sub$miRNA==mirs[i]][1:min(table.items,length(which(sub$miRNA==mirs[i])))],collapse=", ")
				}
				total.regulated<-unique(c(total.regulated,as.character(sub$mRNA[sub$miRNA==mirs[i]])))
				perc.regulated[i]<-length(total.regulated)
			}
			res<-data.frame(freq=resnum,names=freq,perc.regulated=perc.regulated/length(obj@sig.mRNA)*100)
			rownames(res)<-mirs
		} else {
			#res<-resnum
			#names(res)<-mirs
		}
	}

	if (class=="mRNA") {
		res<-table(sub$mRNA)
		res<-sort(res,decreasing=TRUE)
		mrnas<-names(res)
		resnum<-as.numeric(res)
		if (names) {
			total.regulated<-vector()
			perc.regulated<-vector()
			freq<-vector()
			for (i in 1:length(mrnas)) {
			  if (is.null(table.items)) {
			    freq[i]<-paste(sub$miRNA[sub$mRNA==mrnas[i]],collapse=", ")
			  }
			  else {
			    freq[i]<-paste(sub$miRNA[sub$mRNA==mrnas[i]][1:min(table.items,length(which(sub$mRNA==mrnas[i])))],collapse=", ")
			  }
			  #				freq[i]<-paste(sub$miRNA[sub$mRNA==mrnas[i]],collapse=", ")
				total.regulated<-unique(c(total.regulated,sub$miRNA[sub$mRNA==mrnas[i]]))
				perc.regulated[i]<-length(total.regulated)

			}
			res<-data.frame(freq=resnum,names=freq,perc.regulated=perc.regulated/length(obj@sig.mRNA)*100)
			rownames(res)<-mrnas
			#rownames(res)<-names(res)
		} else {
		#	res<-data.frame(freq=resnum)
		#	names(res)<-mrnas
		}

	}

#	if (names) {
#		res<-res[with(res,order(freq,decreasing=TRUE)),]
#	} else {
#		res<-sort(res,decreasing=TRUE)
#	}
	

	if(plot) {


		nf<-layout(mat=matrix(c(1),ncol=1,nrow=1))
		par(mar=c(5.1, 4.1, 4.1, 2.1))

		if (names) {
			if (is.null(n)) {n<-nrow(res)}
			barplot(res$freq[1:n])
			#lines(res$perc.regulated/100*res[1,1],col="red")
			par(new = TRUE)
			plot(1:n,res$perc.regulated[1:n], col = "red", axes = FALSE, xlab = "miRNAs", ylab = "Number of targets",type="l")
			axis(side = 4)
			mtext("Percentage of regulated",side=4,line=-1)
		} else {
			if (is.null(n)) {n<-length(res)}
			if (remove.names) {names(res)<-NULL;xlab<-class}
			else {xlab<-NA}
			barplot(res[1:n],xlab = xlab)
		}
	} else {
		return(res)
	}

}









plotGO <- function ( obj, type, ontology, fdr=0.05, filename="GO_tree_default" ) {
	#library("RamiGO")
	GOres<-obj@GO.results[[paste(type,ontology,sep=":")]]
	sel<-which(GOres$fdr<fdr)
	goIDs<-rownames(GOres)[sel]
	num.fdr<-log(GOres$fdr[sel])
	#remapejar de vermell a groc
	
	remap <- function(x) { (( x ) / max( abs(x) ) / 2)+0.5 }  # map x onto [0, 1]
	fun.col <- function(x) {rgb(colorRamp(c("red","gold", "yellow"))(remap(x)),
                        maxColorValue = 255)
	}

	col.fdr<-fun.col(as.numeric(num.fdr))
	pngRes<-getAmigoTree(goIDs=goIDs, color=col.fdr, filename=filename,picType="png",saveResult=TRUE)

}


plotDensity <- function (obj, subset, col.color=1, colors=c("turquoise", "violet")) {

	nf<-layout(mat=matrix(c(1),ncol=1,nrow=1))
	par(mar=c(5.1, 4.1, 4.1, 2.1))
#	plot.new()


	if (subset=="miRNA") {
		data<-obj@dat.miRNA

		if (length(levels(as.factor(obj@pheno.miRNA[,col.color])))==2) {
			color<-as.character(factor(as.numeric(as.factor(obj@pheno.miRNA[,col.color])),labels=colors))
		} else {
			color<-as.character(factor(as.numeric(as.factor(obj@pheno.miRNA[,col.color]))+1))
		}

	} else { if (subset=="mRNA") {
			data<-obj@dat.mRNA

			if (length(levels(as.factor(obj@pheno.mRNA[,col.color])))==2) {
				color<-as.character(factor(as.numeric(as.factor(obj@pheno.mRNA[,col.color])),labels=colors))
			} else {
				color<-as.character(factor(as.numeric(as.factor(obj@pheno.mRNA[,col.color]))+1))
			}

		}
	}

	#column<-as.factor(obj@pheno.miRNA[,group])
	#color<-as.numeric(column)


	ymax<-max(density(data[,1])$y)*1.25
	plot(density(data[,1])$x,density(data[,1])$y,ylim=c(0,ymax),main="Per-sample density", col=color[1],type="l",xlab="Expression",ylab="Density")
	for (i in 2:length(color)) {
		lines(density(data[,i]),col=color[i])

	}
	legend("topright",levels(as.factor(obj@pheno.miRNA[,col.color])),col=levels(as.factor(color)),lwd=1)
}



summary.corObject <- function (object, ...) {

	if (validObject(object)) {

	cat("corObject with:\n")
	cat(paste(" miRNA slot with",ncol(object@dat.miRNA),"samples and",nrow(object@dat.miRNA),"probesets\n"))
	cat(paste(" mRNA slot with",ncol(object@dat.mRNA),"samples and",nrow(object@dat.mRNA),"probesets\n"))
	cat("Computations done:\n")
	if (dim(object@diffexp.mRNA)[1]>1 ) {
		cat(paste("- Differential expression mRNA: ",object@info[["mRNA.diffexp.method"]][1],"method used\n"))
		cat(paste("                                ",object@info[["mRNA.diffexp.method"]][2],"comparison used\n"))
		
	}

	if (dim(object@diffexp.miRNA)[1]>1 ) {
		cat(paste("- Differential expression miRNA: ",object@info[["miRNA.diffexp.method"]][1],"method used\n"))
		cat(paste("                                 ",object@info[["miRNA.diffexp.method"]][2],"comparison used\n"))

	}

	if (dim(object@cor)[1]>1 | dim(object@cor)[2]>1) {
		cat(paste("- Correlation:  \"",object@info[["correlation.type"]],"\" method used\n",sep=""))
		cat(paste("                \"",object@info[["correlation.function.used"]],"\" function used\n",sep=""))
		cat(c("               ",length(object@info[["correlation.samples.used"]]),"samples used\n"))
		cat(c("               ",dim(object@cor)[1],"miRNAs used\n"))
		cat(c("               ","  ",paste(object@info[["miRNA.criteria"]],collapse="; "),"\n"))

		cat(c("               ",dim(object@cor)[2],"mRNAs used\n"))
		cat(c("               ","  ",paste(object@info[["mRNA.criteria"]],collapse="; "),"\n"))

	}

	if (!is.null(object@info[["database"]])) {
		cat(paste("- Database:  \"",object@info[["database"]],"\" database used\n",sep=""))
	}

	if (!is.null(object@net$p.comb)) {
		cat(paste("- P.value combination:  \"",object@info[["pcomb.method"]],"\" method used\n",sep=""))
	}

	if (!is.null(object@net$adj.pval)) {
		cat(paste("- P.value adjustment:  \"",object@info[["padjust.method"]],"\" method used\n",sep=""))
	}

	}

}



##### Make a report

mkReport <- function (obj, file, title="Default \\texttt{miRComb} output", dat.sum.table=NULL) {
	
	cat("Doesn't work locally? Try loading your \".RData\" file in our server, and we will make the pdf report for you (up to 100Mb):\n")
	cat("http://bioinfo.ciberehd.org/mircomb/mkreport.html")

	#set random seed
	seed<-paste("_rndm_seed_",sample(1:100000)[1],sep="")

	if (is.null(obj@info[["dat.sum"]])) {
		obj@info[["dat.sum"]]<-0
	}

	n.samp<-25

	sink(paste(file,".tex",sep=""))

	cat("
\\documentclass[a4paper,11pt]{article}
\\usepackage[latin1]{inputenc}
\\usepackage[english]{babel}
%\\selectlanguage{catalanb}
\\usepackage[T1]{fontenc}
%\\usepackage{lucidabr}
%\\usepackage{appendix}
\\usepackage{color}
\\usepackage[table]{xcolor}

\\usepackage{enumerate}
\\usepackage{amsmath}
\\usepackage{amssymb}
%\\usepackage[landscape]{geometry}
%\\usepackage{pdflscape}
\\usepackage{comment}
%\\usepackage{colortbl}
\\usepackage{pgfplots}
%\\pgfplotsset{compat=1.6}
\\usepackage{tikz}
%\\usepackage{fontspec}
%\\usepackage{subfigure}
\\usepackage{caption}
\\usepackage{subcaption}
\\usepackage{hyperref}
\\setlength{\\captionmargin}{30pt}

\\usepackage[paperheight=297mm,paperwidth=210mm,top=25mm,left=25mm,height=255mm,width=165mm]{geometry}

%\\usepackage{bbm}
%\\setmonofont{Lucida Console}

%\\usepackage{Sweave}

")

cat(paste("\\title{",title,"}",sep=""))
cat("
%\\author{CIBEREHD Bioinformatics Platform}
")
name<-gsub("_","\\\\_",getwd())
cat(paste("\\author{",name,"}",sep=""))
cat("
\\begin{document}
\\maketitle
")

cat("\\section{Exploratory analysis of miRNA dataset}")

cat("
\\renewcommand{\\arraystretch}{1.1}
\\definecolor{gray09}{gray}{0.9}
\\definecolor{gray099}{gray}{0.99}
\\rowcolors{1}{gray09}{gray099}
\\begin{table}[!h]
\\centering
\\begin{tabular}{lp{8cm}}
")

cat(paste("Number of miRNAs analysed &", nrow(obj@diffexp.miRNA) ,"\\\\")
)
cat(paste("Number of samples &", ncol(obj@dat.miRNA) ,"\\\\"))
if (length(colnames(obj@dat.miRNA))<n.samp) {
cat(paste("Samples &", gsub("\\_","\\\\_",paste(colnames(obj@dat.miRNA),collapse=", ")) ,"\\\\"))
}
cat("
\\end{tabular}
\\caption{Basic information of the miRNA dataset.}
\\end{table}
")


if (nrow(obj@pheno.miRNA)<=n.samp) {
	print(xtable(obj@pheno.miRNA,caption="Phenotypical information of the miRNA dataset.",latex.environments="center"))
} else {
	print(xtable(summary(obj@pheno.miRNA),caption="Summary of the phenotypical information of the miRNA dataset.",latex.environments="center"))
}



pdf(paste("pcamiRNA",seed,".pdf",sep=""))
plotPca(obj,"miRNA")
dev.off()

pdf(paste("densmiRNA",seed,".pdf",sep=""))
plotDensity(obj,"miRNA")
dev.off()


cat("
\\begin{figure}[!h]
\\centering
")
cat(paste("\\includegraphics[width=0.46\\textwidth]{",paste("pcamiRNA",seed,".pdf",sep=""),"}",sep=""))
cat("\\hspace{0.02\\textwidth}")
cat(paste("\\includegraphics[width=0.46\\textwidth]{",paste("densmiRNA",seed,".pdf",sep=""),"}",sep=""))


cat("
\\caption[]{PCA and density plot for miRNAs.}
\\end{figure}
")




cat("\\newpage")

cat("\\section{Exploratory analysis of mRNA dataset}")

cat("
\\renewcommand{\\arraystretch}{1.1}
\\definecolor{gray09}{gray}{0.9}
\\definecolor{gray099}{gray}{0.99}
\\rowcolors{1}{gray09}{gray099}
\\begin{table}[!h]
\\centering
\\begin{tabular}{lp{8cm}}
")

cat(paste("Number of mRNAs analysed &", nrow(obj@diffexp.mRNA) ,"\\\\")
)
cat(paste("Number of samples &", ncol(obj@dat.mRNA) ,"\\\\"))
if (length(colnames(obj@dat.mRNA))<n.samp) {
cat(paste("Samples &", gsub("\\_","\\\\_",paste(colnames(obj@dat.mRNA),collapse=", ")) ,"\\\\"))
}
cat("
\\end{tabular}
\\caption{Basic information of the mRNA dataset.}
\\end{table}
")


if (nrow(obj@pheno.mRNA)<=n.samp) {
	print(xtable(obj@pheno.mRNA,caption="Phenotypical information of the mRNA dataset.",latex.environments="center"))
} else {
	print(xtable(summary(obj@pheno.mRNA),caption="Summary of the phenotypical information of the mRNA dataset.",latex.environments="center"))
}



pdf(paste("pcamRNA",seed,".pdf",sep=""))
plotPca(obj,"mRNA")
dev.off()

pdf(paste("densmRNA",seed,".pdf",sep=""))
plotDensity(obj,"mRNA")
dev.off()


cat("
\\begin{figure}[!h]
\\centering
")
cat(paste("\\includegraphics[width=0.46\\textwidth]{",paste("pcamRNA",seed,".pdf",sep=""),"}",sep=""))
cat("\\hspace{0.02\\textwidth}")
cat(paste("\\includegraphics[width=0.46\\textwidth]{",paste("densmRNA",seed,".pdf",sep=""),"}",sep=""))


cat("
\\caption[]{PCA and density plot for mRNAs.}
\\end{figure}
")



if (!all(dim(obj@diffexp.miRNA)==0)) {

cat("\\newpage")

cat("\\section{Differentially expressed miRNAs}")

cat("
\\renewcommand{\\arraystretch}{1.1}
\\definecolor{gray09}{gray}{0.9}
\\definecolor{gray099}{gray}{0.99}
\\rowcolors{1}{gray09}{gray099}
\\begin{table}[!h]
\\centering
\\begin{tabular}{lp{8cm}}
")

#cat(paste("Number of miRNAs analysed &", nrow(obj@diffexp.miRNA) ,"\\\\"))
cat(paste("Analysis performed & Comparative used: ",obj@info[["miRNA.diffexp.method"]][2],"; method used: ",obj@info[["miRNA.diffexp.method"]][1],".\\\\",sep=""))
cat(paste("Number of differentially expressed miRNAs &", length(obj@sig.miRNA) , " (",
length(which(obj@diffexp.miRNA[obj@sig.miRNA,"logratio"]>0))," upregulated, ",
length(which(obj@diffexp.miRNA[obj@sig.miRNA,"logratio"]<0))," downregulated)", "\\\\")
,sep="")
cat(paste("Number of samples &", ncol(obj@dat.miRNA) ,"\\\\"))
if (length(colnames(obj@dat.miRNA))<n.samp) {
cat(paste("Samples &", gsub("\\_","\\\\_",paste(colnames(obj@dat.miRNA),collapse=", ")) ,"\\\\"))
}
cat(paste("Criteria for selecting miRNAs &", obj@info$miRNA.criteria ,"\\\\"))
cat("
\\end{tabular}
\\caption{Basic statistics}
\\end{table}
")


if (obj@info[["miRNA.diffexp.method"]][1] != "anova") {
pdf(paste("mamiRNA",seed,".pdf",sep=""))
plot(obj@diffexp.miRNA$logratio,-log10(obj@diffexp.miRNA$pval),xlab="logratio",ylab="-log10(pval)",pch=19,cex=0.75)
points(obj@diffexp.miRNA[obj@sig.miRNA,"logratio"],-log10(obj@diffexp.miRNA[obj@sig.miRNA,"pval"]),col="red",pch=19,cex=0.75)
dev.off()
}

pdf(paste("heatmapmiRNA",seed,".pdf",sep=""))
plotHeatmap(obj,"miRNA")
dev.off()

cat("
\\begin{figure}[!h]
\\centering
")
cat(paste("\\includegraphics[width=0.46\\textwidth]{",paste("heatmapmiRNA",seed,".pdf",sep=""),"}",sep=""))
cat("\\hspace{0.02\\textwidth}")
if (obj@info[["miRNA.diffexp.method"]][1] != "anova") {
cat(paste("\\includegraphics[width=0.46\\textwidth]{",paste("mamiRNA",seed,".pdf",sep=""),"}",sep=""))
}
cat(paste("\\caption[]{A) Heatmap vith the top 50 most significant miRNAs (sorted by adjusted p-value). B) Volcano plot showing the selected miRNAs.}
\\end{figure}
",sep=""))


}

if (!all(dim(obj@diffexp.mRNA)==0)) {


cat("\\newpage")

cat("\\section{Differentially expressed mRNAs}")

cat("
\\renewcommand{\\arraystretch}{1.1}
\\definecolor{gray09}{gray}{0.9}
\\definecolor{gray099}{gray}{0.99}
\\rowcolors{1}{gray09}{gray099}
\\begin{table}[!h]
\\centering
\\begin{tabular}{lp{8cm}}
")

#cat(paste("Number of mRNAs analysed &", nrow(obj@diffexp.mRNA) ,"\\\\"))
cat(paste("Analysis performed & Comparative used: ",obj@info[["mRNA.diffexp.method"]][2],"; method used: ",obj@info[["mRNA.diffexp.method"]][1],".\\\\",sep=""))
cat(paste("Number of differentially expressed mRNAs &", length(obj@sig.mRNA) , " (",
length(which(obj@diffexp.mRNA[obj@sig.mRNA,"logratio"]>0))," upregulated, ",
length(which(obj@diffexp.mRNA[obj@sig.mRNA,"logratio"]<0))," downregulated)", "\\\\")
,sep="")
cat(paste("Number of samples &", ncol(obj@dat.mRNA) ,"\\\\"))
if (length(colnames(obj@dat.mRNA))<n.samp) {
cat(paste("Samples &", gsub("\\_","\\\\_",paste(colnames(obj@dat.mRNA),collapse=", ")) ,"\\\\"))
}
cat(paste("Criteria for selecting mRNAs &", obj@info$mRNA.criteria ,"\\\\"))
cat("
\\end{tabular}
\\caption{Basic statistics}
\\end{table}
")

if (obj@info[["mRNA.diffexp.method"]][1] != "anova") {
pdf(paste("mamRNA",seed,".pdf",sep=""))
plot(obj@diffexp.mRNA$logratio,-log10(obj@diffexp.mRNA$pval),xlab="logratio",ylab="-log10(pval)",pch=19,cex=0.75)
points(obj@diffexp.mRNA[obj@sig.mRNA,"logratio"],-log10(obj@diffexp.mRNA[obj@sig.mRNA,"pval"]),col="red",pch=19,cex=0.75)
dev.off()
}

pdf(paste("heatmapmRNA",seed,".pdf",sep=""))
plotHeatmap(obj,"mRNA")
dev.off()


pdf(paste("pcamRNA",seed,".pdf",sep=""))
plotPca(obj,"mRNA")
dev.off()

cat("
\\begin{figure}[!h]
\\centering
")
cat(paste("\\includegraphics[width=0.46\\textwidth]{",paste("heatmapmRNA",seed,".pdf",sep=""),"}",sep=""))
cat("\\hspace{0.02\\textwidth}")
if (obj@info[["mRNA.diffexp.method"]][1] != "anova") {
cat(paste("\\includegraphics[width=0.46\\textwidth]{",paste("mamRNA",seed,".pdf",sep=""),"}",sep=""))
}
cat(paste("\\caption[]{A) Heatmap vith the top 50 most significant mRNAs (sorted by adjusted p-value). B) Volcano plot showing the selected mRNAs.}
\\end{figure}
",sep=""))


}

cat("\\clearpage")

cat("\\section{Correlation \\& intersection with databases}")


cat("
\\renewcommand{\\arraystretch}{1.1}
\\definecolor{gray09}{gray}{0.9}
\\definecolor{gray099}{gray}{0.99}
\\rowcolors{1}{gray09}{gray099}
\\begin{table}[!h]
\\centering
\\begin{tabular}{lp{8cm}}
")

cat(paste("Number of miRNAs &", nrow(obj@cor) ,"\\\\")
)
cat(paste("Number of mRNAs &", ncol(obj@cor) ,"\\\\"))
cat(paste("Total miRNA-mRNA combinations &", nrow(obj@net) ,"\\\\"))
cat(paste("Number of samples &", length(obj@info$correlation.samples.used) ,"\\\\"))
if (length(obj@info$correlation.samples.used)<n.samp) {
cat(paste("Samples &", gsub("\\_","\\\\_",paste(obj@info$correlation.samples.used,collapse=", ")) ,"\\\\"))
}
cat("
\\end{tabular}
\\caption{Number of miRNAs, mRNAs and samples used for correlation.}
\\end{table}
")

n<-length(obj@info$correlation.samples.used)
t005<-qt(0.05,n-2)
cutoff005<-t005/sqrt(n+t005^2-2)
t001<-qt(0.01,n-2)
cutoff001<-t001/sqrt(n+t001^2-2)

sel<-which(obj@net$adj.pval<=0.05)
cutoff005.corrected<-max(obj@net[sel,"cor"])
sel<-which(obj@net$adj.pval<=0.01)
cutoff001.corrected<-max(obj@net[sel,"cor"])


pdf(paste("cordens",seed,".pdf",sep=""))

meth<-obj@info$correlation.type
if (obj@info$correlation.type=="pearson") meth<-"Pearson"
if (obj@info$correlation.type=="spearman") meth<-"Spearman"
if (obj@info$correlation.type=="kendall") meth<-"Kendall"

plot(density(obj@cor)$x,density(obj@cor)$y,xlab=paste(meth,"Correlation Coefficient",sep="\ "),ylab="Density",main="",type="l")
abline(v=cutoff005,lty=2)
abline(v=cutoff001,lty=3)
abline(v=cutoff005.corrected,lty=2,col="red")
abline(v=cutoff001.corrected,lty=3,col="red")
#legend("topright",c("p$<$0.05","p$<$0.01",""

#),col=c(1,1,2,2),lty=c(2,3,2,3))
dev.off()





cat("
\\renewcommand{\\arraystretch}{1.1}
\\definecolor{gray09}{gray}{0.9}
\\definecolor{gray099}{gray}{0.99}
\\rowcolors{1}{gray09}{gray099}
\\begin{table}[!h]
\\centering
\\begin{tabular}{lrr}
")


cat(" & Number & \\% ","\\\\")

cat(paste("Total correlations &", nrow(obj@net) ,  "&", round(nrow(obj@net)/nrow(obj@net)*100,2) ,"\\\\")
)

cat(paste("Total negative correlations &", length(which(obj@net$cor<0)) ,  "&", round(length(which(obj@net$cor<0))/nrow(obj@net)*100,2) ,"\\\\")
)

cat(paste("Total negative correlations p$<$0.05 &", length(which(obj@net$pval<0.05 & obj@net$cor<0)) ,  "&", round(length(which(obj@net$pval<0.05))/nrow(obj@net)*100,2) ,"\\\\")
)

cat(paste("Total negative correlations p$<$0.01 &", length(which(obj@net$pval<0.01& obj@net$cor<0)) ,  "&", round(length(which(obj@net$pval<0.01))/nrow(obj@net)*100,2) ,"\\\\")
)



cat(paste("Total negative correlations adj.p$<$0.05 &", length(which(obj@net$adj.pval<0.05& obj@net$cor<0)) ,  "&", round(length(which(obj@net$adj.pval<0.05))/nrow(obj@net)*100,2) ,"\\\\")
)

cat(paste("Total negative correlations adj.p$<$0.01 &", length(which(obj@net$adj.pval<0.01& obj@net$cor<0)) ,  "&", round(length(which(obj@net$adj.pval<0.01))/nrow(obj@net)*100,2) ,"\\\\")
)




cat("
\\end{tabular}
")
cat(paste("\\caption{Basic statistics for correlation results. Correlation hypothesis: ",obj@info[["cor.alternative.hypothesis"]],".}",sep=""))
cat("\\end{table}
")



cat("
\\begin{figure}[!h]
\\centering
")
cat(paste("\\includegraphics[width=0.6\\textwidth]{",paste("cordens",seed,".pdf",sep=""),"}",sep=""))
cat(paste("\\caption{Density of a total of ",nrow(obj@net)," miRNA-mRNA pairs. Dashed lines distinguish correlations whose p-value is lower than 0.05, dotted lines for 0.01. Black is for raw p-value and red for adjusted p-value.}",sep=""))
cat("
\\end{figure}
")

if (is.null(dat.sum.table)) {
  dat.sum.table<-obj@info[["dat.sum"]]
}

# top 11 correlations
n<-15
cyto<-obj@net
top11<-cyto[with(cyto, order(cor)),]
top11<-top11[which(top11$dat.sum>=dat.sum.table)[1:n],]
top11$miRNA<-as.character(top11$miRNA)
top11$mRNA<-as.character(top11$mRNA)

for (i in 1:n) {
	pdf(paste("cor",i,seed,".pdf",sep=""))
	plotCorrelation(obj,miRNA=top11$miRNA[i],mRNA=top11$mRNA[i])
	dev.off()

}


a<-paste(paste("\\includegraphics[width=0.3\\textwidth]{cor",1:n,seed,".pdf}",sep=""),collapse="\n")

cat("
\\begin{figure}[!h]
\\centering
")
cat(a)
cat(paste("\\caption{Plot of ",n," top correlations, sorted by adjusted p-value. Databases used: ",gsub("_","\\\\_",paste(obj@info[["database"]],collapse=", "))," (each miRNA-mRNA pair has to appear at least ",dat.sum.table," times).}",sep=""))
cat("
\\end{figure}
")



#cat("\\section{Intersection with databases}")


## venn
#require(VennDiagram)
obj@info[["cor_criteria"]]["pval"]<-0.05

corrs<-rownames(obj@net)[which(obj@net$adj.pval<=as.numeric(obj@info[["cor_criteria"]]["pval"]))]
corrs<-rownames(obj@net)[which(obj@net$adj.pval<=0.05)]
tarrs<-rownames(obj@net)[which(obj@net$dat.sum>=obj@info[["dat.sum"]])]

    venn.diagram(list(correlation = corrs, targets = tarrs),fill = c("red", "green"),
  alpha = c(0.5, 0.5), cex = 2,cat.fontface = 4,lty =2, fontfamily =3, 
   filename = paste("venn",seed,".tiff",sep=""))


system(paste("convert venn",seed,".tiff venn",seed,".pdf",sep=""))


cat("
\\begin{figure}[h]
\\centering
")
cat(paste("\\includegraphics[width=0.35\\textwidth]{",paste("venn",seed,".pdf",sep=""),"}",sep=""))

if (obj@info[["cor.alternative.hypothesis"]]!="both") {
#cat(paste("\\caption{Venn Diagram. Database(s) selected: ",gsub("_","\\\\_",paste(obj@info[["database"]],collapse=", "))," (minimum coincidences across databases: ",obj@info[["dat.sum"]],"), Pval-adjusted cutoff: ",obj@info[["cor_criteria"]]["pval"],"}",sep=""))
cat(paste("\\caption{Venn Diagram. Left (red): number of miRNA-mRNA pairs with ajdusted p-value$<$0.05. Right (green): number of all the theoretical miRNA-mRNA pairs reported at least ",obj@info[["dat.sum"]]," times in the following databases: ", gsub("_","\\\\_",paste(obj@info[["database"]],collapse=", ")),". Intersection: miRNA-mRNA pairs that fulfill both conditions.}",sep=""))
	} else {
cat(paste("\\caption{Venn Diagram. Left (red): number of miRNA-mRNA pairs with ajdusted p-value$<$0.05. Right (green): number of all the theoretical miRNA-mRNA pairs reported at least ",obj@info[["dat.sum"]]," times in the following databases: ", gsub("_","\\\\_",paste(obj@info[["database"]],collapse=", ")),". Intersection: miRNA-mRNA pairs that fulfill both conditions.}",sep=""))
	}
cat("
\\end{figure}
")


if (is.null(dat.sum.table)) {
  dat.sum.table<-obj@info[["dat.sum"]]
}


pairs.good<-obj@net[which(obj@net$dat.sum>=dat.sum.table),]
pairs.good<-pairs.good[order(pairs.good$pval),]


##### The content of the table will vary depending on which available information we have
if (!is.null(pairs.good$logratio.miRNA) & !is.null(pairs.good$logratio.mRNA)) {
pairs.good$FC.miRNA<-logratio2foldchange(pairs.good$logratio.miRNA)
pairs.good$FC.mRNA<-logratio2foldchange(pairs.good$logratio.mRNA)

	print(xtable(pairs.good[1:45,c("miRNA","mRNA","cor","adj.pval","FC.miRNA","FC.mRNA","dat.sum")],caption=paste("Top 45 miRNA-mRNA pairs (sorted by adjusted p-value) that have: pval-corrected$<$0.05 and appear at least ",dat.sum.table," times in the following databases: ",gsub("_","\\\\_",paste(obj@info[["database"]],collapse=", ")),".",sep=""),align="rrrrrrrr",display=c("s","s","s","f","e","f","f","d")),table.placement="!h",
include.rownames=FALSE,latex.environments="center")
}



if (!is.null(pairs.good$logratio.miRNA) & is.null(pairs.good$logratio.mRNA)) {
#pairs.good$FC.miRNA<-logratio2foldchange(pairs.good$logratio.miRNA)
pairs.good$FC.mRNA<-logratio2foldchange(pairs.good$logratio.mRNA)

	print(xtable(pairs.good[1:45,c("miRNA","mRNA","cor","adj.pval","FC.miRNA","dat.sum")],caption=paste("Top 45 miRNA-mRNA pairs(sorted by adjusted p-value) that have: pval-corrected$<$0.05 and appear at least ",dat.sum.table," times in the following databases: ",gsub("_","\\\\_",paste(obj@info[["database"]],collapse=", ")),".",sep=""),align="rrrrrrr",display=c("s","s","s","f","e","f","d")),table.placement="!h",
include.rownames=FALSE,latex.environments="center")
}



if (is.null(pairs.good$logratio.miRNA) & !is.null(pairs.good$logratio.mRNA)) {
pairs.good$FC.miRNA<-logratio2foldchange(pairs.good$logratio.miRNA)
#pairs.good$FC.mRNA<-logratio2foldchange(pairs.good$logratio.mRNA)

	print(xtable(pairs.good[1:45,c("miRNA","mRNA","cor","adj.pval","FC.mRNA","dat.sum")],caption=paste("Top 45 miRNA-mRNA pairs(sorted by adjusted p-value) that have: pval-corrected$<$0.05 and appear at least ",dat.sum.table," times in the following databases: ",gsub("_","\\\\_",paste(obj@info[["database"]],collapse=", ")),".",sep=""),align="rrrrrrr",display=c("s","s","s","f","e","f","d")),table.placement="!h",
include.rownames=FALSE,latex.environments="center")
}



if (is.null(pairs.good$logratio.miRNA) & is.null(pairs.good$logratio.mRNA)) {
#pairs.good$FC.miRNA<-logratio2foldchange(pairs.good$logratio.miRNA)
#pairs.good$FC.mRNA<-logratio2foldchange(pairs.good$logratio.mRNA)

	print(xtable(pairs.good[1:45,c("miRNA","mRNA","cor","adj.pval","dat.sum")],caption=paste("Top 45 miRNA-mRNA pairs(sorted by adjusted p-value) that have: pval-corrected$<$0.05 and appear at least ",dat.sum.table," times in the following databases: ",gsub("_","\\\\_",paste(obj@info[["database"]],collapse=", ")),".",sep=""),align="rrrrrr",display=c("s","s","s","f","e","d")),table.placement="!h",
include.rownames=FALSE,latex.environments="center")

}










pdf(paste("circos",seed,".pdf",sep=""))
plotCircos(obj,pval.cutoff=1,n=45)
dev.off()

cat("
\\begin{figure}[!h]
\\centering
")
cat(paste("\\includegraphics[width=0.75\\textwidth]{",paste("circos",seed,".pdf",sep=""),"}",sep=""))

cat(paste("\\caption{Circos plot for the first 45 miRNA-mRNA pairs (sorted by adjusted p-value) that have: pval-corrected$<$0.05 and appear at least ",obj@info[["dat.sum"]]," times in the following databases: ",gsub("_","\\\\_",paste(obj@info[["database"]],collapse=", ")),". Blue: miRNAs, Orange: target mRNAs}",sep=""))
cat("
\\end{figure}
")





cat("\\clearpage")


cat("\\section{Functional analysis}")
cat("\\subsection{Network analysis}")


if (length(which(obj@net$adj.pval<0.05 & obj@net$dat.sum>=obj@info[["dat.sum"]]))>0 ) {

pdf(paste("network",seed,".pdf",sep=""))
plotNetwork(obj,names=FALSE,pval.cutoff=0.05,node.size=1)
dev.off()

} else {

pdf(paste("network",seed,".pdf",sep=""))
plotNetwork(obj,names=FALSE,n=45,pval.cutoff=1,node.size=1.5)
dev.off()


}



cat("
\\begin{figure}[!h]
\\centering
")
cat(paste("\\includegraphics[width=0.9\\textwidth]{",paste("network",seed,".pdf",sep=""),"}",sep=""))

cat(paste("\\caption{Network for all the miRNA-mRNA pairs that have: pval-corrected$<$0.05 and appear at least ",obj@info[["dat.sum"]]," times in the following databases: ",gsub("_","\\\\_",paste(obj@info[["database"]],collapse=", ")),". Circles represent the miRNAs, and squares the mRNAs. Red fill means upregulated miRNAs/mRNAs, while green fill means downregulated miRNAs/mRNAs in comparative ",obj@info[["miRNA.diffexp.method"]][2],"; lines indicate the miRNA-mRNA pairs, red line means positive score and green line means negative score.}",sep=""))
cat("
\\end{figure}
")


if (length(which(obj@net$adj.pval<0.05 & obj@net$dat.sum>=obj@info[["dat.sum"]]))!=0) {


pdf(paste("barplot_miRNA",seed,".pdf",sep=""),height=7/2)
topTable(obj,"miRNA",names=TRUE,plot=TRUE)
dev.off()

cat("
\\begin{figure}[!h]
\\centering
")
cat(paste("\\includegraphics[width=0.9\\textwidth]{",paste("barplot_miRNA",seed,".pdf",sep=""),"}",sep=""))

cat(paste("\\caption{Barplot showing the number of mRNA targets per each miRNA (each bar represents a miRNA and they are sorted by number of targets). MiRNA-mRNA interactions have pval-corrected$<$0.05 and predicted at least ",obj@info[["dat.sum"]]," time on the following databases: ",gsub("_","\\\\_",paste(obj@info[["database"]],collapse=", ")),". Red line (and right axis) represents the percentage of deregulated mRNAs that are cumulatively targeted by the miRNAs.}",sep=""))
cat("
\\end{figure}
")


topmiRs<-topTable(obj,"miRNA",names=TRUE)[1:10,]
topmiRs$names<-as.character(topmiRs$names)

cat("
\\renewcommand{\\arraystretch}{1.1}
\\definecolor{gray09}{gray}{0.9}
\\definecolor{gray099}{gray}{0.99}
\\rowcolors{1}{gray09}{gray099}
\\begin{table}[!h]
\\centering
\\begin{tabular}{lrrp{9cm}}
")

cat("miRNA & \\#targets & cum. \\% & targets (top 20)\\\\")


for (i in 1:10) {
list.targ<-strsplit(topmiRs[i,"names"],", ")[[1]]

if (!is.null(obj@net$logratio.miRNA)) {
	cor.targ<-obj@net[which(obj@net$miRNA==rownames(topmiRs)[i] & obj@net$mRNA %in% list.targ),c("mRNA","cor","logratio.miRNA")]
} else {
	cor.targ<-obj@net[which(obj@net$miRNA==rownames(topmiRs)[i] & obj@net$mRNA %in% list.targ),c("mRNA","cor")]
	cor.targ$logratio.miRNA<-0
}


if (dim(cor.targ)[1]>0) {

ordere<-cor.targ[order(cor.targ$cor),"mRNA"]


if (length(list.targ)>20) {
	list.targ<-paste(ordere[1:20],collapse=", ")
} else {list.targ<-paste(ordere,collapse=", ")}

miRNA.name <- rownames(topmiRs)[i]
miRNA.name <- paste("\\textbf{",miRNA.name,"}",sep="")

if (cor.targ$logratio.miRNA[1]>0) {
miRNA.name <- paste("\\textcolor{red}{",miRNA.name,"}",sep="")
}
if (cor.targ$logratio.miRNA[1]<0) {
miRNA.name <- paste("\\textcolor{green}{",miRNA.name,"}",sep="")
}
if (cor.targ$logratio.miRNA[1]==0) {
miRNA.name <- paste("\\textcolor{black}{",miRNA.name,"}",sep="")
}


cat(paste(miRNA.name," & ",topmiRs[i,"freq"]," & ", round(topmiRs[i,"perc.regulated"],2), " & ", list.targ, "\\\\"   ,sep=""))

}

}

cat("
\\end{tabular}
\\caption{")
cat(paste("Top 10 miRNA with more targets (each miRNA-mRNA pair has pval-corrected$<$0.05 and appears at least ",obj@info[["dat.sum"]]," times in the following databases: ",gsub("_","\\\\_",paste(obj@info[["database"]],collapse=", ")),"). MiRNAs in red are upregulated in ",obj@info[["miRNA.diffexp.method"]][2],", miRNAs in green are downregulated in ",obj@info[["miRNA.diffexp.method"]][2],".",sep=""))
cat("}
\\end{table}
")





topmRs<-topTable(obj,"mRNA",names=TRUE)[1:10,]
topmRs$names<-as.character(topmRs$names)

cat("
\\renewcommand{\\arraystretch}{1.1}
\\definecolor{gray09}{gray}{0.9}
\\definecolor{gray099}{gray}{0.99}
\\rowcolors{1}{gray09}{gray099}
\\begin{table}[!h]
\\centering
\\begin{tabular}{lrp{12.5cm}}
")

cat("mRNA & \\#miRNAs  & miRNAs (top 20)\\\\")


for (i in 1:10) {
list.targ<-strsplit(topmRs[i,"names"],", ")[[1]]

if (!is.null(obj@net$logratio.mRNA)) {
	cor.targ<-obj@net[which(obj@net$mRNA==rownames(topmRs)[i] & obj@net$miRNA %in% list.targ),c("miRNA","cor","logratio.mRNA")]
} else {
	cor.targ<-obj@net[which(obj@net$mRNA==rownames(topmRs)[i] & obj@net$miRNA %in% list.targ),c("miRNA","cor")]
	cor.targ$logratio.mRNA<-0
}


if (dim(cor.targ)[1]>0) {

#cor.targ<-obj@net[which(obj@net$mRNA==rownames(topmRs)[i] & obj@net$miRNA %in% list.targ),c("miRNA","cor","logratio.mRNA")]
ordere<-cor.targ[order(cor.targ$cor),"miRNA"]

if (length(list.targ)>20) {
	list.targ<-paste(ordere[1:20],collapse=", ")
} else {list.targ<-paste(ordere,collapse=", ")}


mRNA.name <- rownames(topmRs)[i]
mRNA.name <- paste("\\textbf{",mRNA.name,"}",sep="")

if (cor.targ$logratio.mRNA[1]>0) {
mRNA.name <- paste("\\textcolor{red}{",mRNA.name,"}",sep="")
}
if (cor.targ$logratio.mRNA[1]<0) {
mRNA.name <- paste("\\textcolor{green}{",mRNA.name,"}",sep="")
}
if (cor.targ$logratio.mRNA[1]==0) {
mRNA.name <- paste("\\textcolor{black}{",mRNA.name,"}",sep="")
}


cat(paste(mRNA.name," & ",topmRs[i,"freq"], " & ", list.targ, "\\\\"   ,sep=""))

}

}



cat("
\\end{tabular}
\\caption{")
cat(paste("Top 10 mRNA with more miRNAs targeting them (each miRNA-mRNA pair has pval-corrected$<$0.05 and appears at least ",obj@info[["dat.sum"]]," times in the following databases: ",gsub("_","\\\\_",paste(obj@info[["database"]],collapse=", ")),"). MRNAs in red are upregulated in ",obj@info[["mRNA.diffexp.method"]][2],", mRNAs in green are downregulated in ",obj@info[["mRNA.diffexp.method"]][2],".",sep=""))
cat("}
\\end{table}
")


obj@info[["mRNA.diffexp.method"]][2]



#pdf(paste("barplot_mRNA",seed,".pdf",sep=""),height=7/2)
#topTable(obj,"mRNA",names=FALSE,plot=TRUE,remove.names=TRUE)
#dev.off()



# piechart
all<-topTable(obj,"mRNA",names=FALSE)
table(all)



pdf(paste("piechart_mRNA",seed,".pdf",sep=""))
a<-topTable(obj,"mRNA",names=FALSE)
a[which(a>5)]<-">5"


#pie(table(a),labels=paste(names(table(a))," (",table(a),", ",round(table(a)/sum(table(a))*100,1),"%)",sep=""),col=topo.colors(10))
pie(table(a),labels=paste(names(table(a))," (",table(a)," mRNAs)",sep=""),col=topo.colors(10))


dev.off()




#cat("
#\\begin{figure}[!h]
#\\centering
#")
#cat(paste("\\includegraphics[width=0.9\\textwidth]{",paste("barplot_mRNA",seed,".pdf",sep=""),"}",sep=""))

#cat(paste("\\caption{Barplot for mRNAs, pval-corrected$<$0.05 and Targets=",gsub("_","\\\\_",paste(obj@info[["database"]],collapse=", ")),"(minimum coincidences between databases:",obj@info[["dat.sum"]],").}",sep=""))
#cat("
#\\end{figure}
#")



cat("
\\begin{figure}[!h]
\\centering
")
cat(paste("\\includegraphics[width=0.7\\textwidth]{",paste("piechart_mRNA",seed,".pdf",sep=""),"}",sep=""))

cat(paste("\\caption{Pie chart representing the number of miRNAs targeting the mRNAs, pval-corrected$<$0.05 and Targets=",gsub("_","\\\\_",paste(obj@info[["database"]],collapse=", ")),"(minimum coincidences between databases:",obj@info[["dat.sum"]],").}",sep=""))
cat("
\\end{figure}
")


}


cat("\\vspace{2cm}")


#cat("\\newpage")
cat("\\clearpage")
cat("\\subsection{GO analysis}")

make.options <- function (ontology) {

options<-". Options used: mRNAs that are present in a mRNA-mRNA pair that has "
options<-paste(options," adjusted-pval cutoff <",obj@info[[ontology]][["pval.cutoff"]],sep="")
options<-paste(options,"; that also appears at least ",obj@info[[ontology]][["dat.sum"]]," times (databases: ",gsub("_","\\\\_",paste(obj@info[["database"]],collapse=", ")),")",sep="")
if (!is.null(obj@info[[ontology]][["miRNA"]])) {
	options<-paste(options,"; restricted to the following miRNAs: ",paste(obj@info[[ontology]][["miRNA"]],collapse=", "),sep="")
}
if (!is.null(obj@info[[ontology]][["mRNA"]])) {
	options<-paste(options,"; restricted to the following mRNAs: ",paste(obj@info[[ontology]][["mRNA"]],collapse=", "),sep="")
}
if (!is.null(obj@info[[ontology]][["exclude.miRNA"]])) {
	options<-paste(options,"; excluding the following miRNAs: ",paste(obj@info[[ontology]][["exclude.miRNA"]],collapse=", "),sep="")
}
if (!is.null(obj@info[[ontology]][["exclude.mRNA"]])) {
	options<-paste(options,"; excluding the following mRNAs: ",paste(obj@info[[ontology]][["exclude.mRNA"]],collapse=", "),sep="")
}
if (obj@info[[ontology]][["up"]]) {
	options<-paste(options,"; restricted to upregulated mRNAs",sep="")
}
if (obj@info[[ontology]][["dw"]]) {
	options<-paste(options,"; restricted to downregulated mRNAs",sep="")
}
options<-paste(options,"; organism: ",obj@info[[ontology]][["organism"]],".",sep="")

	return(options)
}



n.func<-10
if (!is.null(obj@GO.results[["GO:BP"]])) {
	opt<-make.options("GO:BP")
	print(xtable(obj@GO.results[["GO:BP"]][1:n.func,c(2,8,6,7,5,4,9,3)],caption=paste("Biological Process",opt,collapse=". "),align="lrp{4cm}rrrrrr",display=c("s","s","s","d","d","f","f","e","e")),include.rownames=FALSE,latex.environments="center")
}

if (!is.null(obj@GO.results[["GO:CC"]])) {
	opt<-make.options("GO:CC")
	print(xtable(obj@GO.results[["GO:CC"]][1:n.func,c(2,8,6,7,5,4,9,3)],caption=paste("Cellular Component",opt,collapse=". "),align="lrp{4cm}rrrrrr",display=c("s","s","s","d","d","f","f","e","e")),include.rownames=FALSE,latex.environments="center")
}

if (!is.null(obj@GO.results[["GO:MF"]])) {
	opt<-make.options("GO:MF")
	print(xtable(obj@GO.results[["GO:MF"]][1:n.func,c(2,8,6,7,5,4,9,3)],caption=paste("Molecular Function",opt,collapse=". "),align="lrp{4cm}rrrrrr",display=c("s","s","s","d","d","f","f","e","e")),include.rownames=FALSE,latex.environments="center")
}

if (!is.null(obj@GO.results[["KEGG:KEGG"]])) {
	opt<-make.options("KEGG:KEGG")
	print(xtable(obj@GO.results[["KEGG:KEGG"]][1:n.func,c(2,8,6,7,5,4,9,3)],caption=paste("Kegg Pathways",opt,collapse=". "),align="lrp{4cm}rrrrrr",display=c("s","s","s","d","d","f","f","e","e")),include.rownames=FALSE,latex.environments="center")
}


#align
#display


#ltable = sub(sprintf("\\begin{%s}", TABULAR),
#             sprintf("\\rowcolors{2}{gray!25}{white}\n\\begin{%s}", TABULAR),
#             ltable, fixed=TRUE)






cat("\\end{document}")


	sink()

	system(paste("pdflatex ",file,".tex",sep=""))
	system(paste("pdflatex ",file,".tex",sep=""))
	system(paste("rm ",file,".out",sep=""))
	system(paste("rm ",file,".aux",sep=""))
	system(paste("rm ",file,".log",sep=""))

	system(paste("rm *",seed,".pdf",sep=""))
	system(paste("rm *",seed,".tiff",sep=""))
	#system(paste("rm *",seed,".jpg",sep=""))
	

	#system("rm mamiRNA.pdf")
	#system("rm mamRNA.pdf")


}




