options(bitmapType='cairo')

args <- commandArgs(TRUE)

dir<-args[1]
inputFile<-args[2]
headerSize<-as.numeric(args[3])
headerSpacing<-as.numeric(args[4])
imageWidth<-as.numeric(args[5])
imageHeight<-ceiling(imageWidth/8)
yBound<-as.numeric(args[6])
transparentBGFlag<-as.numeric(args[7])

setwd(dir)

directionData<-read.table(inputFile,header=T,sep="\t")

xStart<-min(directionData$start)
xEnd<-max(directionData$end)
xRange<-(xEnd-xStart)
xBinStart<-min(directionData$binStart)
xBinEnd<-max(directionData$binEnd)
xBinRange<-(xBinEnd-xBinStart)

directionData<-subset(directionData,directionData$directionScore!="NA") # remove NAs from the scores

directionData<-directionData[with(directionData, order(start)), ]

x<-directionData$binMidpoint
y<-directionData$directionScore

yStart<-min(y)
yEnd<-max(y)

if(yBound == 0) {
	yBound<-ceiling(max(abs(yStart),abs(yEnd)))
	if(yBound < 1) {	
		yBound<-1
	}
}

yStart<--yBound
yEnd<-yBound
yRange<-(yEnd-yStart)

pos_yBound <- yBound
neg_yBound <- -yBound

# saturate plot at specified yBound(s)
y[y > pos_yBound] <- pos_yBound
y[y < neg_yBound] <- neg_yBound

# set minimum image height
if(imageHeight < 325) {
	imageHeight<-325
}

# png file
pngfile<-paste(inputFile,".png",sep='')
if(transparentBGFlag == 1) {
	png(pngfile,height=imageHeight,width=imageWidth,bg="transparent")
} else {
	png(pngfile,height=imageHeight,width=imageWidth)
}

par(mar=c(4, 4, 4, 4) + 0.1)

plot(x,y,main=paste(inputFile,sep=""),cex=0.5,col=rgb(0.25,0.25,0.25,0.5),xlab="genomic coordinates",ylab="directionality index",xlim=c(xBinStart,xBinEnd),ylim=c(-yBound,yBound),type="n",xaxt="n",yaxt="n")
axis(2,seq(from=-yBound,to=yBound,length.out=11))

if(length(x) > 1) {

	for (i in 1:(length(x)-1) ) {
		tmpX<-x[i]
		tmpY<-y[i]
		nextX<-x[i+1]
		nextY<-y[i+1]
		
		if((i == 1) && ((tmpX-xBinStart) != 0)) {
			rect(xBinStart, yStart, tmpX, yEnd, col=rgb(0.75,0.75,0.75,0.5),border=FALSE,lwd=0)
		}
		
		if((nextX-tmpX) == 1) {
			segments(tmpX,tmpY,nextX,nextY,col="black",lwd=1)
		} else {
			rect(tmpX, yStart, nextX, yEnd, col=rgb(0.75,0.75,0.75,0.5),border=FALSE,lwd=1)
		}
		
		if((i == (length(x)-1)) && ((xBinEnd-tmpX) != 0)) {
			rect(tmpX, yStart, xBinEnd, yEnd, col=rgb(0.75,0.75,0.75,0.5),border=FALSE,lwd=1)
		}
	}
	
	# draw saturation blobs
	for (i in 1:(length(x)-1) ) {
		tmpX<-x[i]
		tmpY<-y[i]
		nextX<-x[i+1]
		nextY<-y[i+1]
		
		if((nextX-tmpX) == 1) {
			if(nextY >= yBound) {
				segments(tmpX,yBound,nextX,yBound,col="purple",lwd=4)
			} else if(nextY <= -yBound) {
				segments(tmpX,-yBound,nextX,-yBound,col="purple",lwd=4)
			}
		}
		
	}

	abline(h=0,lwd=1,lty=2,col="black")
	abline(v=xBinStart,col="red",lwd=2,lty=2)
	abline(v=xBinEnd,col="red",lwd=2,lty=2)
	
}

axis(1,at=seq(from=xBinStart,to=xBinEnd,length.out=11),labels=seq(from=xStart,to=xEnd,length.out=11))

dev.off()



# pdf file
pdffile<-paste(inputFile,".pdf",sep='')
pdf(pdffile)

par(mar=c(4, 4, 4, 4) + 0.1)

plot(x,y,main=paste(inputFile,sep=""),cex=0.5,col=rgb(0.25,0.25,0.25,0.5),xlab="genomic coordinates",ylab="directionality index",xlim=c(xBinStart,xBinEnd),ylim=c(-yBound,yBound),type="n",xaxt="n",yaxt="n")
axis(2,seq(from=-yBound,to=yBound,length.out=11))

if(length(x) > 1) {

	for (i in 1:(length(x)-1) ) {
		tmpX<-x[i]
		tmpY<-y[i]
		nextX<-x[i+1]
		nextY<-y[i+1]
		
		if((i == 1) && ((tmpX-xBinStart) != 0)) {
			rect(xBinStart, yStart, tmpX, yEnd, col=rgb(0.75,0.75,0.75,0.5),border=FALSE,lwd=0)
		}
		
		if((nextX-tmpX) == 1) {
			segments(tmpX,tmpY,nextX,nextY,col="black",lwd=1)
		} else {
			rect(tmpX, yStart, nextX, yEnd, col=rgb(0.75,0.75,0.75,0.5),border=FALSE,lwd=1)
		}
		
		if((i == (length(x)-1)) && ((xBinEnd-tmpX) != 0)) {
			rect(tmpX, yStart, xBinEnd, yEnd, col=rgb(0.75,0.75,0.75,0.5),border=FALSE,lwd=1)
		}
	}
	
	# draw saturation blobs
	for (i in 1:(length(x)-1) ) {
		tmpX<-x[i]
		tmpY<-y[i]
		nextX<-x[i+1]
		nextY<-y[i+1]
		
		if((nextX-tmpX) == 1) {
			if(nextY >= yBound) {
				segments(tmpX,yBound,nextX,yBound,col="purple",lwd=4)
			} else if(nextY <= -yBound) {
				segments(tmpX,-yBound,nextX,-yBound,col="purple",lwd=4)
			}
		}
		
	}

	abline(h=0,lwd=1,lty=2,col="black")
	abline(v=xBinStart,col="red",lwd=2,lty=2)
	abline(v=xBinEnd,col="red",lwd=2,lty=2)
	
}

axis(1,at=seq(from=xBinStart,to=xBinEnd,length.out=11),labels=seq(from=xStart,to=xEnd,length.out=11))

dev.off()