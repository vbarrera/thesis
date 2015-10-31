#################
### Figure 26 ###
#################


WORKINGDIR<-paste(Sys.getenv("HOME"),"/projects/Evaluation_of_single_CpG_sites_as_proxies_of_CpG_island_methylation_states_at_the_genome_scale/oldOrganization/Cstatus/output/",sep="")
CGIH1=paste(WORKINGDIR,"h1/CGi/",sep="")
CGIIMR90=paste(WORKINGDIR,"imr90/CGi/",sep="")

set.seed(2012)

x<-read.table(paste(CGIH1,"cgiCGListH1_5READS.txt",sep=""),sep=";")
cgTH1<-subset(x,x$V2!=-1)

z<-read.table(paste(CGIH1,"methMaps/methMap_hg18_H1_covered5READ_usefulCGi.txt",sep=""),sep="\t")
cgiTH1<-subset(z,z$V3 != -1)

x<-read.table(paste(CGIIMR90,"cgiCGListIMR90_5READS.txt",sep=""),sep=";")
cgTIMR90<-subset(x,x$V2!=-1)

z<-read.table(paste(CGIIMR90,"methMaps/methMap_hg18_IMR90_covered5READ_usefulCGi.txt",sep=""),sep="\t")
cgiTIMR90<-subset(z,z$V3 != -1)


randomCGi<-function(cgiTable,cgTable) {
  usefulCG<-cgiTable$V2
  counter=1
  results<-data.frame()
  for (cgInf in usefulCG){
    cgList<-cgTable[round(runif(cgInf,1,dim(cgTable)[1]),0),2]
    restrEnz<-selectRE(cgList)
    meanCG<-obtainMean(cgList)
    sdCG<-obtainSD(cgList)
    results[counter,1]<-restrEnz
    results[counter,2]<-meanCG
    results[counter,3]<-sdCG
    counter=counter+1
  }
  return (results)
}

selectRE<-function(cgList){
  selectedRE<-cgList[round(runif(1,1,length(cgList)),2)]
  return (selectedRE)
}

obtainMean<-function (cgList){
  return (mean(cgList))
}

obtainSD<-function (cgList){
  return (sd(cgList))
}


## H1

png(paste(resultsDIR,"figure26H1.png",sep=""),units="cm",height=15,width=15,res=300)

def.par <- par(no.readonly = TRUE)

par(lwd=3)
par(cex.axis=1.2)

randomH1_1<-randomCGi(cgiTH1,cgTH1)

colnames(randomH1_1)<-c("RE","mean","sd")

def.par <- par(no.readonly = TRUE)
xhist <- hist(randomH1_1$RE, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(randomH1_1$mean, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))

plot (jitter(randomH1_1$RE,factor=25),jitter(randomH1_1$mean,factor=25),pch=16,col=rgb(0,0,0,0.2),xlab="",
      ylab="",main="",cex=1)

# plot (randomH1_1$RE,randomH1_1$mean,pch=16,col=rgb(1,0,0,0.2),xlab="",
#        ylab="",main="",cex=1)


lmRndH1_1<-lm(randomH1_1$mean~randomH1_1$RE)
abline(lmRndH1_1,lty=2,lwd=1.5,col="red")

par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,
        xlab="")
par(mar=c(5,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0,
        horiz=TRUE)

par(def.par)

dev.off()

## IMR90

png(paste(resultsDIR,"figure26IMR90.png",sep=""),units="cm",height=15,width=15,res=300)

def.par <- par(no.readonly = TRUE)

par(lwd=3)
par(cex.axis=1.2)

randomIMR90_1<-randomCGi(cgiTIMR90,cgTIMR90)
colnames(randomIMR90_1)<-c("RE","mean","sd")

def.par <- par(no.readonly = TRUE)
xhist <- hist(randomIMR90_1$RE, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(randomIMR90_1$mean, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))

plot (jitter(randomIMR90_1$RE,factor=25),jitter(randomIMR90_1$mean,factor=25),pch=16,col=rgb(0,0,0,0.2),xlab="",
      ylab="",main="",cex=1)
# plot (randomIMR90_1$RE,randomIMR90_1$mean,pch=16,col=rgb(1,0,0,0.2),xlab="",
#       ylab="",main="",cex=1)


lmRndIMR90_1<-lm(randomIMR90_1$mean~randomIMR90_1$RE)
abline(lmRndIMR90_1,lty=2,lwd=1.5,col="red")

par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,
        xlab="SD")
par(mar=c(5,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0,
        horiz=TRUE)

par(def.par)

dev.off()