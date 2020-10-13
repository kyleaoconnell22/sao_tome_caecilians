#install.packages('HIest')
library('HIest')
#install.packages('gdata')
library("gdata")
setwd('~indir')

full_data <- read.table("HI.in.41")

dim(full_data) # looks like you have one row per individual and two columns per locus
 
full_data <- replace(full_data,full_data=="0",NA) #and R doesn’t read “0” as NA

freq_out <- read.table("freq_out_41.freq",header=TRUE)

dim(freq_out)

# convert the full dataset:
Allele.1 <- full_data[,seq(from=1,to=81,by=2)] #change to the number of loci in your dataset *2 -1
dim(Allele.1)
Allele.2 <- full_data[,seq(from=2,to=82,by=2)] #number of loci in your dataset *2 (same on next line -1)
dim(Allele.2)
names(Allele.1) <- names(Allele.2) <- freq_out[seq(from=1,to=81,by=2),1] # gotta have matching names

full <- numeric()

for(i in 1:nrow(Allele.1)){
  
  full <- rbind(full,rbind(Allele.1[i,],Allele.2[i,]))
  
}
write.csv(full,'Hi.in.41.csv')
full.HI = HIest(full,freq_out,"codominant",method = "SANN",surf=TRUE)
#full.HI.long = HIest(full,freq_out,"codominant",method = "SANN",iterations=500,surf=TRUE,startgrid=50)
full.HI
write.csv(full.HI,'full.HI.41.csv')
#write.csv(full.HI.long,'full.HI_long.txt')
#2 version is just the redo
# # optional plot
full.HI = read.csv('full.HI.41.csv')
head(full.HI)

pdf("full.HI.pdf",5,5)
plot(c(0,.5,1,0),c(0,1,0,0),type="l",xlab=expression(italic(S)),
     ylab=expression(italic(H[I])),lwd=2,cex.lab=1.5,cex.axis=1.5,bty="n")
points(full.HI$S,full.HI$H,cex=1.5,lwd=2)
axis(1,labels=FALSE,lwd=2);axis(2,labels=FALSE,lwd=2)
dev.off()

full.HIC <- HIC(full.HI)
dim(full.HIC)
# # calculate likelihoods for early generation hybrid classes
full.class <- HIclass(full,freq_out,"codominant")
write.csv(full.class,'full.class.csv')
# # compare classification with maximum likelihood estimates
full.comp <- HItest(full.class,full.HI)
write.csv(full.comp,'full.comp.csv')
table(full.comp$c1)
full.comp=read.csv('full.comp.csv')

# # all 41 are TRUE, meaning the best classification is at least 2 log-likelihood units
# # better than the next best

table(full.comp$c2)
# # 2 are TRUE, meaning the MLE S and H are within 2 log-likelihood units of the best
# # classification, i.e., the simple classification is rejected in all but 2 cases

table(full.comp$Best.class,full.comp$c2)
# # individuals were classified as F2-like (class 3) or backcross to CTS (class 4), but
# # only two of the F2's were credible 

full.comp[full.comp$c2,]
# # in only one case was the F2 classification a better fit (based on AIC) than the
# # continuous model.

full_withloc = read.csv('full.HI.locs.filt.csv',header=T)
#max(full_withloc$H)

#color by locality
pal=c("#97C2A5","#9583A8","#F0975B")
col = adjustcolor(col = pal,alpha.f = .5)

# # optional plot
full.HI = read.csv('full.HI.locs.csv',header=T)

#from genomic data excel file
full.HI=read.xls('genomic_data.xlsx')
full.HI=as.data.frame(full.HI)
#assign north to purple
full.HI <-c(lapply(full.HI, function(x) {
  gsub("ST","#9583A8", x)
  }))
#assign south to green
full.HI <- c(lapply(full.HI, function(x) {
  gsub("SE","#97C2A5", x)
}))
#assign admix to cream
full.HI <- c(lapply(full.HI, function(x) {
  gsub("admixed","#F0975B", x)
}))
full.HI
#morpho colors
#assign north to cream
full.HI <-c(lapply(full.HI, function(x) {
  gsub("yellow","#FFEEAA", x)
}))
#assign south to gray
full.HI <- c(lapply(full.HI, function(x) {
  gsub("flecked","#4D4D4D", x)
}))
full.HI                                      
full.HI$Morpho
#pdf("full.HI_mtDNA_color.pdf",5,5)
#pdf("full.HI_morpho_color.pdf",5,5)
pdf("full.HI_nucmaj_color.pdf",5,5)
#pdf("full.HI_pure-admix-loc.pdf",5,5)
#pdf("full.HI_pure-admix-ind.pdf",5,5)
plot(c(0,.5,1,0),c(0,1,0,0),type="l",xlab=expression(italic(S)),
     ylab=expression(italic(H[I])),lwd=2,cex.lab=1.5,cex.axis=1.5,bty="n")
points(full.HI$AdmixQ2.SE.,full.HI$H,cex=1.8,lwd=2,col=full.HI$Admixed)
#points(full.HI$S,full.HI$H,cex=1.8,lwd=2,pch = 23,col='black',bg = full.HI$Morpho)
#points(full.HI$S,full.HI$H,cex=1.8,lwd=2,pch = 23,col = full.HI$Morpho)
axis(1,labels=FALSE,lwd=2);axis(2,labels=FALSE,lwd=2)
dev.off()
?plot
#lat vs. H with color
pdf("latvsH.pdf",5,5)
plot(full.HI$Lat,full.HI$H,cex=1.8,lwd=2,pch = 23,bg =full.HI$Admixed)
axis(1,labels=FALSE,lwd=2);axis(2,labels=FALSE,lwd=2)
dev.off()
#lat vs. S with color
pdf('latvsS.pdf',5,5)
plot(full.HI$lat,full.HI$S,cex=1.8, cex.axis = 10,lwd=2,pch = 23,col='black',bg =full.HI$PvsA,cex.lab=1.5,cex.axis=1.5)
dev.off()
pdf('SvH.pdf',5,5)
plot(full.HI$S,full.HI$H,cex=1.8,lwd=2,pch = 23,col='black',bg =full.HI$PvsA, cex.lab=1.5,cex.axis=1.5)
dev.off()


####BARPLOTS
#n77_maf5
setwd('~indir')
col=c("#9583A8","#97C2A5")
tbl2=read.table("k2_average_HIest.txt")
#tbl2=read.table('admixture_values.txt')
samples=tbl2$V1
### Need this function
barNaming <- function(vec) {
  retVec <- vec
  for(k in 2:length(vec)) {
    if(vec[k-1] == vec[k])
      retVec[k] <- ""
  }
  return(retVec)
}
pdf("n75_k2_HIest.pdf",8,8)
barplot(t(as.matrix(tbl2[2:3])), col=col, space=0, names.arg=barNaming(samples), las=2, cex.names = 0.5)
dev.off()
