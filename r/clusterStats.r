
# reading performance .csv files
lfiles <- list.files(path="../csv/fieldData/", pattern="Performance_Table")

# loop over files
for (i in 1:length(lfiles))
{
   drt<-strsplit(lfiles[[i]],"_")[[1]][2]
   clData <- read.csv(paste("../csv/fieldData/",lfiles[[i]],sep=""),header = TRUE,sep = ",")
   
   # getting high- and low-performance
   hp <- which(clData$performance == 1)
   lp <- which(clData$performance == 0)
   
   # percentage
   perc <- c(length(hp)/length(clData$performance)*100,
             length(lp)/length(clData$performance)*100)
   
   # labels  
   lbls <- c("HP","LP")
   lbls <- paste(lbls, format(round(perc,2),nsmall=2))
   lbls <- paste(lbls,"%",sep="") # add %
   tit <- paste("Cluster performance DRT =",drt,sep=" ")
   
   # plot and save figure
   figname <- paste(paste("performance-DRT",drt,sep=""),".png",sep="")
   png(paste("../img/",figname,sep=""))
   pie(perc,labels=lbls,main=tit,col=c("green","red"))
   dev.off()
}

