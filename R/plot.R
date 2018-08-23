#!/usr/bin/Rscript

####G 

 argv<-commandArgs(TRUE);
 
 g1<-read.table(argv[1]);
 colnames(g1)<-c("chr","pos","R01_REF","R01_ALT","R02_REF","R02_ALT","R01_REF_F","R01_ALT_F","R02_REF_F","R02_ALT_F","G","GV","GM","G750","G500","G250");
 g2<-read.table(argv[2]);
 colnames(g2)<-c("chr","pos","R01_REF","R01_ALT","R02_REF","R02_ALT","R01_REF_F","R01_ALT_F","R02_REF_F","R02_ALT_F","G","GV","GM","G750","G500","G250");
 g3<-read.table(argv[3]);
 colnames(g3)<-c("chr","pos","R01_REF","R01_ALT","R02_REF","R02_ALT","R01_REF_F","R01_ALT_F","R02_REF_F","R02_ALT_F","G","GV","GM","G750","G500","G250");
 
 g1<-g1[g1$chr == argv[4],];
 g2<-g2[g2$chr == argv[4],];
 g3<-g3[g3$chr == argv[4],];
 
 max<-max(g1$GM,g2$GM,g3$GM) +7;
 name=paste(argv[4],"pdf",sep=".");
 main=paste(argv[4]," smoothed G value profile",sep="");
 out<-as.character(name);
 main<-as.character(main);
 pdf(file=out,width=14,height=7);
 plot(g1$pos,g1$GM,type="l",ylim=c(0,max),col="blue",ylab="G'",xlab="positon",main=main,lwd=3);
 points(g1$pos,g1$GV,cex=0.1,col="blue");
 lines(g3$pos,g3$GM,type="l",col="black",lwd=3);
 points(g3$pos,g3$GV,cex=0.1,col="black");
 lines(g2$pos,g2$GM,type="l",col="red",lwd=3);
 points(g2$pos,g2$GV,cex=0.1,col="red");
 abline(h=argv[5],col="blue",lwd=3,lty=3);
 abline(h=argv[6],col="red",lwd=3,lty=3);
 abline(h=argv[7],col="black",lwd=3,lty=3);
 dev.off();
     
 

 

 
      
