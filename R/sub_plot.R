#!/usr/bin/Rscript

####G 

 argv<-commandArgs(TRUE);

  nu<-as.numeric(argv[3]);                        
  n <-as.numeric(argv[4]);


 
 g1<-read.table(argv[1]);

 colnames(g1)<-c("chr","pos","R01_REF","R01_ALT","R02_REF","R02_ALT","R01_REF_F","R01_ALT_F","R02_REF_F","R02_ALT_F","G","GV","GM","G750","G500","G250");
 
 g1<-g1[g1$chr == argv[2],];
 
 g1<-g1[g1$pos >= nu & g1$pos <=n,];  
 
 
 max<-max(g1$GM,g1$G750,g1$G500,g1$G250) +5;
 name=paste(argv[5],"pdf",sep=".");
 main=paste(argv[5]," mutiple G value profile",sep="");
 out<-as.character(name);
 main<-as.character(main);
 pdf(file=out);
 plot(g1$pos,g1$GM,type="l",ylim=c(0,max),col="red",ylab="G'",xlab="positon",main=main,pch="*");
 lines(g1$pos,g1$G750,type="l",col="orange",pch="*");
 lines(g1$pos,g1$G500,type="l",col="yellow",pch="*");
 lines(g1$pos,g1$G250,type="l",col="brown",pch="*");
 legend("topleft", c("1M", "0.75M", "0.5M","0.25M"), col = c("red","orange","yellow","brown"),pch = "*", ncol = 4, cex = 0.8);

 dev.off();
     
 

 

 
      
