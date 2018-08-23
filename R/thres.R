#!usr/bin/Rscript

 argv<-commandArgs(TRUE);
 
  x<-read.table(argv[1]);
  y1<-apply(x[12],2,median)+2*apply(x[12],2,sd)
  y2<-apply(x[13],2,median)+2*apply(x[13],2,sd)
  y3<-apply(x[14],2,median)+2*apply(x[14],2,sd)
  y4<-apply(x[15],2,median)+2*apply(x[15],2,sd)
  cat(y1,sep="\n");
  cat(y2,sep="\n");
  cat(y3,sep="\n");
  cat(y4,sep="\n");  
