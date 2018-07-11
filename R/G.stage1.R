#!/usr/bin/Rscript

library("Deducer")

####G 

 argv<-commandArgs(TRUE);
 
 x<-read.table(argv[1]);

 s<-data.frame(x[1],x[2],x[5],x[6],x[7],x[8]);

 y<-s[s[1] == argv[2],];
  rm(x,s);
  colnames(y)<-c("chr","pos","R01_REF","R01_ALT","R02_REF","R02_ALT")
   G<-0;  
 for(i in 1:nrow(y)){
      m<-matrix(c(y[i,3],y[i,5],y[i,4],y[i,6]),nc=2);
       t<-likelihood.test(m) 
       G[i]<-as.data.frame(t$statistic)[1,1]
        } 

  z<-data.frame(y,G);
    
  
   name<-paste(argv[2],"stage1.Rdata",sep=".");
   out<-as.character(name);
  data<-data.frame(y,G)
  save(data,file=name)       


