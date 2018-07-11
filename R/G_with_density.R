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
    
#####Gp

ke<-function(x){
          H<-0;
          density<-0;
          all<-0;
   for(i in 1:nrow(x)){
           min<-x[i,2]-500000;
           max<-x[i,2]+500000;
          subset<-x[x[,2]>min & x[,2]<max,];
          
          D<-0;
          K<-0;
          A<-0;           
        
           D=abs(subset[2]-x[i,2])/1000000;
         
           A<-sum((1-D^3)^3);  
           K<-(1-(D)^3)^3/A;  
                             
         
           n<-subset[7]*K;
           H[i]<-sum(n);
           density[i]<-nrow(subset);
         
     }         
                all<-data.frame(H,density);  
                 return(all);
    }      
    
  Gp<-ke(z)  
  
   name<-paste(argv[2],"Rdata",sep=".");
   out<-as.character(name);
  data<-data.frame(y,G,Gp)
  save(data,file=name)       


