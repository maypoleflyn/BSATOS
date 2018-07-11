#!/usr/bin/Rscript



    
#####Gp


  argv<-commandArgs(TRUE);
  load(file=argv[1]);     
  
  
ke<-function(x){
          H<-0;
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
         
     }
                 return(H);
    }      
    
  Gp<-ke(data)  
  
   name<-paste(argv[2],"stage2.Rdata",sep=".");
   out<-as.character(name);
  data<-data.frame(data,Gp)
  save(data,file=name)       


