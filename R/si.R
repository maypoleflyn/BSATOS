#!/usr/bin/Rscript 

 argv <- commandArgs(trailingOnly=T)
 cat("Importing files.\n")
 
 data <- read.table(argv[1],header=F);
 
 colnames(data)<-c("chr","pos","R01_REF","R01_ALT","R02_REF","R02_ALT")
 chr<-argv[2];
 data<-data[data[1]==chr,]
 data$R01_REF_F<-(data$R01_REF)/(data$R01_REF + data$R01_ALT);
 data$R01_ALT_F<-(data$R01_ALT)/(data$R01_REF + data$R01_ALT);
 data$R02_REF_F<-(data$R02_REF)/(data$R02_REF + data$R02_ALT);
 data$R02_ALT_F<-(data$R02_ALT)/(data$R02_REF + data$R02_ALT);

 data$si<-abs(data$R01_ALT_F-data$R02_ALT_F);
 
   

ke<-function(x,w){
          H<-NULL;
   for(i in 1:nrow(x)){
           min<-x[i,2]-w/2;
           max<-x[i,2]+w/2;
          subset<-x[x[,2]>min & x[,2]<max,];
          D<-NULL;
          K<-NULL;
          A<-NULL;           
        
           D=abs(subset[2]-x[i,2])/w;
           A<-sum((1-D^3)^3) 
           K<-(1-(D)^3)^3/A;  
           n<-subset[11]*K;
           H[i]<-sum(n)
         
     }
                 return(H);
    }

    G1M<-ke(data,1000000);
    G750<-ke(data,750000);
    G500<-ke(data,500000);
    G250<-ke(data,250000);
 
   name=paste(argv[2],argv[3],"afd",sep="_");
   out<-as.character(name); 
   data<-data.frame(data,S1M,S750,S500,S250)
   write.table(data,file=name,quote =F,sep = "\t",row.names = F,col.names = F) 
  


