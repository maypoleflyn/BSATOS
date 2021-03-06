#!/usr/bin/Rscript
 
   argv<-commandArgs(TRUE);

   data<-read.table(argv[1]);
   data<-data[,c(1:16)];
   colnames(data)<-c("chr","pos","R01_REF","R01_ALT","R02_REF","R02_ALT","R01_REF_F","R01_ALT_F","R02_REF_F","R02_ALT_F","G","GV","GM","G750","G500","G250");

     nu<-as.numeric(argv[2]);                        
     n <-as.numeric(argv[3]);

        
   sub<-data[data$pos >=nu & data$pos <= n,];
     
   max<-sub[which.max(sub$GM),];   
   max$GM;
   max$pos;



   
