#!/usr/bin/Rscript

 argv<-commandArgs(TRUE); 
 d<-read.table(argv[1]);
 d<-d[,c(1:16)];
 colnames(d)<-c("chr","pos","R01_REF","R01_ALT","R02_REF","R02_ALT","R01_REF_F","R01_ALT_F","R02_REF_F","R02_ALT_F","G","GV","GM","G750","G500","G250");
  E<-0;
  x<-d[d$chr == argv[2],]; 
  nu<-as.numeric(argv[3]);                        
  n <-as.numeric(argv[4]);
for(i in 1:nrow(x)){ min<-x[i,2]-500000; max<-x[i,2]+500000; subset<-x[x[,2]>min & x[,2]<max,]; te<-max(subset$GM)-min(subset$GM); if( (x[i,]$pos == subset[which.max(subset$GM),]$pos) & (max(subset$GM)> nu) & (te >=n) ){ E[i]<-1 }else{E[i]<-0;}}   
  y<-data.frame(x,E);
  max<-y[y$E==1,];

  write.table(max, file = "max", append = FALSE, quote = F, sep = "\t", row.names =F,col.names = TRUE);
 
  
