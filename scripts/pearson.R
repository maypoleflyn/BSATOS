#!/usr/bin/Rscript

#####the script to conduct pearson analysis 

args <- commandArgs(T);
x<-read.table(args[1]);
c<-cor(x[2],x[3]);
print(as.matrix(c)[1,1]);
x<-apply(x,2,as.numeric) 
pv<-cor.test(x[,2],x[,3])[3]$p.value;
print(pv);





           
   
