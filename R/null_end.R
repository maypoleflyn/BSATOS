#!/usr/bin/Rscript 

##Pvalue function for

 null<-function(x){
              ln<-log(x[8]);
              m1<-median(ln[,1]);
             colnames(ln)<-"ln";
             all<-data.frame(x,ln); 
              y<-ln[ln[,1]<=m1,];
              z<-m1-y
              m2<-median(z);
              values<-5.2*m2+m1;
             trimmed<-all[all[,9]<=values,] 
             return(trimmed);
          }
        trimmed<-null(data);
         m3<-as.data.frame(apply(trimmed[8],2,median))[1,1];
         u<-log(m3)
         u1<-log(u)
        qplot(trimmed[8],data=trimmed)

         M<-


         q<-(log(u-log(M)))/2;
        
         PV<-function(x){plnorm(x, meanlog = u1, sdlog = q, lower.tail = F, log.p = FALSE)}
          
         P<-apply(data[8],2,PV);
         colnames(P)<-"Pvalues";

        
           
        ajust<-function(x){
                 p.adjust(x, method = "fdr")
                 }
        p.ad<-apply(P,2,ajust);
         
        colnames(p.ad)<-"p.ajust";
         
        data_all<-data.frame(data,P,p.ad);  
        
        data_FDR<-data_all[data_all[,10]<=0.01,];
           
       data_FDR[data_FDR[10]==max(data_FDR[10]),][8]
        
       
