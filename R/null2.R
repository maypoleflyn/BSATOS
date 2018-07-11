#!/usr/bin/Rscript 

##Pvalue function for

 null<-function(x){
              x<-log(data[8]);
              m1<-median(x[,1]);
             colnames(x)<-"x";
             all<-data.frame(data,x); 
              y<-x[x[,1]<=m1,];
              m2<-m1-median(y);
              values<-5.2*m2+m1;
            trimmed<-all[all[,9]<=values,] 
             return(trimmed);
          }
        trimmed<-null(data);
         m3<-as.data.frame(apply(trimmed[8],2,median))[1,1];
         u<-log(m3)
        
        qplot(trimmed[8],data=trimmed)

        u1<-log(Mode);
        r<-sqrt(u-u1);
       
        
         PV<-function(x){pnorm(x, mean = u, sd=r , lower.tail = F)}
          
         P<-apply(data[9],2,PV);
         colnames(P)<-"Pvalues";

        
           
        ajust<-function(x){
                 p.adjust(x, method = "fdr")
                 }
        p.ad<-apply(P,2,ajust);
         
        colnames(p.ad)<-"p.ajust";
         
        data_all<-data.frame(data,P,p.ad);  
        
        data_FDR<-data_all[data_all[,11]<=0.01,];
           
       data_FDR[data_FDR[11]==max(data_FDR[11]),][8]
        
       
