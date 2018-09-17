#!usr/bin/Rscript


getPvals <-function(Gprime)
    {
        
  
           
            lnGprime <- log(Gprime);
            
            medianLogGprime <- median(lnGprime);
            
            # calculate left median absolute deviation for the trimmed G' prime set
            MAD <-median(medianLogGprime - lnGprime[lnGprime <= medianLogGprime]);
            
            # Trim the G prime set to exclude outlier regions (i.e. QTL) using Hampel's rule
            trimGprime <-Gprime[lnGprime - median(lnGprime) <= 5.2 * MAD];
        
        
        medianTrimGprime <- median(trimGprime);
        
        # estimate the mode of the trimmed G' prime set using the half-sample method
      
        modeTrimGprime <-modeest::mlv(x = trimGprime, bw = 0.5, method = "hsm")$M
        
        muE <- log(medianTrimGprime);
        varE <- abs(muE - log(modeTrimGprime));
        #use the log normal distribution to get pvals
  
        pval <- 1 - plnorm(q = Gprime,
                meanlog = muE,
                sdlog = sqrt(varE));
        
        return(pval);
    }


 getFDRThreshold <- function(pvalues, alpha = 0.01)
{
    sortedPvals <- sort(pvalues, decreasing = FALSE);
    pAdj <- p.adjust(sortedPvals, method = "BH");
    if (!any(pAdj < alpha)) {
        fdrThreshold <- 0
    } else {
    fdrThreshold <- sortedPvals[max(which(pAdj < alpha))];
    }
    return(fdrThreshold);
}


 argv<-commandArgs(TRUE);
 
 x<-read.table(argv[1]);
 alpha<-as.numeric(argv[2]);
  
  z1<-getPvals(x[,12]);
  z2<-getPvals(x[,13]);
  z3<-getPvals(x[,14]);
  z4<-getPvals(x[,15]);

  y1<-getFDRThreshold(z1,alpha);
    
  if(y1 == 0){
  y1<-apply(x[12],2,median)+2*apply(x[12],2,sd);
                }
 
  y2<-getFDRThreshold(z2,alpha);

  if(y2 == 0){  
  y2<-apply(x[13],2,median)+2*apply(x[13],2,sd);
                }
    
  y3<-getFDRThreshold(z3,alpha);

  if(y3 == 0){
  y3<-apply(x[14],2,median)+2*apply(x[14],2,sd);
                }

  y4<-getFDRThreshold(z4,alpha); 

  if( y4 == 0){
  y4<-apply(x[15],2,median)+2*apply(x[15],2,sd);
                 }
 
  cat(y1,file=as.character(argv[3]),sep="\n",append = T);
  cat(y2,file=as.character(argv[3]),sep="\n",append = T);
  cat(y3,file=as.character(argv[3]),sep="\n",append = T);
  cat(y4,file=as.character(argv[3]),sep="\n",append = T);  









