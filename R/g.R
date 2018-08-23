#!/usr/bin/Rscript

G.test<-function (x, y = NULL, conservative = FALSE) 
    {
    DNAME <- deparse(substitute(x))
    if (is.data.frame(x)) 
        x <- as.matrix(x)
    if (is.matrix(x)) {
        if (min(dim(x)) == 1) 
            x <- as.vector(x)
    }
    if (!is.matrix(x) && !is.null(y)) {
        if (length(x) != length(y)) 
            stop("x and y must have the same length")
        DNAME <- paste(DNAME, "and", deparse(substitute(y)))
        OK <- complete.cases(x, y)
        x <- as.factor(x[OK])
        y <- as.factor(y[OK])
        if ((nlevels(x) < 2) || (nlevels(y) < 2)) 
            stop("x and y must have at least 2 levels")
        x <- table(x, y)
    }
    if (any(x < 0) || any(is.na(x))) 
        stop("all entries of x must be nonnegative and finite")
    if ((n <- sum(x)) == 0) 
        stop("at least one entry of x must be positive")
    if (!is.matrix(x)) 
        stop("Could not make a 2-dimensional matrix")
    nrows <- nrow(x)
    ncols <- ncol(x)
    sr <- apply(x, 1, sum)
    sc <- apply(x, 2, sum)
    E <- outer(sr, sc, "*")/n
    g <- 0
    for (i in 1:nrows) {
        for (j in 1:ncols) {
            if (x[i, j] != 0) 
                g <- g + x[i, j] * log(x[i, j]/E[i, j])
        }
    }
    q <- 1
    if (conservative) {
        row.tot <- col.tot <- 0
        for (i in 1:nrows) {
            row.tot <- row.tot + 1/(sum(x[i, ]))
        }
        for (j in 1:ncols) {
            col.tot <- col.tot + 1/(sum(x[, j]))
        }
        q <- 1 + ((n * row.tot - 1) * (n * col.tot - 1))/(6 * 
            n * (ncols - 1) * (nrows - 1))
    }
    STATISTIC <- G <- 2 * g/q
    PARAMETER <- (nrow(x) - 1) * (ncol(x) - 1)
    PVAL <- 1 - pchisq(STATISTIC, df = PARAMETER)
    if (!conservative) 
         METHOD <- "Log likelihood ratio (G-test) test of independence without correction"
    else METHOD <- "Log likelihood ratio (G-test) test of independence with Williams' correction"
    names(STATISTIC) <- "Log likelihood ratio statistic (G)"
    names(PARAMETER) <- "X-squared df"
    names(PVAL) <- "p.value"
    structure(list(statistic = STATISTIC, parameter = PARAMETER, 
        p.value = PVAL, method = METHOD, data.name = DNAME, observed = x, 
        expected = E), class = "htest")
}

  argv<-commandArgs(TRUE);
 
  x<-read.table(argv[1]);

 #s<-data.frame(x[1],x[2],x[5],x[6],x[7],x[8]);

 y<-x[x[1] == argv[2],];
 rm(x);
 colnames(y)<-c("chr","pos","R01_REF","R01_ALT","R02_REF","R02_ALT")
 G<-0;  

 data<-y;
 data$R01_REF_F<-(data$R01_REF)/(data$R01_REF + data$R01_ALT);
 data$R01_ALT_F<-(data$R01_ALT)/(data$R01_REF + data$R01_ALT);
 data$R02_REF_F<-(data$R02_REF)/(data$R02_REF + data$R02_ALT);
 data$R02_ALT_F<-(data$R02_ALT)/(data$R02_REF + data$R02_ALT);

  for(i in 1:nrow(data)){
      m<-matrix(c(data[i,3],data[i,5],data[i,4],data[i,6]),nc=2);
       t<-G.test(m)
       G[i]<-as.data.frame(t$statistic)[1,1]
        }

   
  data<-data.frame(data,G);
    
#####

ke<-function(x,w){
          H<-0;
   for(i in 1:nrow(x)){
           min<-x[i,2]-w/2;
           max<-x[i,2]+w/2;
          subset<-x[x[,2]>min & x[,2]<max,];
          D<-0;
          K<-0;
          A<-0;           
        
           D=abs(subset[2]-x[i,2])/w;
         
           A<-sum((1-D^3)^3);  
           K<-(1-(D)^3)^3/A;  
                                      
           n<-subset[11]*K;
           H[i]<-sum(n);
         
     }
                 return(H);
    }      


  
  size<-as.numeric(argv[4]);
  size1<-0.75 * size;
  size2<-0.5 * size;
  size3<-0.25*size;

 if(argv[5] == 1){    
  G1M<-ke(data,size); 
  G750<-ke(data,size1);
  G500<-ke(data,size2);
  G250<-ke(data,size3);            
  data<-data.frame(data,G1M,G750,G500,G250);
                  }else{
    G1M<-ke(data,size);
  data<-data.frame(data,G1M);     
                       }
                       
        name=paste(argv[2],argv[3],"afd",sep="_");
        out<-as.character(name);
        write.table(data,file=name,quote =F,sep = "\t",row.names = F,col.names = F) 


