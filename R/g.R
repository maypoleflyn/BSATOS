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


####################################################################################################
###################################################################################################
ke<-function(x,w,fn){

#####fn bigger, slower,but low memory   

  if(fn==1){
           
           sm<-smooth(x,w);
           return(sm[,3]);
           
          }else{
           

            start<-x[1,2];
            end<-x[nrow(x),2];			
            n<-ceiling(end/w);
					
        if(fn>n){
		    fn<-n;      
			    }else{
			fn<-fn;	 		  
			    }
				
				
	     thes<-n/(sqrt((end/w)-(1/4)));
                    
          if(fn<thes){
		      fn<-2*ceiling((ceiling(fn)/2));
		    
		       }else{ 
			   fn<-fn;
				  }		 
								
          W=w*round(n/fn);
			
          splice<-seq(0,end,by=W);
		  
          array<-c(splice,end);    			
           
      for(i in c(1:(length(array)-1))){
	  		
	        k<-i+1;
			if(i==1){
			   s<-array[i];
			   stim<-array[i]; 
			       }else{
	           s<-array[i]-w/2;
			   stim<-array[i]; 
			            }
			if(i==length(array)-1){
			 
              e<-array[k];
			  etrim<-array[k]; 
			                      }else{
              e=array[k]+w/2;
              etrim<-array[k];			  
 								   
							     	   }
			
	        sub<-x[x[,2]>=s & x[,2] <=e,];
	        sm<-smooth(sub,w);
            sm_s<-sm[sm[,2] >=stim &sm[,2] <=etrim,];  	
			if(i==1){
			   res<-sm_s;
			        }else{ 
            res<-rbind(res,sm_s);
                         }			
                  	             }

            }
            return(res[,3]);
                                   										  				     
				        }


smooth<-function(x,w){
 
        pos<-x[,2];
		g<-x[,11];
        l<-length(pos);
        vec<-t(matrix(rep(1,l)));
        A<-pos%*%vec;
		B<-t(A);
        X<-A-B;
        X<-abs(X);
        
        f<-matrix(rep(1,l^2),nrow=l,ncol=l);
        f[X<=w/2]=1;              
	    f[X>w/2]=0; 
		D<-X*f;
		D<-D/w;
		D<-(1-D^3)^3;
		D<-D*f;
		D_sum<-apply(D,1,sum)
		E<-D_sum%*%vec;
		K<-D/E;
		
		G<-g%*%vec;
		H<-t(G);
		
		N<-apply(H*K,1,sum);
		dat<-data.frame(x[,1],x[,2],N);
		return(dat);
		 
		}

##########################################################################
##########################################################################

ke2<-function(x,w){
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









############################# 1 data; 2 chromosome; 3 snp_type 4 window size 5 mutiple or not 6 fn number; 

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

#########################################################################

  size<-as.numeric(argv[4]);
  fn<-as.numeric(argv[6]);
  size1<-0.75 * size;
  size2<-0.5 * size;
  size3<-0.25*size;

 if(argv[5] == 1){    
  G1M<-ke(data,size,fn); 
  G750<-ke(data,size1,fn);
  G500<-ke(data,size2,fn);
  G250<-ke(data,size3,fn);            
  data<-data.frame(data,G1M,G750,G500,G250);
                  }else{


    G1M<-ke2(data,size);
    G750<-ke2(data,size1);
    G500<-ke2(data,size2);
    G250<-ke2(data,size3);

  data<-data.frame(data,G1M,G750,G500,G250);



                       }
                       
        name=paste(argv[2],argv[3],"afd",sep="_");
        out<-as.character(name);
        write.table(data,file=name,quote =F,sep = "\t",row.names = F,col.names = F) 






