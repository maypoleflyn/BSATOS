#!/usr/bin/Rscript 

 argv <- commandArgs(trailingOnly=T)
 cat("Importing files.\n")
 
 data <- read.table(argv[1],header=T);
 
 colnames(data)<-c("chr","pos","R01_REF","R01_ALT","R02_REF","R02_ALT")
 chr<-argv[2];

 data$R01_REF_F<-(data$R01_REF)/(data$R01_REF + data$R01_ALT);
 data$R01_ALT_F<-(data$R01_ALT)/(data$R01_REF + data$R01_ALT);
 data$R02_REF_F<-(data$R02_REF)/(data$R02_REF + data$R02_ALT);
 data$R02_ALT_F<-(data$R02_ALT)/(data$R02_REF + data$R02_ALT);

 data<-data[data[1]==chr,]   
 euc<-sqrt((data$R01_REF_F-data$R02_REF_F)^2 + (data$R01_ALT_F - data$R02_ALT_F)^2);
 data$euc4<-(euc)^4;
 

############################################################################################################## 
##############################################################################################################
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

            return(res[,3]);
                                    }           										  				     
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

#######################################################################################################
#######################################################################################################

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



#######################################################################################################
#######################################################################################################



  size<-as.numeric(argv[4]);
  
  fn<-as.numeric(argv[6]);
  size1<-(0.75)*size;
  size2<-(0.5)*size;
  size3<-(0.25)*size;

 if(argv[5] == 1){    
  D1M<-ke(data,size,fn); 
  D750<-ke(data,size1,fn);
  D500<-ke(data,size2,fn);
  D250<-ke(data,size3,fn);            
  data<-data.frame(data,D1M,D750,D500,D250);
                  }else{
    D1M<-ke2(data,size);
    D750<-ke2(data,size1);
    D500<-ke2(data,size2);
    D250<-ke2(data,size3);   
  data<-data.frame(data,D1M,D750,D500,D250);   
                       }
                       
        name=paste(argv[2],argv[3],"afd",sep="_");
        out<-as.character(name);
        write.table(data,file=name,quote =F,sep = "\t",row.names = F,col.names = F) 

   







 


  

 
