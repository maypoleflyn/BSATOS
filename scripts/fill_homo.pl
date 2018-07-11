#!/usr/bin/perl
 
   open V,"<$ARGV[1]";
   open P,"<$ARGV[0]";
    my %hash1;
    my %hash2;
    my @a;
    my $name1;
    my $name2;
   while(<V>){
        chomp;
     @a=split/\t/,$_;  
    unless($_=~/#/){
     if($_=~/DP=(.*?);/){
           if($1 >20 && $a[5] >30 ){
            if($a[9]=~/0\/0/){
                   $name1="0";
           $hash1{$a[0]}{$a[1]}=$name1;
                            }
             if($a[9]=~/1\/1/){   
                   $name1="1";    
            $hash1{$a[0]}{$a[1]}=$name1;
                            }
             if($a[10]=~/0\/0/){   
                   $name2="0";    
           $hash2{$a[0]}{$a[1]}=$name2;
                            }
              if($a[10]=~/1\/1/){   
                   $name1="1";  
             $hash2{$a[0]}{$a[1]}=$name2;  
                             }
                 }
            }
            }
           }                     
   
       my @b;
      while(<P>){
             chomp; 
         @b=split/\t/,$_;
        if(defined($hash1{$b[0]}{$b[1]}) && $b[4] eq "-" ){
                 $b[4]=$hash1{$b[0]}{$b[1]};
                 $b[5]=$hash1{$b[0]}{$b[1]};
                      
                                         }
        if(defined($hash2{$b[0]}{$b[1]}) && $b[7] eq "-" ){
                 $b[7]=$hash2{$b[0]}{$b[1]};
                 $b[8]=$hash2{$b[0]}{$b[1]}; 
                                         }
        my $c=join("\t",@b);
          print "$c\n";
               }  
     
            
        
                    
 
      
      
        
