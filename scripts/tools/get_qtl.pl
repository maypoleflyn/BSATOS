#!/usr/bin/perl


  use warnings;
  use strict;
  
  my $script1=$ARGV[0]."_"."1R";
  my $outfile=$ARGV[0]."_Gvale";
  open LOL,">$script1";
  print LOL "load(file=\"$ARGV[0]\")\nwrite.table(data,file=\"$outfile\",append = FALSE, quote =F, sep =\"\\t\",row.names = F,col.names = F)\n";
  
  system("Rscript $script1 $ARGV[0]");
  


  open LOD,"<$outfile";
   my @a;
   my $th=$ARGV[1]; 
   my @st=();
   my @st2=();
   my @b;
   my $chr;
  while(<LOD>){
          chomp;
      @a=split/\t/,$_;
      $chr=$a[0];
       push@st,$a[7];
      push@st2,$a[1];
    if($#st==0){
           next;
       }
    if($#st==1){
        if($st[0] >$th && $st[1] < $th){
                push@b,$st2[0];
             }
        if($st[0] <$th && $st[1] > $th){
                push@b,$st2[1];
             }
        shift@st;
        shift@st2;
      } 
    }
           
##############################
 
#    print "@b\n";
    my $count=0;
    my @tmp;
   foreach my $index(@b){
             $count++;
              push@tmp,$index;             
         if($count==2){
              my $value=`Rscript get_max.R $ARGV[0] $tmp[0] $tmp[1]`;
              my @w=split/\n/,$value;
              my @j=split/\s+/,$w[0];
              my @k=split/\s+/,$w[1];                        

         
              $count=0;
             
         print "$chr\t$tmp[0]\t$tmp[1]\t$k[1]\t$j[1]\n";
               @tmp=();
                 }
           }
    
  system("rm $script1 $outfile");   
 


  
 
