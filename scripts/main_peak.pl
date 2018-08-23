#!/usr/bin/perl

=head1 NAME

   main_peak.pl ------ the script to get the main peak from the polished or profile data

=head1 Usage

  perl main_peak.pl --th1 <threshold value of G1 data> --th2 <threshold value of G2 data> --th3 <the threshold value of G3 data> --g1 <the G1 data from polish step> --g2 <the G2 data from polish step> --g3 <the G3 data from polish step> --outprefix <out prefix>

  --th1: the threshold value of G1 data
  --th2: the threshold value of G2 data
  --th3: the threshold value of G3 data 
  --g1: the G1 data from polish step
  --g2: the G2 data from polish step 
  --g3: the G3 data from polish step
  --outprefix <out prefix>

=head1 Example

=cut

   use warnings;
   use strict;
   use Getopt::Long;
   use Cwd; 
   use File::Basename;
   use FindBin qw($Bin); 
   my ($th1,$th2,$th3,$g1,$g2,$g3,$o,$h);
 
   GetOptions(
           "g1=s"=>\$g1,
           "g2=s"=>\$g2,
           "g3=s"=>\$g3,
           "th1=s"=>\$th1,
           "th2=s"=>\$th2,
           "th3=s"=>\$th3,
           "outprefix=s"=>\$o, 
           "help"=>\$h);

 die `pod2text $0` if ($h);


 my $gn1="P";
 my $gn2="M";
 my $gn3="H";
 
 &get_peak($g1,$th1,$gn1);
 &get_peak($g2,$th2,$gn2);
 &get_peak($g3,$th3,$gn3); 



 sub get_peak {
  
   my @par=@_;      
   open LOD,"<$par[0]";
   my @a=();
   my $th=$par[1]; 
   my @st=();
   my @st2=();
   my @b=();
   my $chr;
  while(<LOD>){
          chomp;
      @a=split/\t/,$_;
      $chr=$a[0];
      push@st,$a[12];
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
           
    my $count=0;
    my @tmp;
    my $s_q;
    my $e_q; 
    my $acc;
    
   foreach my $index(@b){
             $count++;
              push@tmp,$index;             
         if($count==2){
                      
              my $get_max="$Bin/../R/get_max.R";       
              my $value=`Rscript  $get_max  $par[0] $tmp[0] $tmp[1]`;
              my @w=split/\n/,$value;
              my @j=split/\s+/,$w[0];
              my @k=split/\s+/,$w[1];                        
              $count=0;
             $s_q=$k[1]-500000;
             $e_q=$k[1]+500000;
            if($s_q <0){
                  $s_q=0;
                      }
           if($s_q < $tmp[0]){
                    $s_q=$tmp[0];
                             }
           if($e_q > $tmp[1]){
                    $e_q=$tmp[1];
           
                  } 

           if(($tmp[1]-$tmp[0]) >250000){
                  $acc++;
                  my $name=$par[2].$acc;
           print "$par[2]\t$chr\t$s_q\t$e_q\t$k[1]\t$j[1]\t$name\n";
                                        }


               @tmp=();
                 }
           }
           }
    
 


  
 

