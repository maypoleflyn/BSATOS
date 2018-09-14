#!/usr/bin/perl


=head1 NAME
  filter_raw_sv.pl ------ the script to filter the raw SV VCF file

=head1 Usage

  perl filter_raw_sv.pl --d <depth> --mq <quality> -l <sv length> --precise yes

  --d: the theshold of coverage depth [10]
  --mq: the mapQ quality  [20]
  --l: the minimum SV length [1000]
  --precise: only kept precise SVs [yes]
 
=head1 Example

=cut

   use warnings;
   use strict;
   use Getopt::Long;
   use Cwd;
   my ($d,$mq,$h,$l,$pre);
  GetOptions(
           "d=i"=>\$d,
           "mq=i"=>\$mq,
           "l=i"=>\$l,
           "precise=s" =>\$pre,
           "help"=>\$h);

 die `pod2text $0` if (@ARGV==0 || $h);
  
     unless(defined($d)){
                     $d=10;
                        }
     unless(defined($mq)){
                     $mq=20;
                         }
 
     unless(defined($pre)){
                    $pre="yes";
                          }
     unless(defined($l)){
                    $l=1000;
                        }
       
    
    my @a; 
  while(<>){
        chomp;
     if($_=~/#/){
       print "$_\n"; 
               }else{
          @a=split/\t/,$_;      
         if($_=~/\.\/\./){
                    next;
                         }
         if(defined($pre)){               
             if($_ =~ /IMPRECISE/){
                 next;        
                                  }
                          }
          my $start=$a[1];
          $_=~/END=(\d+);/;
          my $end =$1;
          my $le=$end-$start;
          $_=~/MAPQ=(\d+);/;
          my $vq=$1;
          $_=~/PE=(\d+);/;
          my $pe=$1;
          $_=~/SR=(\d+)/;
          my $sr=$1;
          my $da=$sr+$pe;
         if($da >=$d && $vq >=$mq && $le <=$l){
            print "$_\n";
                                            }
              }
             }        
          
          



         

       
