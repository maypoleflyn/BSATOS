#!/usr/bin/perl


#Copyright <2019> <CHINA AGRICULTURAL UNIVERSITY>
#
##Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
##The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
##THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
#












=head1 NAME

   gs.pl--------the script to conduct gs analysis

=head1 Usage

   gs.pl  -r <1000> -m <t/m> -p <precentage>  <genotype file> <phenotype file>
 

   -r: rounds(default:1000)
   -m: T test or just Mean method (t)
   -p: precentage (0.6)

=head1 Example

=cut

   use warnings;
   use strict;
   use Statistics::Basic qw(:all);
   use Statistics::PointEstimation; 
   use Getopt::Long;
   use Cwd;

   my $r;
   my $h;
   my $tm;
   my $p;

  GetOptions(
    "help"=>\$h,
    "m=s"=>\$tm,
    "p=f"=>\$p,
    "r=i" =>\$r);

   if($h || @ARGV==0){
            die `pod2text $0`;
             }
 
   unless(defined($r)){
 
          $r=1000;
                      }
   unless(defined($tm)){

         $tm="t";
                        }

   unless(defined($p)){            
            $p=0.6;
                      }
#####################read the data



###############################################################

my $time=localtime();
 
print STDERR "

\#\#\#\#\#


gs.pl--------the script to conduct gs analysis\n

...............................................$time\n";             
                                                   

$time=localtime();


print STDERR "

The fist stage reading data and generate trainning set and validation set\n 
...............................................$time\n";


$time=localtime();


       open LOK,"<$ARGV[1]";

        my %ph;
        my @sep;
        my @all;
        my @all_s;
 
     while(<LOK>){
             chomp;
             @sep=split/\t/,$_;
             $ph{$sep[0]}=$sep[1];
             push@all,$sep[1];
             push@all_s,$sep[0];
                 }

      close(LOK);

##################################################################
##################################################################

         my $max=$#all_s+1;    
         my $tr=int($p*$max);
         my $valid=$max-$tr;
         my @train_set;
         my @valid_set;  
         my $mb=0;
         my %gb;
         my %tr_c;
   
     for(my $ni=0;$ni<$tr;$ni++){  
    
        my $num=int(rand($#all_s));
        $gb{$num}++;
        if($gb{$num}==1){
         push@train_set,$all_s[$num];
             $tr_c{$all_s[$num]}++;  
                        }else{
                  $ni=$ni-1;   
                             }
                            }


print STDERR "

TAINNING SET SAMPLES\t$#train_set\n";


     foreach my $v_i(@all_s){   

            unless(defined($gb{$v_i})){
           
                push@valid_set,$v_i;  
                                      }
                          }  

print STDERR "

VALIDATION SET SAMPLES\t$#valid_set\n"; 
      
################################################################
#############################################################


  open LOP,"<$ARGV[0]";
 
   my %hash;
   my %h;
   my $c=0;
   my @name;
   my @t;
   my $tp_name;
   my %core;
   my @ak;
  while(<LOP>){  
         chomp;
         $c++;
      if($c==1){    
           @name=split/\t/,$_;          
               }else{
          @t=split/\t/,$_;
              my $cl=0;  
           foreach my $r(@t){      
                     
                      if($cl>=1){ 
                          $hash{$t[0]}{$name[$cl]}=$r;
                          $core{$name[$cl]}{$t[0]}=$r;
                          push@ak,$t[0];   
                          # print "$r\n";
                           unless($r eq "NA"){

                            if(defined($tr_c{$name[$cl]})){ 

                          push@{$h{$t[0]}{$r}},$ph{$name[$cl]};
                                                       
                                                      }
                                              }
                                }
                           $cl++;
                            }
             
                    }

              }
              
            
 close(LOP);

###############################################################

### calulate marker effect  
################################################################
  

 open MF,">marker_effect.txt";
 
 print MF "Marker\tAllele\tNum\tEffect\n";

    my @alle=sort keys %hash;
    my $me;
    my $va;
 
   if($tm eq "t"){
     $me=&tmean(@all); 
                }else{
     $me=&mean(@all);
                     }
    my %eff;
   # print "@alle\n"; 
   # print "$me\n";
   foreach my $index(@alle){   
                  
            my @k1=sort keys%{$h{$index}};   
     #       print "@k1\n";         
       foreach my $i2(@k1){  

                 my @res=@{$h{$index}{$i2}};
                 my $nb=$#res+1;
           
                 if($tm eq "t"){
                    $va=&tmean(@res)-$me; 
                               }else{
                    $va=&mean(@res)-$me;                    
                                    }
                  $eff{$index}{$i2}=$va;  
                   print MF  "$index\t$i2\t$nb\t$va\n"; 
                          }

                           }

############################################################
 
     # my @sample=sort keys%core;
      my @sample=@valid_set;
      my %gpv;
      my %epv;
      my $et;
      my $geno;
      my @res_a;
      my @res_b;

      open LOM,">EPV_file";

      foreach my $sindex(@sample){ 

               my @mk=sort keys %{$core{$sindex}};
  
         foreach my $ins(@mk){   
                 $geno=$hash{$ins}{$sindex};
                 if(defined($eff{$ins}{$geno})){   
                 $et=$eff{$ins}{$geno};        
                                               }else{
                 $et=0;                                 
                                                    }
                 $gpv{$sindex}+=$et;           
                         
                                
                             }

                 $epv{$sindex}=$me+$gpv{$sindex};
                  
                print LOM "$sindex\t$epv{$sindex}\t$ph{$sindex}\n";  
                
                
                                 }


####################################################################   

###################calculate the correlation
              
              my @pear=split/\n/,`Rscript ./pearson.R EPV_file`;
              my @p_cor=split/\s+/,$pear[0];
              my @p_pv=split/\s+/,$pear[1];
               print "$p_cor[1]\t$p_pv[1]\n";  

#############################################################

 sub man {

    my @m=@_;
    my $a=0;
    my $n=$#m+1;
    foreach my $i(@m){
               $a+=$i;
                      }
      my $av=$a/($n); 

      return($av);

          } 
#############################################

sub tmean {

  my @fm=@_; 
  my $stat = new Statistics::PointEstimation;
  $stat->set_significance(99); #set the significance(confidence) level to 95%
  $stat->add_data(@fm);
  
  return($stat->mean()); 
         
        }  

##############################################
#############################################

sub cmp {  
   
    (my $ar1,my $ar2)=@_;
     my @r1=@{$ar1};
     my @r2=@{$ar2};
    my $ttest = new Statistics::TTest;  
       $ttest->set_significance(90);
       $ttest->load_data(\@r1,\@r2);  
    my $sg=$ttest->significance;
     return($sg);
        }
   
###################################################
###################################################


sub pearson {

    (my $pr1,my $pr2)=@_;
    my $v1  = vector(@$pr1);
    my $v2  = vector(@$pr2);
    my $cor = corr($v1,$v2);  
    return($cor);
    
            } 


