#!/usr/bin/perl 


=head1 NAME
 
 get_bed.pl  ------ the script to get genes in the QTL region

=head1 Usage

  perl get_bed.pl <gene bed file> <chromosome> <start> <end>  

        
 

=head1 Example

=cut

   use warnings;
   use strict;
 
   use Getopt::Long;
   use Cwd;
   my ($h,$header);
        GetOptions(    
	 "header=s" =>\$header,      
	 "help"=>\$h);
 die `pod2text $0` if ( @ARGV==0 || $h);
 
   
   unless(defined($header)){
	    
	     $header="F";
	          
	                   }


  my @a; 

  open LOD,"<$ARGV[0]";
 
 while(<LOD>){
     chomp;
     @a=split/\t/,$_;

   if($_=~/#/){

	if($header eq "F"){

		   next; 
	                  }else{
			   
           print "$_\n";
				  
			       }	  
	      }

   if($a[0] eq $ARGV[1] && ($a[1] >= $ARGV[2]) && ($a[1] <=$ARGV[3])){
           print "$_\n";
       }
       }  
  
