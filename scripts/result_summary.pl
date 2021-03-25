#!/usr/bin/perl

=head1 NAME
 
 perl result_summary.pl  ------ the script to summary all the data in the end

=head1 Usage

 perl result_summary.pl <result dir>

       
=head1 Example

 perl result_summary.pl result_dir

=cut


      use warnings;
      use strict;
      use Getopt::Long;
      use Cwd;
      my ($h,$d);
        GetOptions(
           "d=i" =>\$d,
           "help"=>\$h);
 die `pod2text $0` if ($h);


 






