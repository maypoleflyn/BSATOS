#!/usr/bin/perl


 system("sort -k 1\,1 -k 2\,2n $ARGV[0] > tmp1;sort -k 1\,1 -k 2\,2n  $ARGV[1] > tmp2;filter -k cn tmp1 tmp2 > res;cut -f 1-4,7-8 res > res_end;rm tmp1 tmp2 res"); 
