#!/bin/bash
ls *Rdata |grep -- "Ch" |xargs -i echo "get_qtl.pl2 {} $1  " |sh 
