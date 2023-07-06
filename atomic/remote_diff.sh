#!/bin/tcsh -f
set hope="cat /Users/jdh/atomic/"$1
ssh jdh@ibis.phyast.pitt.edu "$hope" | diff -sqy - /nfs-share/home/hillier/atomic/$1 
