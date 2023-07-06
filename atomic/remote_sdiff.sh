#!/bin/tcsh -f
set hope="cat /Users/jdh/atomic/"$1
ssh jdh@ibis.phyast.pitt.edu "$hope" | sdiff -s  - /nfs-share/home/hillier/atomic/$1 
