#!/usr/bin/perl
#
# changes all local files named *_in.dat to *_IN
#
@outlist = split(" ",`ls *.dat`);
foreach $outfile (@outlist) {
    $infile = $outfile;
    $infile =~ s/\.dat//;             #chop off .DAT
    `mv $outfile $infile`;            #rename file on disk
}
