#!/usr/bin/perl
use List::MoreUtils qw{
        uniq distinct minmax 
    };

##set bin number = 100
$bin_num = 100;
open FH , "$ARGV[0]";

while(<FH>){
	my @inline = split("\t");
	my($TE,$num,$pos,$freq,$length) = $inline[0],$inline[1],$inline[2],$inline[3];
	my $index = $TE.$num;
	$TE{$index}{$pos} = $freq;
	$TE{$index}{'length'} = $length;
}

foreach(sort keys(%TE)){
	$id_tmp = $_;
	my @poss = keys($TE{$id_tmp});
	$min_pos,$max_pos = minmax(@poss);
	##
	for($i=0; $i <= $length / $bin_num ; $i++){
		
	} 
}