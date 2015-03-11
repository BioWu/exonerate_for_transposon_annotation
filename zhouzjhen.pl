while(<STDIN>){
	
	@inline = split("\t");
	if($inline[5] eq $inline[-2]){
		$flag=1;
	}else{
		$flag=1;
	}
	$TE{$inline[-3]}{$inline[10]}{$inline[1]-$inline[10]}{'1'} = $inline[4] * $flag; 
	$TE{$inline[-3]}{$inline[10]}{$inline[1]-$inline[10]}{'2'} = $inline[11] - $inline[10]; 
}

foreach(sort keys(%TE)){
	
	$key_tmp=$_;
	$i=0;
	foreach(sort keys($TE{$key_tmp})){
		
		$i++;
		
		$key_tmp2=$_;
		
		foreach(sort keys($TE{$key_tmp}{$key_tmp2})){
			print "$key_tmp\t$i\t$_\t$TE{$key_tmp}{$key_tmp2}{$_}{'1'}\t$TE{$key_tmp}{$key_tmp2}{$_}{'2'}\n";
		}
	}
}