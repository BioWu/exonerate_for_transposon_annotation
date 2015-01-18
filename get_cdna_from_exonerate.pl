#!usr/bin/perl
#file = humrep_HERV4_I_Transposable_Element_FT_HERV4_I_2p_14357_2.calJac

#use Bio::SearchIO;
#my $searchio = Bio::SearchIO->new(-file => 'humrep_HERV4_I_Transposable_Element_FT_HERV4_I_2p_14357_2.calJac',
#                                 -format => 'exonerate');
# 
# 
#while( my $r = $searchio->next_result ) {
#  print $r->query_name, "\n";
#  print $i++;
#  
#}

open(FH,"<","$ARGV[0]");

#open(FH,"<","humrep_HERV4_I_Transposable_Element_FT_HERV4_I_2p_14357_2.calJac");
my $frame_shift = 0;
my $cdna_turn = -1;
while(<FH>){
	if(/exonerate:protein2genome:bestfit/&&/\t/){	#gff2 turn and frameshifts turn
		
		if(/frameshifts/){
			$frame_shift = 1;
		}
		if(/\tgene\t/){
#			print $i++,"\n";
		@gff = split("\t");
		$gff[0] =~ /(\w*?):(\d*?)-(\d*?)\(/;
		$gff[0] = $1;
		$gff[1] = $range[1]*3/($gff[4]-$gff[3]);
		$gff[3] +=  $2;
		$gff[4] +=  $2;
		}
		
		
	}elsif(/^>>/){#cdna-start end TURNING
		$cdna_turn *= -1;
	}elsif(/Query range:/){
		$tmp=$_;
		chomp($tmp);
		@range = split(">",$tmp);
	}elsif(/^ >/){#start output
		if($frame_shift == 1){
			undef @gff;
			undef @cdna;
			undef @range;
			$cdna_turn = -1;
		   $frame_shift = -1;
		}else{
		#print gff
		print join("\t",@gff[0..6]);
#		$cdna[0] =~ /(\w*?);/;
	@id = split(";",$cdna[0]);
	$id[0] =~ s/>>//;
		print "\t",$id[0]."_".$id[3],"\t",@cdna,"\n";
		undef  @gff;
		undef @cdna;
		undef @range;
		$cdna_turn = -1;
		$frame_shift = -1;
		}
	}
	if($cdna_turn == 1 ){
		$dna_line = $_;
		chomp($dna_line);
		push(@cdna,$dna_line);
	}
	
}
close(FH);