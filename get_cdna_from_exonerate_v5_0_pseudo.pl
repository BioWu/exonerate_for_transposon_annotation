#use warnings;
#use strict;
use String::Index qw( cindex ncindex crindex ncrindex );
use List::MoreUtils ':all';
use Bio::Seq;
use Bio::SeqIO;

#my $file = "hg19.Step_2.get_cdna_v2.0.out.test";
my $file = "$ARGV[0]";
my %te_num;
my %te_num_pseudo;

sub from_bed_to_fa(@bed_record){
	$tem = shift(@_);
	@bed_record = split("\t",$tem);
	
	$pseudo_flag = $bed_record[3];
	
   $id = $bed_record[4];
	
	$seq = $bed_record[6];
	
	$pos = $bed_record[0]."_".$bed_record[1]."_".$bed_record[2]."_".$bed_record[5];
	
		#### for sequence output
	if(! (defined($te_num{$id})||defined($te_num_pseudo{$id}))){
#		print $id;
			$fasta_file_name = $id.".actived.fa_test";
			$fasta_file_name_2 = $id.".pseudo.fa_test";
			$te_num{$id} =0 ;
			open(NT,">>",$fasta_file_name);
			close(NT);
			open(NT,">>",$fasta_file_name_2);
			close(NT);
		   }
	   $fasta_file_name = $id.".actived.fa_test";
		$fasta_file_name_2 = $id.".pseudo.fa_test";
	my $seqout_file_obj = Bio::SeqIO->new(-file => ">>$fasta_file_name",-format => 'Fasta');
	my $seqout_file_obj_2 = Bio::SeqIO->new(-file => ">>$fasta_file_name_2",-format => 'Fasta');
	
		if($pseudo_flag){
		$te_num_pseudo{$id} += 1;
		my $seq_obj = Bio::Seq->new(
									-seq => $seq,
									-id  => $id."__".$pos."__".$te_num_pseudo{$id},
									-alphabet => 'dna');
		$seqout_file_obj_2->write_seq($seq_obj);
	}else{
		$te_num{$id} += 1;
		my $seq_obj = Bio::Seq->new(
									-seq => $seq,
									-id  => $id."__".$pos."__".$te_num{$id},
									-alphabet => 'dna');
		$seqout_file_obj->write_seq($seq_obj);
	}
}
sub parse_bed($file_bed){
	open BED , "<" , shift(@_);
	while(<BED>){
		my @beds = split();
		if(/\&/){
			chomp;
			# has overlap
			my @chrs = split("&",$beds[3]);
			my @starts = split("&",$beds[4]);
			my @ends = split("&",$beds[5]);
			my @pseud_flags = split("&",$beds[6]);
			my @TE_ids = split("&",$beds[7]);
			my $strand = $beds[8];
			my @seqs = split("&",$beds[9]);
			my @lengths = split("&",$beds[10]);
			## Activty first then length
#			my @elements_will_be_removed_t = indexes {$_ =~ /\#|-|\+|\{|\}/} @target_ele;
			my @active_TE = indexes {$_ =~ /0/} @pseud_flags;
			if($#active_TE){
				# has active TEs
				my($length_min,$length_max) = minmax(@lengths[@active_TE]);
				my(@max_length_index) = indexes {$_ =~ $length_max} @lengths[@active_TE];
				if($#max_length_index > 1){print STDERR "at least two same-length TE overlaped and active. \n "}else{
					$bed = join("\t",$chrs[$max_length_index[0]],$starts[$max_length_index[0]],$ends[$max_length_index[0]],
											$pseud_flags[$max_length_index[0]],$TE_ids[$max_length_index[0]],$strand,$seqs[$max_length_index[0]]
											),"\n";
					# $bed ::: chrY	25291721	25294365	1	prirep_L1-1_Cja_32538_1	+	agggggaa
					&from_bed_to_fa($bed);
				}
			}else{
				# has no active TE
				my($length_min,$length_max) = minmax(@lengths);
				my(@max_length_index) = indexes {$_ =~ $length_max} @lengths;
				if($#max_length_index > 1){print STDERR "at least two same-length TE overlaped and active. \n "}else{
					$bed = join("\t",$chrs[$max_length_index[0]],$starts[$max_length_index[0]],$ends[$max_length_index[0]],
											$pseud_flags[$max_length_index[0]],$TE_ids[$max_length_index[0]],$strand,$seqs[$max_length_index[0]]
											),"\n";
					# $bed ::: chrY	25291721	25294365	1	prirep_L1-1_Cja_32538_1	+	agggggaa
					&from_bed_to_fa($bed);
				}
			}
			 
		}else{
			# no_overlap
			chomp;
			$bed = join("\t",@beds[3..9]);
			# $bed ::: chrY	25291721	25294365	1	prirep_L1-1_Cja_32538_1	+	agggggaa
			&from_bed_to_fa($bed);
		}
	}
}
sub get_bed_file($chr,$start,$end,$pseudo_flag,$TE_id,$strand,$sequence){
	## get one seq_obj then return one bed-record
	my @bed = @_;
	$bed[7] = length($bed[6]);
	return join("\t",@bed);
}
sub get_one_record($file){
	local $/ = " >\n\n";
	my $file_handle = shift(@_);
	open(FH,$file_handle) or die 'can not open $file\n';
	my @records;
	while(<FH>){
		$tmp = $_;
		$tmp =~ s/^\n//;
		push @records,$tmp;
	}
#	our @records = <FH>;
	return @records;
}
sub parse_exonerate_per_record($singlerecord){
	my @sequence;
	my @record_lines = split('\n',shift(@_));
	#Isertion => delete
	#Frameshift => delete
	#stopcodon => changes to NNN
	#with N => defined to pseudo-transposon
	#if it is an psudo gene please note it on id
	## Query: humrep_L1PREC2;Transposable_Element;FT;L1PREC2_2p;14809;2
	$record_lines[2] =~ /:\s(.*?);.*?;(\w*);(\w*$)/;
	my $id=join("_",$1,$2,$3);
	#0 for actived ; 1 for frameshift ; 2 for stop_codon ; 3 for N ; 4 for frameshift&stopcodon  
	my $pseudo_flag = 0;
	## offset Target range: 1 -> 3635
	if($record_lines[7] =~ /:\s(\d*)\s\S+\s(\d*)/){
		$tmp_offset[0] = $1;
		$tmp_offset[1] = $2;
		#if complement to target sequence, then throw this sequence . 
		if($tmp_offset[0] > $tmp_offset[1]){
			print $record_lines[7],"\n";
			return 0;
		}
		($offset_1,$offset_2)=minmax(@tmp_offset);
		print $record_lines[7],"\n",join("\t",$tmp_offset[0],$tmp_offset[1]),"\n";
	}else{
		die "Can not get target-ranges.\n"; 
	}
	## Target: chr1:73900-77690(+)
	$record_lines[3] =~ /:\s(.*?):(\w*)-(\w*)\(([+|-])\)/;
	my ($chr,$start,$end,$strand) = ($1,$2+$offset_1,$2+$offset_2,$4);
	print join("\t",$chr,$start,$end,$strand),"\n";
	my $pos = join('_',$chr,$start,$end,$strand);
	# from(12,13) with stepsize 5；while( ! ($record_lines[12+5i] =~ /^#/ && $record_lines[13+5i] =~ /^#/) )
	my $i = 0;
	while( ! ($record_lines[11+5*$i] =~ /^#/ && $record_lines[12+5*$i] =~ /^#/) ){
		
#		my $original_target = $record_lines[10+5*$i];

		my $target = $record_lines[11+5*$i];
		$target =~ /^\s*(.*?)\s*$/;
		$target = $1;
		my $query = $record_lines[12+5*$i];
		$query =~ /:\s*(.*?)\s*:/;
		$query = $1;
		$i++;
		$add = abs(length($target) - length($query));
		#for frameshift or stop codon or intron
		if(length($target) < length($query) && ($query =~ /^[a|t|c|g|n|N]/)){
			#"-" start from left
			$target = $target.('-'x$add);
		}else{
			# '-'start from right
			$target = ('-'x$add).$target;
		}
		if(length($target)==length($query)&& ($target =~ /\*|\#|-|\+/i || $query =~ /n|-|\.|\{|\}|N/i)){
			$pseudo_flag = 1;
			my @target_ele = split('',$target);
			my @query_ele = split('',$query);
			#get removal elements indexs
			#for frame_shift || N 
			# Notice:: for N sites , we choose to keep
			my @elements_will_be_removed_t = indexes {$_ =~ /\#|-|\+|\{|\}/} @target_ele;
			my @elements_will_be_changed_t = indexes {$_ =~ /\*/i} @target_ele;
			# I choose to keep N sites(unknown sites).
			my @elements_will_be_removed_q = indexes {$_ =~ /\.|-|\{|\}/i} @query_ele;
			push my @elements_will_be_removed_mix ,@elements_will_be_removed_t,@elements_will_be_removed_q;
			@elements_will_be_removed_mix_uniq = uniq @elements_will_be_removed_mix;
#			print join("\t",@elements_will_be_removed_mix_uniq,@elements_will_be_changed_t),"\n";
			foreach(@elements_will_be_changed_t){
				$query_ele[$_]='N';
			}
			foreach(@elements_will_be_removed_mix_uniq){
					$query_ele[$_]='';
			}
#			print "Raw:";
#			print $query,"\nNow:";
#			print join('',@query_ele),"\n";
			$query = join('',@query_ele);
		}
		push @sequence,$query;
		
	}
	my $seq = join('',@sequence),"\n";
	##### for bed output
		return &get_bed_file($chr,$start,$end,$pseudo_flag,$id,$strand,$seq);

}
##
	
## main ##
my $file_b4sort = $file.".b4sort";
my $file_b4merge = $file.".b4merge";
my $file_bed = $file.".bed";

my @records = &get_one_record($file);

`rm $file_b4sort`;

open B4SORT ,">>",$file_b4sort;
print STDERR "getting bed file .\n";
foreach(@records){
	# for debug only
	if(&parse_exonerate_per_record($_)){
	print B4SORT &parse_exonerate_per_record($_),"\n";
	}else{
		print STDERR "\nComplement to target sequence !\n";
	}
}
`bedtools sort -i $file_b4sort >$file_b4merge`;
`bedtools merge -i $file_b4merge  -s -c 1,2,3,4,5,6,7,8 -o collapse,collapse,collapse,collapse,collapse,distinct,collapse,collapse -delim "&" >$file_bed`;
`rm *.fa_test`;
print STDERR 'pasrsing bed .\n';
&parse_bed($file_bed);

