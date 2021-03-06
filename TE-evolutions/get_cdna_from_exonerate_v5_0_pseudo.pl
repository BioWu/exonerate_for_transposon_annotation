#use warnings;
#use strict;

use String::Index qw( cindex ncindex crindex ncrindex );
use List::MoreUtils ':all';
use Bio::Seq;
use Bio::SeqIO;

my $file = "calJac.Step_2.get_cdna_v2.0.out.test";
#my $file = "$ARGV[0]";
my %te_num;
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
	$record_lines[7] =~ /:\s(\d*?)\s->\s(\d*)/;
	my $offset_1 = $1;
	my $offset_2 = $2;
	## Target: chr1:73900-77690(+)
	$record_lines[3] =~ /:\s(.*?):(\w*)-(\w*)\(([+|-])\)/;
	if($4 eq '+'){
		my ($chr,$start,$end,$strand) = ($1,$2+$offset_1,$2+$offset_2,$4);
	}else{
		my ($chr,$start,$end,$strand) = ($1,$3-$offset_2,$3-$offset_1,$4);
	}
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
		if(length($target) < length($query)){
			$target = ('-'x$add).$target;
#			print $target;
		}
		if(length($target)==length($query)&& ($target =~ /\*|\#|-|\+/i || $query =~ /n|-|\.|\{|\}|N/i)){
			$pseudo_flag = 1;
			my @target_ele = split('',$target);
			my @query_ele = split('',$query);
			#get removal elements indexs
			#for frame_shift || N 
			my @elements_will_be_removed_t = indexes {$_ =~ /\#|-|\+|\{|\}/} @target_ele;
			my @elements_will_be_changed_t = indexes {$_ =~ /\*/} @target_ele;
			my @elements_will_be_removed_q = indexes {$_ =~ /n|\.|-|\{|\}/i} @query_ele;
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
	if(! defined($te_num{$id})){
#		print $id;
			$fasta_file_name = $id."actived.fa_test";
			$fasta_file_name_2 = $id."pseudo.fa_test";
			$te_num{$id} =0 ;
			open(NT,">>",$fasta_file_name);
			close(NT);
			open(NT,">>",$fasta_file_name_2);
			close(NT);
		   }
		$te_num{$id} += 1;
		$te_num_pseudo{$id} += 1;
	my $seqout_file_obj = Bio::SeqIO->new(-file => ">>$fasta_file_name",-format => 'Fasta');
	my $seqout_file_obj_2 = Bio::SeqIO->new(-file => ">>$fasta_file_name_2",-format => 'Fasta');
	
	if($pseudo_flag){
		my $seq_obj = Bio::Seq->new(
									-seq => $seq,
									-id  => $id."__".$pos."__".$te_num_pseudo{$id},
									-alphabet => 'dna');
		$seqout_file_obj_2->write_seq($seq_obj);
	}else{
		my $seq_obj = Bio::Seq->new(
									-seq => $seq,
									-id  => $id."__".$pos."__".$te_num{$id},
									-alphabet => 'dna');
		$seqout_file_obj->write_seq($seq_obj);
	}
}
##
	
## main ##

my @records = &get_one_record($file);

foreach(@records){
	&parse_exonerate_per_record($_);
}

