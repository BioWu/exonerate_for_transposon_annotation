#use warnings;
#use strict;

# should notice "** (process:25273): WARNING **: Exhaustively generating suboptimal alignments will be VERY SLOW"
use String::Index qw( cindex ncindex crindex ncrindex );
use List::MoreUtils ':all';
use Bio::Seq;
use Bio::SeqIO;

#my $file = "hg19.Step_2.get_cdna_v2.0.out.test";
my $file = "$ARGV[0]";
my %te_num;
my %te_num_pseudo;

sub from_bed_to_fa(@bed_record) {
	$tem = shift(@_);
	my @bed_record = split( "\t", $tem );

	my $pseudo_flag = $bed_record[3];

	my $id = $bed_record[4];

	my $seq = $bed_record[6];

	my $pos =
	    $bed_record[0] . "_"
	  . $bed_record[1] . "_"
	  . $bed_record[2] . "_"
	  . $bed_record[5];

	my $raw_score = $bed_record[7];
	#### for sequence output
	if ( !( defined( $te_num{$id} ) || defined( $te_num_pseudo{$id} ) ) ) {

		#		print $id;
		$fasta_file_name   = $id . ".actived.fa";
		$fasta_file_name_2 = $id . ".pseudo.fa";
		$te_num{$id}       = 0;
		open( NT, ">>", $fasta_file_name );
		close(NT);
		open( NT, ">>", $fasta_file_name_2 );
		close(NT);
	}
	$fasta_file_name   = $id . ".actived.fa";
	$fasta_file_name_2 = $id . ".pseudo.fa";
	my $seqout_file_obj =
	  Bio::SeqIO->new( -file => ">>$fasta_file_name", -format => 'Fasta' );
	my $seqout_file_obj_2 =
	  Bio::SeqIO->new( -file => ">>$fasta_file_name_2", -format => 'Fasta' );

	if ($pseudo_flag) {
		$te_num_pseudo{$id} += 1;
		my $seq_obj = Bio::Seq->new(
			-seq      => $seq,
			-id       => $id . "__" . $pos . "__" . $te_num_pseudo{$id},
			-alphabet => 'dna',
			-desc     => $raw_score
		);
		$seqout_file_obj_2->write_seq($seq_obj);
	}
	else {
		$te_num{$id} += 1;
		my $seq_obj = Bio::Seq->new(
			-seq      => $seq,
			-id       => $id . "__" . $pos . "__" . $te_num{$id},
			-alphabet => 'dna',
			-desc     => $raw_score
		);
		$seqout_file_obj->write_seq($seq_obj);
	}
}

sub deal_with_same_length_in_overlap(@selected_TE,@max_length_index,@scores) {
	my ( @selected_TE, @max_length_index, @scores ) = @_;
	my ( $score_min, $score_max ) =
	  minmax( @scores[ @active_TE[@max_length_index] ] );
	my @max_score_index =
	  indexes { $_ =~ $length_max } @scores[ @active_TE[@max_length_index] ];
	return @max_score_index;
}

sub parse_bed($file_bed) {

	#2015-3-15:Add Raw score to correct the selection of overlap regions.
	open BED, "<", shift(@_);
	while (<BED>) {
		my @beds = split();
		if (/\&/) {
			chomp;

			# has overlap
			my @chrs        = split( "&", $beds[3] );
			my @starts      = split( "&", $beds[4] );
			my @ends        = split( "&", $beds[5] );
			my @pseud_flags = split( "&", $beds[6] );
			my @TE_ids      = split( "&", $beds[7] );
			my $strand      = $beds[8];
			my @seqs        = split( "&", $beds[9] );
			my @scores      = split( "&", $beds[10] );
			my @lengths     = split( "&", $beds[11] );
			## Activty first then CDS length then best score
#			my @elements_will_be_removed_t = indexes {$_ =~ /\#|-|\+|\{|\}/} @target_ele;
			my @active_TE = indexes { $_ =~ /0/ } @pseud_flags;
			if ( any { $_ >= 0 } @active_TE ) {
				print STDERR "has active TEs\n";
				my ( $length_min, $length_max ) =
				  minmax( @lengths[@active_TE] );
				my @max_length_index =
				  indexes { $_ =~ $length_max } @lengths[@active_TE];

				if ( defined $max_length_index[1] ) {
					print STDERR
					  "at least two same-length TE overlaped and active. ",
					  "\n";
					### I choose the first one
					my @max_score_index =
					  &deal_with_same_length_in_overlap( @active_TE,
						@max_length_index, @scores );
					if ( defined $max_score_index[1] ) {
						print STDERR
"at least two same-length && same-score TE overlaped and active. ",
						  "\n";
					}
					$bed_index = $active_TE[ $max_score_index[0] ];
					$bed       = join( "\t",
						$chrs[$bed_index],   $starts[$bed_index],
						$ends[$bed_index],   $pseud_flags[$bed_index],
						$TE_ids[$bed_index], $strand,
						$seqs[$bed_index],   $scores[$bed_index] );

		  # $bed ::: chrY	25291721	25294365	1	prirep_L1-1_Cja_32538_1	+	agggggaa
					&from_bed_to_fa($bed);
					print BED_OUT_FINAL $bed,"\n";
				}
				else {
					$bed_index = $active_TE[ $max_length_index[0] ];
					$bed       = join( "\t",
						$chrs[$bed_index],   $starts[$bed_index],
						$ends[$bed_index],   $pseud_flags[$bed_index],
						$TE_ids[$bed_index], $strand,
						$seqs[$bed_index],   $scores[$bed_index] );

		  # $bed ::: chrY	25291721	25294365	1	prirep_L1-1_Cja_32538_1	+	agggggaa
					&from_bed_to_fa($bed);
					print BED_OUT_FINAL $bed,"\n";
				}
			}
			else {
				print STDERR "has Pseudo TEs\n";
				my ( $length_min, $length_max ) =
				  minmax( @lengths[@active_TE] );
				my @max_length_index =
				  indexes { $_ =~ $length_max } @lengths[@active_TE];

				if ( defined $max_length_index[1] ) {
					print STDERR
					  "at least two same-length TE overlaped and active. ",
					  "\n";
					### I choose the first one
					my @max_score_index =
					  &deal_with_same_length_in_overlap( @active_TE,
						@max_length_index, @scores );
					if ( defined $max_score_index[1] ) {
						print STDERR
"at least two same-length && same-score TE overlaped and active. ",
						  "\n";
					}
					$bed_index = $active_TE[ $max_score_index[0] ];
					$bed       = join( "\t",
						$chrs[$bed_index],   $starts[$bed_index],
						$ends[$bed_index],   $pseud_flags[$bed_index],
						$TE_ids[$bed_index], $strand,
						$seqs[$bed_index],   $scores[$bed_index] ),
					  "\n";

		  # $bed ::: chrY	25291721	25294365	1	prirep_L1-1_Cja_32538_1	+	agggggaa
					&from_bed_to_fa($bed);
					print BED_OUT_FINAL $bed,"\n";
				}
				else {
					$bed_index = $active_TE[ $max_length_index[0] ];
					$bed       = join( "\t",
						$chrs[$bed_index],   $starts[$bed_index],
						$ends[$bed_index],   $pseud_flags[$bed_index],
						$TE_ids[$bed_index], $strand,
						$seqs[$bed_index],   $scores[$bed_index] );

		  # $bed ::: chrY	25291721	25294365	1	prirep_L1-1_Cja_32538_1	+	agggggaa
					&from_bed_to_fa($bed);
					print BED_OUT_FINAL $bed,"\n";
				}
			}

		}
		else {

			# no_overlap
			chomp;
			$bed = join( "\t", @beds[ 3 .. 10 ] );

		  # $bed ::: chrY	25291721	25294365	1	prirep_L1-1_Cja_32538_1	+	agggggaa
		  #			print $bed,"\n";
			&from_bed_to_fa($bed);
			print BED_OUT_FINAL join( "\t",
				@beds[ 3 .. 5 ],
				$beds[7], $beds[6], @beds[ 8 .. 10 ] ),
			  "\n";
		}
	}
}

sub get_bed_file($chr,$start,$end,$pseudo_flag,$TE_id,$strand,$sequence) {
	## get one seq_obj then return one bed-record
	my @bed = @_;
	$bed[8] = length( $bed[6] );
	return join( "\t", @bed );
}

sub get_one_record($file) {
	local $/ = " >\n\n";
	my $file_handle = shift(@_);
	open( FH, $file_handle ) or die 'can not open $file\n';
	my @records;
	while (<FH>) {
		$tmp = $_;
		$tmp =~ s/^\n//;
		push @records, $tmp;
	}

	#	our @records = <FH>;
	return @records;
}

sub parse_exonerate_per_record($singlerecord) {
	my @sequence;
	my @record_lines = split( '\n', shift(@_) );

	#Isertion => delete
	#Frameshift => delete
	#stopcodon => changes to NNN
	#with N => defined to pseudo-transposon
	#if it is an psudo gene please note it on id
	## Query: humrep_L1PREC2;Transposable_Element;FT;L1PREC2_2p;14809;2
	$record_lines[2] =~ /:\s(.*?);.*?;(\w*);(\w*$)/;
	my $id = join( "_", $1, $2, $3 );
	$raw_score = 0;
	if ( $record_lines[5] =~ /:\s(\d+)/ ) {
		$raw_score = $1;
	}
	else {
		die "Can not get the Raw score!\n";
	}

#0 for actived ; 1 for frameshift ; 2 for stop_codon ; 3 for N ; 4 for frameshift&stopcodon
	$pseudo_flag = 0;
	## offset Target range: 1 -> 3635
	if ( $record_lines[7] =~ /:\s(\d*)\s\S+\s(\d*)/ ) {
		$tmp_offset[0] = $1;
		$tmp_offset[1] = $2;

		#if complement to target sequence, then throw this sequence .
		if ( $tmp_offset[0] > $tmp_offset[1] ) {

#which means this sequence is another TE BLATed ,but is reverse complemented by the same one;
			print STDERR join( "\n", @record_lines ), "\n";
			return 'Complement';
		}
	}
	else {
		die "Can not get target-ranges.\n";
	}
	## Target: chr1:73900-77690(+)
	$record_lines[3] =~ /:\s(.*?):(\w*)-(\w*)\(([+|-])\)/;
	## Should  noticed $strand if '-'
	if ( $4 eq '+' ) {
		( $chr, $start, $end, $strand ) =
		  ( $1, $2 + $tmp_offset[0], $2 + $tmp_offset[1], $4 );
	}
	elsif ( $4 eq '-' ) {
		( $chr, $start, $end, $strand ) =
		  ( $1, $3 - $tmp_offset[1], $3 - $tmp_offset[0], $4 );
	}
	else {
		die "Can not get strand-info. \n";
	}

	#	print join("\t",$chr,$start,$end,$strand),"\n";
	my $pos = join( '_', $chr, $start, $end, $strand );

# from(12,13) with stepsize 5；while( ! ($record_lines[12+5i] =~ /^#/ && $record_lines[13+5i] =~ /^#/) )
	my $i = 0;
	$N_num = 0;
	while (
		!(
			   $record_lines[ 11 + 5 * $i ] =~ /^#/
			&& $record_lines[ 12 + 5 * $i ] =~ /^#/
		)
	  )
	{

		#		my $original_target = $record_lines[10+5*$i];
		my $target = $record_lines[ 11 + 5 * $i ];
		$target =~ /^\s*(.*?)\s*$/;
		$target = $1;
		$target =~ s/" "/"-"/g;
		my $query = $record_lines[ 12 + 5 * $i ];
		$query =~ /:\s*(.*?)\s*:/;
		$query = $1;
		@N_now_query = split("",$query);
		@N_selected = indexes {$_ =~ /n/i} @N_now_query;
		$N_num += $#N_selected + 1; 
		if($N_num > $allowed_N){
			return 'Too many unknown nucliteds';
		}
		$i++;
		$add = abs( length($target) - length($query) );

		#for frameshift or stop codon or intron

		if ( length($target) < length($query)
			&& ( $query =~ /^[a|t|c|g|n|\-|\{|\}]/i ) )
		{

		#"-" for sth like this "
		# --SerGlnGluLeuArgPheGlyAsnLeuCysLeuAspPheArgArgCysMet+-
		# --agtcaagaattgaggtttgggaacctctgcctagatttcagaagatgtatgga....... : 3249"
		#"-" start from left or may rigth
			print "Start R\t$id\n";

			#only from left
			$target = $target . ( '-' x $add );
		}
		elsif ( length($target) < length($query) ) {
			if ( ( $query =~ /^(\.)/ ) && ( $query =~ /(\.)$/ ) ) {

				#"-" start also from right(etc. both left and right)
				print "Start R && L\t$id\n";
				my @query_ele = split( '', $query );
				my @add_index = indexes { !( $_ =~ /\./i ) } @query_ele;
				my $add_left = $add_index[0];
				$target = ( '-' x $add_left ) . $target;
				$target = $target . ( '-' x ( $add - $add_left ) );
			}
			else {

				# '-'start from right
				print "Start L\t$id\n";
				$target = ( '-' x $add ) . $target;
			}

		}
		if ( length($target) == length($query)
			&& ( $target =~ /\*|\#|\-|\+/i || $query =~ /n|\-|\.|\{|\}|N/i ) )
		{

#small bugs, I choose to correct pseudo_flag accordding target and query sequences
			if ( $target =~ /\*|\#/i || $query =~ /n/i ) {
				$pseudo_flag = 1;
			}
			my @target_ele = split( '', $target );
			my @query_ele  = split( '', $query );

			#get removal elements indexs
			#for frame_shift || N
			# Notice:: for N sites , we choose to keep
			my @elements_will_be_removed_t =
			  indexes { $_ =~ /\#|\-|\+|\{|\}/i } @target_ele;

			my @elements_will_be_changed_t =
			  indexes { $_ =~ /\*/i } @target_ele;

			# I choose to keep N sites(unknown sites).
			my @elements_will_be_removed_q =
			  indexes { $_ =~ /\.|\-|\{|\}/i } @query_ele;
			push my @elements_will_be_removed_mix, @elements_will_be_removed_t,
			  @elements_will_be_removed_q;
			@elements_will_be_removed_mix_uniq =
			  uniq @elements_will_be_removed_mix;

#			print join("\t",@elements_will_be_removed_mix_uniq,@elements_will_be_changed_t),"\n";
			foreach (@elements_will_be_changed_t) {
				$query_ele[$_] = 'N';
			}
			foreach (@elements_will_be_removed_mix_uniq) {
				$query_ele[$_] = '';
			}

			#			print "Raw:";
			#			print $query,"\nNow:";
			#			print join('',@query_ele),"\n";
			$query_2 = join( '', @query_ele );
			if ( $query_2 =~ /\./ ) { print $query; die "BUG\n"; }
		}
		elsif ( length($target) == length($query) ) {
			$query_2 = $query;
		}
		else {
			print @record_lines;
			die "BUG:the Target does not equal to Query lines.\n";
		}
		push @sequence, $query_2;

	}
	my $seq = join( '', @sequence ), "\n";
	##### for bed output
	return &get_bed_file( $chr, $start, $end, $pseudo_flag, $id, $strand, $seq,
		$raw_score );

}
##

## main ##
my $allowed_N = 20;
my $file_b4sort    = $file . ".b4sort";
my $file_b4merge   = $file . ".b4merge";
my $file_bed       = $file . ".bed";
my $file_final_bed = $file . "final.bed";

my @records = &get_one_record($file);

`rm $file_b4sort`;

open B4SORT, ">>", $file_b4sort;
print STDERR "getting bed file .\n";
open BED_OUT_FINAL, ">", $file_final_bed;
foreach (@records) {

	# for debug only
	$return_s = &parse_exonerate_per_record($_);
	if ( $return_s eq 'Too many unknown nucliteds' ) {
		print STDERR "\nToo many unknown nucliteds !\n";
	}
	elsif($return_s eq 'Complement') {
		print STDERR "\nComplement to target sequence !\n";
	}else{
		print B4SORT &parse_exonerate_per_record($_), "\n";
	}
}
`bedtools sort -i $file_b4sort >$file_b4merge`;
`bedtools merge -i $file_b4merge  -s -c 1,2,3,4,5,6,7,8,9 -o collapse,collapse,collapse,collapse,collapse,distinct,collapse,collapse,collapse -delim "&" >$file_bed`;
`rm *.fa `;
print STDERR 'pasrsing bed .\n';

&parse_bed($file_bed);

