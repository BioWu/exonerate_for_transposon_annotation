## exclude frame-shift || with percentage of N less than 10 % || without stop_codon in it [20..-20]aas.

open( FH, "<", "$ARGV[0]" );

open( FH_2, "<", "$ARGV[1]" );

use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::GenBank;
use Bio::PrimarySeq;

#open( FH_2, "<", "primate_all.usearch99.chr.size" );
#open( FH, "<","humrep_HERV4_I_Transposable_Element_FT_HERV4_I_2p_14357_2.calJac" );

sub has_stop_codon($seq,$f_length) {

	#get one seq try to test if it has stop codon
	my $seq                      = shift(@_);
	my $protein_identical_length = shift(@_);
	my @sequence                 = split( '#', $seq );

	my $seq_obj = Bio::Seq->new(
		-seq      => $sequence[1],
		-id       => $sequence[0],
		-alphabet => 'dna'
	);
	my $seq_length = $seq_obj->length();

	#suppose 1 based;subseq(1,$seq_length-60)->translate();
	my $protien = $seq_obj->trunc( 1, $seq_length - 60 )->translate()->seq();
	if ( $protien =~ /\*/ ) {
		return 1;
	}
	else {
		return 0;
	}

}

sub parse_merge_file() {
	while (<OUT_MERGE>) {

#		chomp;
#chr start end N_percentages true_postions scores strands TEs lengths sequences.
		my @inlines        = split("\t");
		my @true_positions = split( "&", $inlines[4] );    #1-based
		my @scores         = split( "&", $inlines[5] );
		chomp $inlines[9];
		my @sequences     = split( "&", $inlines[9] );
		my @TEs           = split( "&", $inlines[7] );
		my @length        = split( "&", $inlines[8] );
		my @N_percentages = split( "&", $inlines[3] );

		#find the index of the biggest element in an array.
		$scores[$idxMax] > $scores[$_] or $idxMax = $_ for 0 .. $#scores;
		my @true_postions_array = split( "@", $true_positions[$idxMax] );

		#		print @true_postions_array[1..3],"\n\n";
		$scores_all = join( ",", @scores );
		if ( $#scores >= 1 ) {
			print OUT_FINAL_BED join( "\t", @true_postions_array[ 1 .. 3 ] ),
			  "\t", $N_percentages[$idxMax], "\t", $scores[$idxMax], "\t",
			  $inlines[6], "\t", $TEs[$idxMax], "\t", $sequences[$idxMax], "\t",
			  $scores_all, "\n";
		}
	}
}

sub get_cdna_fasta() {
	my %te_num;

	#chr start end N_percentage best_score strand TEs fasta scores
	while (<IN_FINAL_BED>) {
		$i++;
		my @inline_bed = split("\t");
		my @fasta_info = split( /\s+/, $inline_bed[7] );

		#only coordinate info remained.
		my $fasta_desc = join( "-", @inline_bed[ 0 .. 2 ] );
		$fasta_desc =~ s/\>\>//;
		$fasta_desc =~ s/;/_/g;
		$fasta_desc =~ s/:/-/g;

		#ID description length#sequence
		my @sequence_fasta = split( "#", $fasta_info[-1] );
		$fasta_file_name   = $inline_bed[6] . ".fa";
		$fasta_file_name_2 = $inline_bed[6] . ".aa";
		if ( !defined( $te_num{ $inline_bed[6] } ) ) {
			$te_num{ $inline_bed[6] } = 0;
			open( NT, ">>", $fasta_file_name );
			close(NT);
			open( NT, ">>", $fasta_file_name_2 );
			close(NT);
		}
		$te_num{ $inline_bed[6] } += 1, "\n";

		my $seqout_file_obj =
		  Bio::SeqIO->new( -file => ">>$fasta_file_name", -format => 'Fasta' );
		my $seqout_file_obj_2 = Bio::SeqIO->new(
			-file   => ">>$fasta_file_name_2",
			-format => 'Fasta'
		);

#		my $seq_obj_bed = Bio::Seq->new (
#											  -seq => $sequence_fasta[-1],
#                                   -id  => $inline_bed[6].".".$te_num{$inline_bed[6]},
#                                   -accession_number => 'TE_'.$i,
#                                   -alphabet => 'dna',
#                                   -desc => $fasta_desc);
#       my $seq_obj_bed_2 = Bio::Seq->new (
#											  -seq => $seq_obj_bed->translate()->seq(),
#                                   -id  => $inline_bed[6].".".$te_num{$inline_bed[6]},
#                                   -accession_number => 'TE_'.$i,
#                                   -alphabet => 'dna',
#                                   -desc => $fasta_desc);
		my $seq_obj_bed = Bio::Seq->new(
			-seq => $sequence_fasta[-1],
			-id  => $inline_bed[6] . "."
			  . $te_num{ $inline_bed[6] } . "."
			  . $fasta_desc,
			-accession_number => 'TE_' . $i,
			-alphabet         => 'dna',
		);
		my $seq_obj_bed_2 = Bio::Seq->new(
			-seq => $seq_obj_bed->translate()->seq(),
			-id  => $inline_bed[6] . "."
			  . $te_num{ $inline_bed[6] } . "."
			  . $fasta_desc,
			-accession_number => 'TE_' . $i,
			-alphabet         => 'dna',
		);
		my $where = index( $seq_obj_bed_2->seq(), "*" );
		print $where, "\n";
		if ( $where > 0 ) {

			#				print $seq_obj_bed_2->seq();
			$seq_obj_bed_trunc = Bio::Seq->new(
				-seq => $seq,
				-id  => $inline_bed[6] . "."
				  . $te_num{ $inline_bed[6] } . "."
				  . $fasta_desc,
				-accession_number => 'TE_' . $i,
				-alphabet         => 'dna',
			);
			$seq_obj_bed_2_trunc = Bio::Seq->new(
				-seq => $seq_obj_bed_trunc->translate()->seq(),
				-id  => $inline_bed[6] . "."
				  . $te_num{ $inline_bed[6] } . "."
				  . $fasta_desc,
				-accession_number => 'TE_' . $i,
				-alphabet         => 'dna',
			);
		}
		if ( $where > 0 ) {
#			print "OOOOOOOO";
			$seqout_file_obj->write_seq($seq_obj_bed_trunc);
			$seqout_file_obj_2->write_seq($seq_obj_bed_2_trunc);

		}
		else {
			$seqout_file_obj->write_seq($seq_obj_bed);
			$seqout_file_obj_2->write_seq($seq_obj_bed_2);
		}
		undef $where;
	}
}

sub get_before_sort {
	open( OUT_SORT, ">", $ARGV[0] . ".b4sort" );
	my $frame_shift = 0;
	my $cdna_turn   = -1
	  ; #print "#chr\tN_percentage\tgene\tstart\tend\tscore\tstrand\tTE_ID\tTE_length\tcdna\n";
	while ( $in = <FH_2> ) {

		#index a hash for chr size
		chomp($in);
		@chr_sizes = split( "\t", $in );
		$chr_sizes_hash{ $chr_sizes[0] } = $chr_sizes[1];
	}

	while (<FH>) {
		if ( /exonerate:protein2genome:bestfit/ && /\t/ )
		{    #gff2 turn and frameshifts turn

			if (/frameshifts/) {
				$frame_shift = 1;
			}
			if (/\tgene\t/) {

				#			print $i++,"\n";
				@gff = split("\t");
				$gff[0] =~ /(\w*?\.?\w*?):(\d*?)-(\d*?)\(/;
				$gff[0] = $1;

				#calculate coverage of query sequences.
				$gff[3] += $2;
				$gff[4] += $2;
			}

		}
		elsif (/^>>/) {    #cdna-start end TURNING
			$cdna_turn *= -1;
		}
		elsif (/Query range:/) {
			$tmp = $_;
			chomp($tmp);
			@range = split( " -> ", $tmp );
		}
		elsif (/^ >/) {    #start output
			               #calculate N percentage
			$cdna_lines_num = scalar(@cdna) - 1;
			$cdna_seq = join( "", @cdna[ 1 .. $cdna_lines_num ] );
			my @N_matches = $cdna_seq =~ m/(N)/ig;
			my $N_count = scalar(@N_matches);
			$N_percentage = $N_count / ( length($cdna_seq) + 1 );

		 #      print "$N_percentage\t$cdna_seq\n";
		 #set N_percentage less than 10%
		 #set stop codon dosenot appear in seq except in last or begining 20 aa.
			$codon_turnning = &has_stop_codon( $cdna_seq, $fasta_length );

			# change to with frameshift and has stop codon
			if ( 1 )
			{

				#			print @gff,"\n";
				#		$gff[1] = ($range[1]-$range[0])/$fasta_length*1.0;
				$gff[1] = $N_percentage;
				( $gff[1], $gff[3] ) = ( $gff[3], $gff[1] );
				( $gff[2], $gff[4] ) = ( $gff[4], $gff[2] );
				( $gff[5], $gff[6] ) = ( $gff[6], $gff[5] );
				$gff[4] =
				  $gff[4] . "@" . $gff[0] . "@" . $gff[1] . "@" . $gff[2];
				print OUT_SORT join( "\t", @gff[ 0 .. 6 ] );

				#		$cdna[0] =~ /(\w*?);/;
				@id = split( ";", $cdna[0] );
				$id[0] =~ s/>>//;
				print OUT_SORT "\t", $id[0] . "_" . $id[3], "\t", $fasta_length,
				  "\t",, "\t", @cdna, "\n";
				undef @gff;
				undef @cdna;
				undef @range;
				undef $fasta_length;
				undef $N_percentage;
				$cdna_turn   = -1;
				$frame_shift = -1;

			}
			else {
				undef @gff;
				undef @cdna;
				undef @range;
				undef $fasta_length;
				undef $N_percentage;
				$cdna_turn   = -1;
				$frame_shift = -1;

			}
		}
		if ( $cdna_turn == 1 ) {
			if (/^>>/) {
				if ( $_ =~ />>(\S*)/ ) {
					$fasta_length = $chr_sizes_hash{$1};
					$dna_line     = $_;
					chomp($dna_line);
					push( @cdna, $dna_line );
					$bewteen_dna_seq = '#';
					push( @cdna, $bewteen_dna_seq );
				}
			}
			else {
				$dna_line = $_;
				chomp($dna_line);
				push( @cdna, $dna_line );
			}
		}

	}
	close(FH);
}

####main###

&get_before_sort();

$before_sorted_file = $ARGV[0] . ".b4sort";
$sorted_file        = $ARGV[0] . ".sort";
$merged_file        = $ARGV[0] . ".merge";
`bedtools sort -i $before_sorted_file >$sorted_file`;
`bedtools merge -i $sorted_file  -s -c 4,5,7,6,8,9,11 -o collapse,collapse,collapse,distinct,collapse,collapse,collapse -delim "&" >$merged_file`;

print STDERR "parsing merged file ! \n";

open( OUT_MERGE, "<", $ARGV[0] . ".merge" );

#open(OUT_FINAL_FASTA,">",$ARGV[0].".fa");
open( OUT_FINAL_BED, ">", $ARGV[0] . ".bed" );

&parse_merge_file();
close(OUT_FINAL_BED);
open( IN_FINAL_BED, "<", $ARGV[0] . ".bed" );
print STDERR "parsing bed file ! \n";

&get_cdna_fasta();

print STDERR "Done !\n";
