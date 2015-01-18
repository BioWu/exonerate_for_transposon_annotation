#" perl ./run_exonerate_temp.pl TE.ID Protein_path Dna_path out_path"
open ID, "$ARGV[0]";

while (<ID>) {
	$id = $_;
	chomp($id);
	print $id;
	$query  = $ARGV[1] . "/" . $id . ".fsa";
	$target = $ARGV[2] . "/" . $id . ".blat.sort.bed.bed_2.fa";
	@path   = split( "/", $target );
	$out    = $ARGV[3] . "/" . $id . $path[-1];
`exonerate  --model protein2genome:bestfit --exhaustive y -g y --verbose 0 --showtargetgff y $query $target  --ryo ">>%qi %ti %tcb %tce %tcl\n%tcs >\n" > $out`;
}

# * print out the 'gene' line in GFF3
# * prints out the 'exon' lines in GFF3, link to parent gene
#
# INPUT FILE: exonerate should have been run with (at least) the following options:
# --model coding2genome OR --model est2genome
# --showtargetgff yes
#
# it is OK to have other output formats in the exonerate output
#
#end

### for coordinater :: start+a-1;start+b then use bedtools get fasta get true sequence
##angapp_RTE-2_AG;Transposable_Element;FT;RTE-2_AG_2p;21;2    445    Aech_gn2.0_scaffold185    gene    379259    382419    127    -    .    gene_id 1 ; sequence angapp_21_2 ; gene_orientation +    angapp_21_2_1
use Getopt::Long;
my $help;
GetOptions( "help" => \$help );
usage() if ($help);
my $est_acc;

# my $alignment_serial_number = 0;
open FH, "$ARGV[0]";
while (<FH>) {
	chomp;
	my @cols = split /\t/, $_;
	if ( /exonerate:\w+2genome/ && !/^"#"/ ) {
		if ( /\tgene\t/ && /sequence (\S+)/ ) {
			$cols_2[0] = $1;

			if ( !defined $turn{ $cols[0] } ) {
				$turn{ $cols[0] } = 1;
			}
			else {
				$turn{ $cols[0] } = 0;
			}
		}
		if($turn{ $cols[0] } ==1){
		$cols[0] =~ /(\w*\d*):(\d+)-(\d+)/;

		$cols_2[1] = $1;
		$cols_2[2] = $cols[2];
		$cols_2[3] = $2 + $col[3] - 1;
		$cols_2[4] = $3 + $col[4];
		$cols_2[5] = $col[5];
		$cols_2[6] = $col[6];
		$cols_2[7] = $col[7];
		print join "\t", @cols_2;
		print "\n";
	}
	}
}

sub usage {
	if ( open( ME, $0 ) ) {
		my $ignore = <ME>;
		while (<ME>) {
			last if (/ ^ "#" + "end" /);
			if (s/^"#"+//) {
				print;
			}
		}
		close(ME);
	}
	exit;
}
