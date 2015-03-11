
use List::Util qw(sum);

sub re_to_array($x) {
	my $b = shift(@_);
	my @array;
	while ( $b =~ /(\d+)[M|D]/ ) {
		push @array, $1;
		$b = $';
	}
	return @array;
}

sub re_to_array_N($x) {
	my $b = shift(@_);
	my @array;
	while ( $b =~ /(\d+)N/ ) {
		push @array, $1;
		$b = $';
	}
	return @array;

}

sub parse_CIGAR($a) {
	my $a = shift(@_);
	if ( $a =~ /N/ ) {
		@Ns = &re_to_array_N($a);
		@b = split( /\d*N/, $a );
		foreach (@b) {
			$i++;
			push @CIGAR_2, sum &re_to_array($_);
			push @CIGAR_2, $Ns[ $i - 1 ];

		}
	}else{
	push @CIGAR_2, sum &re_to_array($a);
	}
	return @CIGAR_2;
}

## MAIN ##
open(FH,"$ARGV[0]");
#$a = "45S16M568N16M7D566N16M3D45I16M568N16M7D566N16M3S";
#$b = "564M4D7I6576M";
while(<FH>){
	chomp;
@inlines =split("\t");
@CIGAR = &parse_CIGAR($inlines[2]);
print $_,"\t",join( "\t", @CIGAR ),"\n";
}


