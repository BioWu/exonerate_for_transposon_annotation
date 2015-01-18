
use strict;
use warnings;
use Bio::SearchIO;
use Getopt::Long;

my ($format, $file,$output) = ('blast');

#GetOptions(
#	   'f|format:s'   => \$format,
#	   'i|input:s'    => \$file,
#	   'o|output:s'   => \$output);

#if( @ARGV ) { 
#    $file = shift;
#}
#    
#my $in = Bio::SearchIO->new(-format => $format,
#			    -file   => $file);

my $in = Bio::SearchIO->new(-file => 'humrep_HERV4_I_Transposable_Element_FT_HERV4_I_2p_14357_2.calJac',
                                 -format => 'exonerate');
my $out;
if( $output ) { 
    open $out, '>', $output or die "Could not write file '$output': $!\n";
} else { 
    $out = \*STDOUT;
}

while( my $r = $in->next_result ) { 
    while( my $hit = $r->next_hit ) {
	while( my $hsp = $hit->next_hsp ) {
	    my $mismatchcount = $hsp->length('total') - 
		($hsp->num_conserved + $hsp->gaps('total'));
	    print $out join("\t", ( $r->query_name,
				    $hit->name,
				    sprintf("%.2f",$hsp->percent_identity),
				    $hsp->length('total'),
				    $mismatchcount,
				    $hsp->gaps('total'),
				    # flip start/end on rev strand
				    $hsp->query->strand < 0 ?
				    ( $hsp->query->end,
				      $hsp->query->start ) :
				    ( $hsp->query->start,
				      $hsp->query->end ),
				    $hsp->hit->strand < 0 ?
				    ( $hsp->hit->end,
				      $hsp->hit->start ) :
				    ( $hsp->hit->start,
				      $hsp->hit->end ),

				    $hsp->evalue,
				    # chance this to $hsp->sw_score 
				    # if you would rather have that
				    # it will only work for FASTA parsing though!
				    $hsp->bits)),"\n";
	}
    }
}