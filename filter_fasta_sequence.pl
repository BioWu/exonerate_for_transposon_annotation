use Bio::SeqIO;
use List::Util qw(sum);
use POSIX qw(ceil);
$fasta_file = $ARGV[0];
$fasta_file_2 = $ARGV[1];
my $seq_file = Bio::SeqIO->new(-file => $fasta_file,-format => fasta);
my @lengths;
while($seq_obj = $seq_file->next_seq ){
	if($seq_obj->length()%3 != 0){print $seq_obj->id,  "\tPartial\n";}
	push @lengths , $seq_obj->length();
}

$ave_len = ceil(sum(@length)/($#lengths+1));

$seq_file = Bio::SeqIO->new(-file => $fasta_file,-format => 'fasta');
$seq_file_out = Bio::SeqIO->new(-file => ">$fasta_file_2",-format => 'fasta');
while($seq_obj = $seq_file->next_seq ){
	if($seq_obj->length() > $ave_len-1){
		$seq_file_out->write_seq($seq_obj);
	}
}