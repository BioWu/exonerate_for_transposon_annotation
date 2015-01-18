#use warnings;
#use strict;

use String::Index qw( cindex ncindex crindex ncrindex );
use List::MoreUtils ':all';
 
my $file = "calJac.Step_2.get_cdna_v2.0.out.test";
#my $file = $ARGV[0];

sub get_one_record($file){
	local $/ = " >\n\n";
	my $file_handle = shift(@_);
	open(FH,$file_handle);
	my @records;
	while(<FH>){
		$_ =~ s/^\n//;
		push @records,$_;
		print "a\n";
	}
#	our @records = <FH>;
	return @records;
}
sub parse_exonerate_per_record($singlerecord){
	#Isertion => delete
	#Frameshift => delete
	#stopcodon => changes to NNN
	#with N => defined to pseudo-transposon
	
	my @record_lines = split('\n',shift(@_));

	# from(12,13) with stepsize 5；while( ! ($record_lines[12+5i] =~ /^#/ && $record_lines[13+5i] =~ /^#/) )
	my $i = 0;
	while( ! ($record_lines[12+5*$i] =~ /^#/ && $record_lines[13+5*$i] =~ /^#/) ){
		$i++;
#		my $original_target = $record_lines[10+5*$i];
		my $target = $record_lines[12+5*$i];
		my $query = $record_lines[13+5*$i];
		#for frameshift or stop codon or intron
		if($target =~ /\*|\#|-|\+/i,$query =~ /n|-|\.|\{|\}/i){
			my @target_ele = split('',$target);
			my @query_ele = split('',$query);
			#get removal elements indexs
			print "b\n";
			my @elements_will_be_removed_t = indexes {$_ =~ /\*|\#|-|\+/} @target_ele;
			my @elements_will_be_removed_q = indexes {$_ =~ /n|\.|-/i} @query_ele;
			
		}
		
	}
}

## main ##

my @records = &get_one_record($file);
foreach(@records){
	&parse_exonerate_per_record($_);
}

