# " perl blast_sperate.pl TE_list.txt *.blast"
#!/usr/bin/perl
use FileHandle;

open(TE,"<","$ARGV[0]");

open(BLAST,"<","$ARGV[1]");

while(<TE>){
	$a=$_;
	chomp($a);
	$a =~ tr/;/_/;
	print $a,"\n";
	open($fh{$a},">>","$a.blat");
}

while(<BLAST>){
	@inline = split("\t");
	$inline[0] =~ tr/;/_/;
	$inline[0] =~ tr/\\s//;
	$fh{$inline[0]}->print($_);
}
