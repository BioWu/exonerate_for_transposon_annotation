#!/usr/bin/perl
 use File::Find::Rule;
 $dir = $ARGV[0];
 
 @subdir = File::Find::Rule->directory->in( $dir);
   my $rule =  File::Find::Rule->new;
  $rule->file;
  $rule->name( '*.pm' );
  my @files = $rule->in(@subdir);
#rm 12-species-pseudo-TEs-nums.txt
#echo "" >12-species-pseudo-TEs-nums.txt
#rm 12-species-actived-TEs-nums.txt
#echo "" >12-species-actived-TEs-nums.txt
#
#for i in `awk -F';' '{print $1"_"$5"_"$6}' /home/wucc-lab/PJ_TE-evolution/primate_all.usearch99.chr.size |awk '{print $1}'`
#	 do
#		echo $i|grep -P "\d" >>12-species-pseudo-TEs-nums.txt
#		echo $i|grep -P "\d" >>12-species-actived-TEs-nums.txt
#	 done;
#
#dir=$(ls -l ./ |awk '/^d/ {print $NF}')
#
#for j in $dir
#do
#    cd $j
#    rm 12-species-pseudo-TEs-nums.txt
#    rm 12-species-actived-TEs-nums.txt
#    for i in `awk -F';' '{print $1"_"$5"_"$6}' /home/wucc-lab/PJ_TE-evolution/primate_all.usearch99.chr.size |awk '{print $1}'`
#	 do 
#		echo  >>`basename $i`."pseudo.fa"
#		echo  >>`basename $i`."actived.fa"	
#	 done;
#    cp ../12-species-pseudo-TEs-nums.txt ../12-species-pseudo-TEs-nums.txt.tmp
#	
#    for i in *.pseudo.fa;do  grep ">" $i |wc -l >>temp.nums ;done;
#    paste  ../12-species-pseudo-TEs-nums.txt.tmp >../12-species-pseudo-TEs-nums.txt temp.nums; 
#    rm temp.nums
#    cd ..
#done
#
#rm 12-species-pseudo-TEs-nums.txt.tmp
#
#for j in $dir
#do
#    cd $j
#    cp ../12-species-actived-TEs-nums.txt ../12-species-actived-TEs-nums.txt.tmp
#    for i in *.actived.fa;do  grep ">" $i |wc -l >>temp.nums ;done;
#    paste ../12-species-actived-TEs-nums.txt.tmp temp.nums  >../12-species-actived-TEs-nums.txt; 
#    rm temp.nums
#    cd ..
#done
#rm 12-species-actived-TEs-nums.txt.tmp 
#
