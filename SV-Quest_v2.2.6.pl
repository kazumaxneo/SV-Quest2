#!/usr/bin/perl
#mamba install -c bioconda mira spades bwa samtools minimap2 seqkit sambamba blast -y
#リードサポート4以下はカウントしないように変更($coverage2 1 => 4)
#v2.1 Reverted back to the old version (;Sam output failed to consider lead orientation).
#v2.2 Fixed to retrieve reads considering direction from F and R respectively at local assembly and to retrieve  the best hit in blastnoutfmt 7.
#v2.2.2 Fixed size cutoff to 200 => 150 (;to rescue very short target seqs like 187-bp) 
#v2.2.3 spades option changed to use small k-mer values.
#v2.2.5 Changed the extraction range from 500 bp to 250 bp. Spades run without error collection



use strict;
use List::Util qw(max);
use Getopt::Long qw(:config posix_default no_ignore_case gnu_compat);

&title;

my $fastq1 = "";
my $fastq2 = "";
my $bam = "";
my $fasta = "";
my $threshould = 0.1;
my $size = 25;
my $call_threshould = "0";
my $coverage = 0.7;
my $coverage2 = 4;
my $outputdir = "summary";
my $output = "InDel.txt";
my $max = "100000";
my $check = "1";
my $wide = "12";
my $Lborder = "300";
my $cpu = "8";
my $seed = "25";
my $seedpenalty  = "1";
my $loose = "275";
my $lowcoverage = "0.1";
my $duplicateRemove = "1";
my $map = "no";
my $compress = "2";
my $radius = "0.65";
my $slope_threshould = "3";
my $Tn = "";
my $samplename = "sample";
my $loop  = "1";

GetOptions('i=s' => \$bam, 'f=s' => \$fasta, 'b=s' => \$Tn, 'm=f' => \$threshould, 'c=i' => \$size, 'h=i' => \$loop, 't=f' => \$call_threshould, 'C=f' => \$coverage, 'y=i' => \$coverage2, 'o=s' => \$output, 'O=s' => \$outputdir, 'n=s' => \$samplename, 'L=i' => \$max, 'r=f' => \$check, 'w=i' => \$wide, 'U=i' => \$Lborder, 'p=i' => \$cpu, 's=i' => \$seed, 'M=i' => \$seedpenalty, 'e=i' => \$loose,'1=s' => \$fastq1,'2=s' => \$fastq2,'cov=f' => \$lowcoverage,'d=i' => \$duplicateRemove,'k=s' => \$map,'z=i' => \$compress, 'R=f' => \$radius, 'S=i' => \$slope_threshould);

die "\nMapped.bam or fastq file is required !\n\n\n" if($bam eq "" && $fastq1 eq "" && $check == 1);
die "\ninput fasta fileis required !\n\n\n" if($fasta eq "");
die "\nThe threshould must be 0-1 real number ! default 0.3\n\n\n" unless($threshould <= 1 or $threshould >= 0);
die "\nThe call_size must be plus Integer ! default 30\n\n\n" unless($size >= 1);
die "\nThe call threshould must be plus Integer ! default 30\n\n\n" unless($call_threshould >= 0);
die "\nThe coverage threshould for calling deletion must be 0-1 real number ! default 0.3\n\n\n" unless($coverage <= 1 or $coverage >= 0);
die "\nThe maximum length for searching deletion must be greater than 101 ! default 100000\n\n\n" unless($max > 101);

print " - Mode is $check\n - Threshould is $threshould\n - Wides is $wide\n - Coverage threshould $coverage\n\n";
system("mkdir temp $outputdir");
open TEMP0, ">temp/temp0" or die "cant save \n"&&exit;
print TEMP0 "$fasta";
open TEMP2, ">temp/temp2" or die "cant save \n"&&exit;
print TEMP2 "$size";
open TEMP8, ">temp/temp8" or die "cant save \n"&&exit;
print TEMP8 "$output";
close TEMP0;close TEMP2;close TEMP8;

&return if($check == 1);
system("zcat $fastq1 > temp/F.fastq") if($fastq1 =~ /.gz$/);
system("zcat $fastq2 > temp/R.fastq") if($fastq2 =~ /.gz$/);
system("cat $fastq1 > temp/F.fastq") unless($fastq1 =~ /.gz$/);
system("cat $fastq2 > temp/R.fastq") unless($fastq2 =~ /.gz$/);
&count if($check == 1);
my @now = localtime;print "INFO $now[2]:$now[1]:$now[0]\tProcessing fasta\n" if($check == 1);
&mapping if($check == 1 && $bam eq "");
&sam_convert if($fastq1 eq "" && $check == 1);
system("rm $output") if($check == 1);
system("rm $output $outputdir/log.txt 6F.txt 6R.txt F_mismatch_mountain.txt R_mismatch_mountain.txt temp/F_mismatch_mountain_all.txt temp/R_mismatch_mountain_all.txt $output") unless($check == 1);
&split if($check == 1);#split the reads by the orientation of the reads
my @now = localtime;print "INFO $now[2]:$now[1]:$now[0]\tStarts samtools sort of F.sam\n" if($check == 1);
system("samtools fastq -@ $cpu temp/F.sam > temp/F.fastq") if($check == 1);
system("samtools fastq -@ $cpu temp/R.sam > temp/R.fastq") if($check == 1);
system("samtools sort -@ $cpu -O BAM temp/F.sam > temp/F_sorted.bam") if($check == 1);
my @now = localtime;print "INFO $now[2]:$now[1]:$now[0]\tStarts samtools sort of R.sam\n" if($check == 1);
system("samtools sort -@ $cpu -O BAM temp/R.sam > temp/R_sorted.bam") if($check == 1);
my @now = localtime;print "INFO $now[2]:$now[1]:$now[0]\tStarts samtools index of F.bam\n" if($check == 1);
system("samtools index temp/F_sorted.bam") if($check == 1);
my @now = localtime;print "INFO $now[2]:$now[1]:$now[0]\tStarts samtools index of R.bam\n" if($check == 1);
system("samtools index temp/R_sorted.bam") if($check == 1);
system("rm temp/F.sam temp/R.sam temp/input.sam") if($check == 1);
my @now = localtime;print "INFO $now[2]:$now[1]:$now[0]\t" if($check == 1);
system("sambamba mpileup temp/F_sorted.bam -t $cpu --samtools -f temp/reference.fa -O -Q 0.1 2> /dev/null > temp/F_mpileup") if($check == 1 or $check == 3);
my @now = localtime;print "INFO $now[2]:$now[1]:$now[0]\t" if($check == 1 or $check == 3);
system("sambamba mpileup temp/R_sorted.bam -t $cpu --samtools -f temp/reference.fa -O -Q 0.1 2> /dev/null > temp/R_mpileup") if($check == 1 or $check == 3);
my @now = localtime;print "INFO $now[2]:$now[1]:$now[0]\tStarts calculation of F&R mismatch\n";

#variant vcf
#system("echo \"\#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample\" |sed s/\' \'/$\'\t\'/g > $outputdir/variant1.vcf");
#system("echo \"\#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample\" |sed s/\' \'/$\'\t\'/g > $outputdir/variant2.vcf");
system("echo \"\#sample chr start end mismatch_Score name identity\" |sed s/\' \'/$\'\t\'/g > $outputdir/target.txt");
#target blastDB
system("mkdir temp/blastDB/ && makeblastdb -dbtype nucl -in $Tn -out temp/blastDB/targetDB 2> /dev/null");

&extract;#calclulate mismatch base from mpileup
system("mv InDel.txt $outputdir") if($map eq "no");

my @now = localtime;print "\nINFO $now[2]:$now[1]:$now[0]\tProgram end.\n";
system("rm temp/temp* temp/intrachromosomal1.fq temp/mixed1.fasta temp/target1.fa temp/targetextension1.fa temp/F_mpileup temp/R_mpileup temp/R.fastq temp/F.fastq .non-redundant.txt temp/*_mpileup temp/mismatch* temp/reference.fa.* temp/fasta* temp/*.fastq* temp/In.txt temp/Del.txt temp/*redundant.txt");
system("rm temp/fasta rm temp/fasta.index rm temp/reference.fa.sa rm temp/reference.fa.pac rm temp/reference.fa.fai rm temp/reference.fa.bwt rm temp/reference.fa.ann rm temp/reference.fa.amb");

system("mkdir temp/local_assembly && mv temp/spades* temp/local_assembly");
system("rm -rf temp/local_assembly/spades*/pipeline_state temp/local_assembly/spades*/assembly_graph_after_simplification.gfa");
system("rm temp/local_assembly/spades*/run_spades.*");
system("mkdir temp/vcf && mv temp/*vcf temp/vcf");
system("zip -r $output.zip */ InDel.$output $outputdir/log.txt circos* *html $output") if($compress == 1);
print "\nRessults are saved as file name \"$output\"\n\n" unless($compress == 1);
print "\nRessults are saved as file name \"$output\.zip\"\n\n" if($compress == 1);
exit;


#######################################################################################################################################################################################
### SUBROUTINES ###

sub title {
	####################################################################################################################################
	#
	#						SV-Quest version 2.0 #仮
	#
	#						Kazuma Uesaka
	#						University of Nagoya
	#						9 February 2018
	#		
	#						A Perl scripts to call SV position from mapped.bam.
	#
	#	
	#		SV Quest: Sensitive Structure Variation detection tool.
	#		Kazuma Uesaka, Hiroshi Kouno, Kazuki Terauchi, Yuichi Fujita, Tatsuo Omata, and Kunio Ihara
	#
	#		
	#						Input:
	#							bam file and reference.fasta for the mapping
	#
	#						Outnput:
	#							Insertion and deletion position printed to STDOUT
	#						
	#						Usage:
	#						perl SV-Quest.pl -f reference.fa -1 forward.fq -2 reverse.fq -b transposon.fasta -n sample1 -y 4
	#
	#				The mapped.bam and it's reference.fasta should be included in the same folder.
	#
	#
	#
	#	ver0.4 2017-02-23 version0.4 solid mapped.bam are supported.illumina paired-end fastq are supported.
	#	ver0.5 2017-05-15 coverage ratio from 0.3 to 0.7
	#	ver0.6 check coverage (dont call low coverage region)
	#	ver0.7 mismatch mountain height supported.
	#	ver0.8 redundant call removed.
	#	ver0.9 sambamba mplileup was introduced.
	#	ver2.0 local assembly supported.
	####################################################################################################################################


	print "\n\n############################################################################################################################################################\n";
	print "Program: SV-Quest\n";
	print "version 2.0\n\n";
	print "\nUsage:	perl SV-Quest.pl -f reference.fa -1 forward.fq -2 reverse.fq -b transposon.fasta <options>\n\n";
	print "Input/output options:\n\n";
	print "\t -1	input fastq (Required if .bam file doesn't assigned)\n";
	print "\t -2	input pair fastq (for paired end reads)\n";
	print "\t -i	input BAM (Required if .fastq file doesn't assigned)\n";
	print "\t -f	input fasta (Required)\n";
	print "\t -b	target squence (Required)\n";
	print "\t -o	output file name (default Indel.txt)\n";
	print "\t -O	output file directory (default summary)\n";
	print "\t -n	sample name (default sample)\n";

	print "\nMismatch calclulation options:\n\n";
	print "\t -y	Threshould of mismatch ratio coverage. Dont count mismatch ratio if supprting reads are smaller than value (default 4)\n";
	print "\t -C	Threshould of coverage  (default 0.7; stop searching when the candidate read coverage are bigger than 70% of genomic average coverage.)\n";
	print "\t -L	Maximum length (bp) for searching deletion  (default 100000)\n";
	print "\t -U	The border width (bp) of Pf and Pr to call insertion or deletion (default average read length)\n";
	print "\t -c	Calculation size (bp) of mismatch value (default 25)\n";
	print "\t -m	Threshould of mismatch ratio (default 0.1; stop searching when the mismatch ratio was bigger than 0.1)\n";
	print "\t -t	Threshould of call (default 0)\n\n";
	print "\t -S	Mismatch mountain slope threshould to call indel (default 3)\n";
	print "\t -r	Run or skip mapping analysis (default 1)\n";
	print "\t			1;runs\n";
	print "\t			2;skip (re-perform SV-Quest with same mapping parameter. use previous mpileup results)\n";
	print "\t			3;skip and perform mpileup (re-perform SV-Quest with same mapping parameter, but do mpileup)\n";	
	print "\t -w	The size of mismatch value (bp)  (default 12: 12 means ± 12bp mismatch at each postion are included)\n";
	print "\t -p	CPU thread (default 20)\n";
	print "\t -s	seed length (default 25)\n";
	print "\t -M	maximum differences in the seed (default 1)\n";
	print "\t -e	do not put an indel within INT bp towards the ends (default 275)\n";
	print "\t -h	number of local assembly (default 1)\n";
	print "\t -d	remove duplicate call (default 1)\n";
	print "\t			1;runs\n";
	print "\t			2;skip\n";
	print "\t -k	draw indel map using circos (default yes)\n";
	print "\t			yes;runs\n";
	print "\t			no;skip\n";
	print "\t -R	diameter of map (default 0.65)\n";
	print "\t -z	compress results and related folder with zip (default 2)\n";
	print "\t			1;runs\n";
	print "\t			2;skip\n";
	#print "\t -cov	remove low coverage region. 0.1 means the possible Pf/Pr that have 1/10 times small than average coverage are removed (default 0.1)\n";
	print "############################################################################################################################################################\n\n";
	my @now = localtime;print "\nINFO $now[2]:$now[1]:$now[0]\t";
	print "Starts SV-Quest\n";
	system("sleep 0.5s");
}

#-----------------------------------------------------------------------------------------------------------------------------------
sub return {
	open INPUT1, "<$fasta" or die "cant open sam file!\n";
	open INPUT2, "<$fasta";
	open (OUT, '>temp/temp0');

	my $line2 = <INPUT2>;
	while (my $line1 = <INPUT1>) {
		my $line2 = <INPUT2>;#1行先読みファイルの入力
		chomp($line1);#改行を除く
		print OUT "$line1\n" if($line1 =~ "\>");#先頭行を出力
		next if($line1 =~ "\>");
		print OUT "$line1";
		print OUT "\n"  if($line2 =~ "\>");#1行先読みファイルを識別に利用している
	}
}

#-----------------------------------------------------------------------------------------------------------------------------------
sub count {
	open INPUT3, "<temp/temp0";
	open (INDEX, '>temp/fasta.index');
	open (REFERENCE, '>temp/reference.fa');
	my $a = 1;
	while (my $line1 = <INPUT3>) {
	chomp($line1);
	$line1 =~ s/\>//;
	my $line2 = <INPUT3>;
	my $contig_size = length($line2);
	print INDEX "$a\t$line1\t$contig_size\n";
	print REFERENCE "\>$a\n$line2";
	$a++;
	}
close INPUT3;close INDEX;
}

#-----------------------------------------------------------------------------------------------------------------------------------
sub mapping {
	die "\nnFastq file is required !\n\n\n" if($fastq1 eq "");
	system("bwa index -p temp/reference.fa -a is temp/reference.fa 2> /dev/null");
	system("bwa aln -l $seed -k $seedpenalty -n $loose -i $loose -t $cpu temp/reference.fa $fastq1 2> /dev/null | bwa samse temp/reference.fa - $fastq1 2> /dev/null > temp/input.sam");
	
	next if($fastq2 eq "");
	system("bwa aln -l $seed -k $seedpenalty -n $loose -i $loose -t $cpu temp/reference.fa $fastq2 2> /dev/null | bwa samse temp/reference.fa - $fastq2 2> /dev/null |grep -v \"^@\" - >> temp/input.sam");
	system("sleep 3s");
}

#-----------------------------------------------------------------------------------------------------------------------------------
sub sam_convert {
	die "Mapped.bam file is required !\n\n\n" if($bam eq "");
	my @now = localtime;print "INFO $now[2]:$now[1]:$now[0]\tStarts samtools view\n" if($check == 1);
	system("samtools view -@ $cpu -h $bam 2> /dev/null > temp/input.sam") if($check == 1);
	my @now = localtime;print "INFO $now[2]:$now[1]:$now[0]\tDivides input.sam into F.sam and R.sam\n\n" if($check == 1);
}

#-----------------------------------------------------------------------------------------------------------------------------------
sub split {
	open INPUT, "<temp/input.sam" or die "cant open sam file!\n";
	open (OUTF, '>temp/F.sam');
	open (OUTR, '>temp/R.sam');
	open (NAME , '>temp/fasta');
	open (LOG , '>$outputdir/log.txt');
	my $totalp = 0;my $totalm = 0;my $unmap = 0;

	while (my $line1 = <INPUT>) {
		my $head = substr($line1,0,1);
		my @array = split(/\t/, $line1);
		if($head =~ /^\@/){
			print OUTF "$line1";
			print OUTR "$line1";
			my @name2 = "";
			my $head2 = substr($array[1],0,2);
			my @array2 = split(/\:/, $array[1]);
			print NAME "$array2[1]\t" if($head2 =~ /SN/);
		}
		next if($head =~ /^\@/);
		
		print OUTF "$line1" if($array[1] == 0);#forward reads
		print OUTF "$line1" if($array[1] == 256);#forward reads
		$totalp++ if($array[1] == 0);#forward reads
		$totalp++ if($array[1] == 256);#forward reads

		print OUTR "$line1" if($array[1] == 16);#reverse reads
		print OUTR "$line1" if($array[1] == 272);#reverse reads
		$totalm++ if($array[1] == 16);#reverse reads
		$totalm++ if($array[1] == 272);#reverse reads
	}
print LOG "Forward reads $totalp\nReverse reads $totalm\n\n";
close OUTF;close OUTR;close LOG;close NAME;close INPUT;
system("sleep 2s");
}#&split end

#-----------------------------------------------------------------------------------------------------------------------------------
sub extract {
	open CONTIG, "<temp/fasta.index" or die;
	my $counterA = 0;
	my $border = 0;
	while(my $cycle = <CONTIG>){#while1 start
		$counterA++;
		chomp($cycle);
		my @box =  split(/\t/,$cycle);
		#system("mkdir $box[1]");
		print "Chromosome name is $box[1]\n";
		my @now = localtime;print "\nINFO $now[2]:$now[1]:$now[0]\tStarts calculation of $box[1] mismatch\n";
		
		open F, "<temp/F_mpileup" or die ("File can\'t find.") && exit;
		open OUTF, ">temp/mismatch_ratio_F" or die "cant save \n"&&exit;
		my $count = 0;my @CoverageBoxF = "";
		while(my $line = <F>) {
			my @array = split(/\t/, $line);
			$CoverageBoxF[$array[1]] = $array[3];
			next unless($array[0] == $box[0]);
			my @seq = split(//, $array[4]);
			my $mis = 0; my $all = 0;my $nomiscov = 0;
			$count++;
			unless ($count == $array[1]) {
				while ($count<$array[1]) {
					print OUTF "$count\t0\t0\n";$count++;
				}
			}
			foreach (@seq) {#count mismatched bases
				$mis++ if ((/A/i) or (/T/i) or (/G/i) or (/C/i) or (/N/i));
				$all++ if ((/\./) or (/\,/) or (/A/i) or (/T/i) or (/G/i) or (/C/i) or (/N/i));
			}
			$CoverageBoxF[$array[1]] = $all - $mis;#coverageF
			my $mismuch_ratioF = 0;
			$mismuch_ratioF = $mis / $all unless($all == "0");
			$mismuch_ratioF = 0 if($mismuch_ratioF < $threshould);
			$mismuch_ratioF = 0 if($all <= $coverage2);
			$nomiscov = $all - $mis;#ver1.7追加部位
			print OUTF "$array[1]\t$nomiscov\t$mismuch_ratioF\n";#position, coverage, mismatch ratio of F #ver1.7修正部位
		}
		close F; close OUTF;


		open R, "<temp/R_mpileup" or die ("File can\'t find.") && exit;
		open OUTM, ">temp/mismatch_ratio_R" or die "cant save \n"&&exit;
		my $count = 0;my @CoverageBoxR = "";
		while(my $line = <R>) {
			my @array = split(/\t/, $line);
			$CoverageBoxR[$array[1]] = $array[3];#coverageR
			next unless($array[0] == $box[0]);
			my @seq = split(//, $array[4]);
			my $mis = 0; my $all = 0;my $nomiscov = 0;
			$count++;
			unless ($count == $array[1]) {
				while ($count<$array[1]) {
					print OUTM "$count\t0\t0\n";$count++;
				}
			}
			foreach (@seq) {#count mismatched bases
				$mis++ if ((/A/i) or (/T/i) or (/G/i) or (/C/i) or (/N/i));
				$all++ if ((/\./) or (/\,/) or (/A/i) or (/T/i) or (/G/i) or (/C/i) or (/N/i));
			}
		$CoverageBoxR[$array[1]] =  $all - $mis;
		my $mismuch_ratioR = 0;
		$mismuch_ratioR = $mis / $all unless($all == "0");
		$mismuch_ratioR = 0 if($mismuch_ratioR < $threshould);
		$mismuch_ratioR = 0 if($all <= $coverage2);
		$nomiscov = $all - $mis;#ver1.7追加部位
		print OUTM "$array[1]\t$nomiscov\t$mismuch_ratioR\n";#position, coverage, mismatch ratio of R #ver1.7修正部位
		}
		close R; close OUTM;
		print "mismatch ratio calculation end\n";
		
		my $half = int($size / 2);
		open F1, "<temp/mismatch_ratio_F" or die;
		open F2, "<temp/mismatch_ratio_F" or die;
		open FOUT1, ">temp/mismatch_value_F" or die;
		my @cov = (); my @misratioF = (); my @pos = ();my $i = 1;
	
		while (my $f1 = <F1>){
			chomp($f1);
			my @array = split(/\t/, $f1);
			$pos[$i] = $i;
			$cov[$i] = $array[1];
			$misratioF[$i] = $array[2];
			$i++;
		}
		my $ii = 1;
		while (my $f2 = <F2>){
			my @array2 = split(/\t/,$f2);
			my $right = $ii + $wide;
			my $left = $ii - $wide;
			$left = 1 if($left < 1);
			my @group = ();my $mismuch_valueF = 0;
			for (my $iii = $left; $iii <= $right;$iii++){#calclulate mismatch value in the loop
				$mismuch_valueF = $mismuch_valueF + $misratioF[$iii];#calclulate mismatch value of F
			}
		print FOUT1 "$pos[$ii]\t$array2[1]\t$mismuch_valueF\n";
		$ii++;
		}
		close F1;close F2;close FOUT1;

		open R1, "<temp/mismatch_ratio_R" or die;
		open R2, "<temp/mismatch_ratio_R" or die;
		open ROUT1, ">temp/mismatch_value_R" or die;
		my @cov = (); my @misratioR = (); my @pos = ();my $i = 1;
		while (my $r1 = <R1>){
			chomp($r1);
			my @array = split(/\t/, $r1);
			$pos[$i] = $i;
			$cov[$i] = $array[1];
			$misratioR[$i] = $array[2];
			$i++;
		}

		my $ii = 1;
		while (my $r2 = <R2>){
			my @array2 = split(/\t/,$r2);
			my $left = $ii - $wide;
			my $right = $ii + $wide;
			$left = 1 if($left < 1);
			my @group = (); my $mismuch_valueR = 0; my $region = 0;
			for (my $iii = $right; $iii >= $left;$iii--){
				$mismuch_valueR = $mismuch_valueR + $misratioR[$iii];#calclulate mismatch value of R
			}
		print ROUT1 "$pos[$ii]\t$array2[1]\t$mismuch_valueR\n";
		$ii++;
		}
		close R1;close R2;close ROUT1;
		print "mismatch value calculation end\n";
		
		
		open (LOG , '>>$outputdir/log.txt');
		open PLUSM, "temp/mismatch_value_F" or die;
		open MINUSM, "<temp/mismatch_value_R" or die;
		my $i = 1; my @plus_cov = "";#average read count of F
		my $p = <PLUSM>;
		while (my $p = <PLUSM>){
			chomp($p); my @boxP = split(/\t/,$p);
			$plus_cov[$i] = $boxP[1] if($i <= $box[2]);
			$i++;
		}
		
		my $i = 1; my @minus_cov = ""; my @cov = "";my $total = 0.1;#average read count of R
		my $m = <MINUSM>;
		while (my $m = <MINUSM>){
			chomp($m); my @boxM = split(/\t/,$m);
			$minus_cov[$i] = $boxM[1] if($i <= $box[2]);
			$cov[$i] = $plus_cov[$i] + $minus_cov[$i]  if($i <= $box[2]);
			$total = $total + $cov[$i];#total read count
			$i++;
		}
		
		
		$total = $total / $box[2];#average read count
		print LOG "average read count in $box[1] is $total\n";
		
		open F3, "<temp/mismatch_value_F" or die;
		open F4, "<temp/mismatch_value_F" or die;
		open FOUT3, ">temp/F_mismatch_mountain.txt" or die;
		open FOUT4, ">>temp/F_mismatch_mountain_all.txt" or die;
		
		my $counter = 0;my $sum = 0;my $start = 0;
		my $f4 = <F4>;
		my $mismatch_wideness = 0;
		while (my $f3 = <F3>){
			my $f4 = <F4>;
			my @array3 = split(/\t/,$f3);
			my @array4 = split(/\t/,$f4);
			chomp($array3[2]);chomp($array4[2]);
			my $CoverageThisPositionF = $CoverageBoxF[$array3[0]] / $total / 2;
			#next if($CoverageThisPositionF < $lowcoverage);
			if ($array3[2] > 0){
				$counter++;
				$sum += $array3[2];
				$start = $array3[0] if($counter == 1);
				$mismatch_wideness++;
			}
			if ($array3[2] > 0 && $array4[2] == 0){
				my $endness = $mismatch_wideness;
				$mismatch_wideness = 0;
				my $end = $array3[0];
				my $center = ($start + $end) / 2;
				#my $sizew = $end - $start;
				$center = int($center);
				my $slope = 0;
				$slope = $sum / $endness unless($endness == 0);
				print FOUT3 "$center\t$sum\t$endness\n" if($sum >= $call_threshould && $slope >= $slope_threshould);#version1.7修正部位
				print FOUT4 "$box[0]\t$center\t$sum\t$endness\n" if($sum >= $call_threshould && $slope >= $slope_threshould);#version1.7修正部位
				$counter = 0;$sum = 0;
			}
		}
		close F3;close F4;


		open R3, "<temp/mismatch_value_R" or die;
		open R4, "<temp/mismatch_value_R" or die;
		open ROUT3, ">temp/R_mismatch_mountain.txt" or die;
		open ROUT4, ">>temp/R_mismatch_mountain_all.txt" or die;

		my $counter = 0;my $sum = 0;my $start = 0;
		my $r4 = <R4>;
		my $mismatch_wideness = 0;
		while (my $r3 = <R3>){
			my $r4 = <R4>;
			my @array3 = split(/\t/,$r3);
			my @array4 = split(/\t/,$r4);
			chomp($array3[2]);chomp($array4[2]);
			my $CoverageThisPositionR = $CoverageBoxR[$array3[0]] / $total / 2;
			#next if($CoverageThisPositionR < $lowcoverage);
			if ($array3[2] > 0){
				$counter++;
				$sum += $array3[2];
				$start = $array3[0] if($counter == 1);
				$mismatch_wideness++;
			}
		
			if ($array3[2] > 0 && $array4[2] == 0){
				my $endness = $mismatch_wideness;
				$mismatch_wideness = 0;
				my $end = $array3[0];
				my $center = ($start + $end) / 2;
				$center = int($center);
				#my $sizew = $end - $start;
				my $slope = 0;
				$slope = $sum / $endness unless($endness == 0);
				print ROUT3 "$center\t$sum\t$endness\n" if($sum >= $call_threshould && $slope >= $slope_threshould);#version1.7修正部位
				print ROUT4 "$box[0]\t$center\t$sum\t$endness\n" if($sum >= $call_threshould && $slope >= $slope_threshould);#version1.7修正部位
				$counter = 0;$sum = 0;
			}
		}
		close R3;close R4;
		system("sleep 5s");
		
		my @now = localtime;print "INFO $now[2]:$now[1]:$now[0]\tMismatch calculation of $box[1] was end. Now starts indel call of $box[1] from mismatch information\n";
		#call Insertion and small deletion candidate
		open PLUS2, "<temp/F_mismatch_mountain.txt" or die;
		open MINUS2, "<temp/R_mismatch_mountain.txt" or die;
		open OUTPUT3I, ">temp/In.txt" or die;
		open OUTPUT3IC, ">>$output" or die;
		open TARGETSITE, ">>$outputdir/target.txt" or die;

		print OUTPUT3I "\#chrom\tstart\tend\tmismatch\n";
		print OUTPUT3IC "\#chrom\ttype\tpositon\tmismatch\n" if($counterA eq 1);
		my @mismuchMP =();my @mismuchMM =();my $i = 0;
		while (my $minus2 = <MINUS2>){
			my @array2 = split(/\t/,$minus2);chomp($array2[1]);
			$mismuchMP[$i] = $array2[0];
			$mismuchMM[$i] = $array2[1];
			$i++;
		}
		my $length = @mismuchMP;
		
		while (my $plus2 = <PLUS2>){
			
			my @array1 = split(/\t/,$plus2); chomp($array1[1]);
			for (my $i = 0; $i < $length; $i++){
			
				my $right =  $mismuchMP[$i] + $Lborder;#≤ 300 bp
				my $left =  $mismuchMP[$i] -  $Lborder;#≤ 300 bp
				$array1[0] =~ s/\s+//g;
				
				if ($array1[0] >= $left && $array1[0] <= $right){#Both Pf and Pr are found (≤ 300 bp)
					my $value = ($array1[1] + $mismuchMM[$i]) / 2;
					$value  = int($value);
					
					if($array1[0] <= $mismuchMP[$i]){
						#intra-chromosomal mappingを取り出すために、ターゲット領域の±250bpの長いバイト配列も取り出す
						my $rightseq = $mismuchMP[$i] + 100;
						my $rightseqex = $mismuchMP[$i] + 250;
						my $leftseq = $array1[0] - 100;
						my $leftseqex = $array1[0] - 250;
						system("samtools faidx $fasta $box[1]:$leftseqex-$rightseqex > temp/targetextension1.fa");
						system("samtools faidx $fasta $box[1]:$leftseq-$rightseq > temp/target1.fa");
						system("cat temp/targetextension1.fa $Tn > temp/mixed1.fasta");
						
						for (my $num = 1; $num <= $loop;$num++){#iterative local assembly
							my $numminus = $num - 1;
							my $numplus = $num + 1;
							system("mirabait -b temp/target1.fa -p temp/F.fastq temp/R.fastq 1> /dev/null") if($num eq 1);
							system("mirabait -b temp/contigs_select.fa -p temp/F.fastq temp/R.fastq 1> /dev/null") if($num > 1);
							last unless (-s "bait_match_F.fastq");
							system("minimap2 -ax sr temp/mixed${num}.fasta bait_match_F.fastq bait_match_R.fastq 2> /dev/null |samtools sort - 2> /dev/null | samtools view -F 14 -h - | awk \'\$7 \!\~ \/=\/\' | samtools fastq > temp/intrachromosomal${num}.fq");
							system("rm bait_match_F.fastq bait_match_R.fastq");
							system("cat temp/intrachromosomal*.fq > temp/${num}_merged.fq") if($num > 1);
							system("seqkit rmdup -n temp/${num}_merged.fq > temp/${num}_removedup.fq 2> temp/seqkitlog") if($num > 1);
							system("spades.py -s temp/intrachromosomal${num}.fq --cov-cutoff 5 -o temp/spades${num} -k 21,31,41 -t 12 --only-assembler 1> temp/log") if($num eq 1);
							system("spades.py -s temp/${num}_removedup.fq --cov-cutoff 5 -o temp/spades${num} -k 21,31,41 -t 12 --only-assembler 1> temp/log") if($num > 1);
							system("seqkit seq -m 150 temp/spades${num}/contigs.fasta > temp/contigs_select.fa");
							system("rm temp/${num}_merged.fq temp/${num}_removedup.fq 2> /dev/null") if($num > 1);
							system("rm temp/contigs_select.fa") unless (-s "temp/contigs_select.fa");
							last unless (-s "temp/contigs_select.fa");
							system("cat temp/contigs_select.fa $Tn > temp/mixed${numplus}.fasta") unless ($num eq $loop);
							system("rm -rf temp/spades${num}/K* temp/spades${num}/*paths temp/spades${num}/input_dataset.yaml temp/spades${num}/dataset.info temp/spades${num}/before_rr.fasta temp/spades${num}/warnings.log temp/spades${num}/assembly_graph_with_scaffolds.gfa temp/spades${num}/params.txt temp/spades${num}/tmp temp/spades${num}/misc temp/spades${num}/mismatch_corrector temp/spades${num}/corrected temp/seqkitlog");
							system("mv temp/spades${num} temp/spades${num}-$array1[0]-$mismuchMP[$i]");
						}
						
						system("blastn -db temp/blastDB/targetDB -query temp/contigs_select.fa -evalue 1e-100 -outfmt 7 -out temp/blast.txt");
						system("cat temp/blast.txt | awk \'/hits found/{getline;print}\' | grep -v \"#\" >> temp/blast-best.txt");
						
						open DB, "<temp/blast-best.txt" or die;
						my $blastdb = <DB>;
						my @blastarray = split(/\t/,$blastdb);
						#identity ≥95%を満たすヒット
						print OUTPUT3I "$samplename\t$box[1]\t$array1[0]\t$mismuchMP[$i]\t$value\t$blastarray[1]\t$blastarray[2]\n" if($blastarray[2] >= 95);
						print OUTPUT3IC "$samplename\t$box[1]\tType_I_SV\t$array1[0]\-$mismuchMP[$i]\t$value\t$blastarray[1]\t$blastarray[2]\n" if($blastarray[2] >= 95);
						print TARGETSITE "$samplename\t$box[1]\t$array1[0]\t$mismuchMP[$i]\t$value\t$blastarray[1]\t$blastarray[2]\n" if($blastarray[2] >= 95);
						#identity ≥95%を満たさない
						print OUTPUT3I "$samplename\t$box[1]\t$array1[0]\t$mismuchMP[$i]\t$value\n" if($blastarray[2] < 95);
						print OUTPUT3IC "$samplename\t$box[1]\tType_I_SV\t$array1[0]\-$mismuchMP[$i]\t$value\n" if($blastarray[2] < 95);
						
						system("mv temp/blast-best.txt temp/blastDB/$array1[0]\-$mismuchMP[$i]_blast-best.txt");
						system("rm temp/contigs_select.fa temp/blast.txt");
						
					}elsif($array1[0] > $mismuchMP[$i]){
						#intra-chromosomal mappingを取り出すために、ターゲット領域の±250bpの長いバイト配列も取り出す
						my $rightseq = $array1[0] + 100;
						my $rightseqex = $array1[0] + 250;
						my $leftseq = $mismuchMP[$i] - 100;
						my $leftseqex = $mismuchMP[$i] - 250;
						system("samtools faidx $fasta $box[1]:$leftseqex-$rightseqex > temp/targetextension1.fa");
						system("samtools faidx $fasta $box[1]:$leftseq-$rightseq > temp/target1.fa");
						system("cat temp/targetextension1.fa $Tn > temp/mixed1.fasta");
						
						for (my $num = 1; $num <= $loop;$num++){#iterative local assembly
							my $numminus = $num - 1;
							my $numplus = $num + 1;
							system("mirabait -b temp/target1.fa -p temp/F.fastq temp/R.fastq 1> /dev/null") if($num eq 1);
							system("mirabait -b temp/contigs_select.fa -p temp/F.fastq temp/R.fastq 1> /dev/null") if($num > 1);
							last unless (-s "bait_match_F.fastq");
							system("minimap2 -ax sr temp/mixed${num}.fasta bait_match_F.fastq bait_match_R.fastq 2> /dev/null |samtools sort - 2> /dev/null | samtools view -F 14 -h - | awk \'\$7 \!\~ \/=\/\' | samtools fastq > temp/intrachromosomal${num}.fq");
							system("rm bait_match_F.fastq bait_match_R.fastq");
							system("cat temp/intrachromosomal*.fq > temp/${num}_merged.fq") if($num > 1);
							system("seqkit rmdup -n temp/${num}_merged.fq > temp/${num}_removedup.fq 2> temp/seqkitlog") if($num > 1);
							system("spades.py -s temp/intrachromosomal${num}.fq --cov-cutoff 5 -o temp/spades${num} -k 21,31,41 -t 12 1> temp/log") if($num eq 1);
							system("spades.py -s temp/${num}_removedup.fq --cov-cutoff 5 -o temp/spades${num} -k 21,31,41 -t 12 1> temp/log") if($num > 1);
							system("seqkit seq -m 150 temp/spades${num}/contigs.fasta > temp/contigs_select.fa");
							system("rm temp/${num}_merged.fq temp/${num}_removedup.fq 2> /dev/null") if($num > 1);
							system("rm temp/contigs_select.fa") unless (-s "temp/contigs_select.fa");
							last unless (-s "temp/contigs_select.fa");
							system("cat temp/contigs_select.fa $Tn > temp/mixed${numplus}.fasta") unless ($num eq $loop);
							system("rm -rf temp/spades${num}/K* temp/spades${num}/*paths temp/spades${num}/input_dataset.yaml temp/spades${num}/dataset.info temp/spades${num}/before_rr.fasta temp/spades${num}/warnings.log temp/spades${num}/assembly_graph_with_scaffolds.gfa temp/spades${num}/params.txt temp/spades${num}/tmp temp/spades${num}/misc temp/spades${num}/mismatch_corrector temp/spades${num}/corrected temp/seqkitlog");
							system("mv temp/spades${num} temp/spades${num}-$array1[0]-$mismuchMP[$i]");
						}
						
						system("blastn -db temp/blastDB/targetDB -query temp/contigs_select.fa -evalue 1e-10 -outfmt 7 -out temp/blast.txt");
						system("cat temp/blast.txt | awk \'/hits found/{getline;print}\' | grep -v \"#\" >> temp/blast-best.txt");
						
						open DB, "<temp/blast-best.txt" or die;
						my $blastdb = <DB>;
						my @blastarray = split(/\t/,$blastdb);
						#identity ≥95%を満たすヒット
						print OUTPUT3I "$samplename\t$box[1]\t$mismuchMP[$i]\t$array1[0]\t$value\t$blastarray[1]\t$blastarray[2]\n" if($blastarray[2] >= 95);
						print OUTPUT3IC "$samplename\t$box[1]\tType_I_SV\t$mismuchMP[$i]-$array1[0]\t$value\t$blastarray[1]\t$blastarray[2]\n" if($blastarray[2] >= 95);
						print TARGETSITE "$samplename\t$box[1]\t$mismuchMP[$i]\t$array1[0]\t$value\t$blastarray[1]\t$blastarray[2]\n" if($blastarray[2] >= 95);
						#identity ≥95%を満たさない
						print OUTPUT3I "$samplename\t$box[1]\t$mismuchMP[$i]\t$array1[0]\t$value\n" if ($blastarray[2] < 95);
						print OUTPUT3IC "$samplename\t$box[1]\tType_I_SV\t$mismuchMP[$i]-$array1[0]\t$value\n" if ($blastarray[2] < 95);
						
						system("mv temp/blast-best.txt temp/blastDB/$array1[0]\-$mismuchMP[$i]_blast-best.txt");
						system("rm temp/contigs_select.fa temp/blast.txt");
						
					}
				}
			}
		}
		
		#call large deletion candidate
		open COVF, "<temp/mismatch_ratio_F" or die;
		open COVR, "<temp/mismatch_ratio_R" or die;
		my @coverageF = "";	
		while(my $f = <COVF>){
		my @fbox = split(/\t/,$f);
		$coverageF[$fbox[0]] = $fbox[1];
		}
		my @coverageR = "";	
		while(my $r = <COVR>){
			my @rbox = split(/\t/,$r);
			$coverageR[$rbox[0]] = $rbox[1];
		}

		open PLUS2, "<temp/F_mismatch_mountain.txt" or die;
		open OUTPUT4D, ">temp/Del.txt" or die;
		open OUTPUT4DC, ">>$output" or die;
		print OUTPUT4D "\#chrom\tstart\tend\trelative_coverage\n";
			while (my $checkp = <PLUS2>){
			chomp($checkp);
			my @inputp = split(/\t/,$checkp);
			my $right =  $inputp[0] + $max;
			my $left = $inputp[0] + $Lborder;##Pf (read length ≥ , ≤ 100 kbp)
			
			for (my $i = 0; $i < $length; $i++){
				if ($left <= $mismuchMP[$i] && $mismuchMP[$i] <= $right){

					my $sumC = 0;my $ii = 0;
					for (my $k = $inputp[0]; $k <= $mismuchMP[$i]; $k++){
						$sumC = $sumC + $coverageF[$k] + $coverageR[$k];
						$ii++;
					}
					my $AveC = $sumC / $ii;
					my $read_ratio = 0;
					$read_ratio = $AveC / $total;
					print OUTPUT4D "$box[1]\t$inputp[0]\t$mismuchMP[$i]\t$read_ratio\n" if($read_ratio < $coverage);
					print OUTPUT4DC "$box[1]\tType_II_SV\t$inputp[0]\-$mismuchMP[$i]\t$read_ratio\n" if($read_ratio < $coverage);
				}
			}
		}
		my @now = localtime;print "INFO $now[2]:$now[1]:$now[0]\t$box[1] end.\n";
		system("sleep 2s");
		if($duplicateRemove == 1){
			open LAST1, "<temp/In.txt" or die;
			open LAST2, "<temp/In.txt" or die;
			open OUTP, ">temp/In_non-redundant.txt" or die;
			open OUTPP, ">>.non-redundant.txt" or die;
			print OUTP "\#chrom\tstart\tend\tmismatch\n";
			print OUTPP "\#chrom\ttype\positon\tmismatch\n";
			my $input1 = <LAST1>;
			my $input2 = <LAST2>;my $input2 = <LAST2>;
			while($input1 = <LAST1>){
				$input2 = <LAST2>;
				next if($input1 eq "");
				chomp($input1);
				chomp($input2);
				my @lbox1 = split(/\t/,$input1);
				my @lbox2 = split(/\t/,$input2);
				if(($lbox1[1] + 300) > $lbox2[1]){
					my $firstP = $lbox1[1];
					my $endP = $lbox2[2];
					$lbox1[3] = $lbox2[3] if($lbox1[3] < $lbox2[3]);
					for(my $e=1; $e < 1000; $e++){
						$input1 = <LAST1>;
						$input2 = <LAST2>;
						chomp($input1);
						chomp($input2);
						my @lbox3 = split(/\t/,$input1);
						my @lbox4 = split(/\t/,$input2);
						$lbox1[3] = $lbox4[3] if(($lbox1[1] + 300) > $lbox4[1] && $lbox1[3] < $lbox4[3]);
						$endP = $lbox4[2] if(($lbox1[1] + 300) > $lbox4[1]);
						print OUTP "$lbox1[0]\t$firstP\t$lbox1[2]\t$lbox1[3]\n" if((($lbox1[1] + 300) < $lbox4[1]) or ($lbox4[1] eq ""));
						print OUTPP "$lbox1[0]\tType_I_SV\t$firstP\-$lbox1[2]\t$lbox1[3]\n" if((($lbox1[1] + 300) < $lbox4[1]) or ($lbox4[1] eq ""));
						last if((($lbox1[1] + 300) < $lbox4[1]) or ($lbox4[1] eq ""));
					}#for
				}elsif((($lbox1[1] + 300) < $lbox2[1]) or ($lbox2[1] eq "")){
					print OUTP "$input1\n";
					print OUTPP "$lbox1[0]\tType_I_SV\t$lbox1[1]\-$lbox1[2]\t$lbox1[3]\n";
				}
			}
		}
		
		if($duplicateRemove == 1){
			open LAST3, "<temp/Del.txt" or die;
			open LAST4, "<temp/Del.txt" or die;
			open OUTD, ">temp/Del_non-redundant.txt" or die;
			open OUTDR, ">>.non-redundant.txt" or die;
			print OUTD "\#chrom\tstart\tend\trelative_coverage\n";
			my $input1 = <LAST3>;
			my $input2 = <LAST4>;my $input2 = <LAST4>;
			while($input1 = <LAST3>){
				$input2 = <LAST4>;
				next if($input1 eq "");
				chomp($input1);
				chomp($input2);
				my @lbox1 = split(/\t/,$input1);
				my @lbox2 = split(/\t/,$input2);
				if($lbox1[1] == $lbox2[1] || $lbox1[2] == $lbox2[2]){
					
					my $firstP = $lbox1[1];
					my $endP = $lbox2[2];
					$lbox1[3] = $lbox2[3] if($lbox1[3] < $lbox2[3]);
					for(my $e=1; $e < 1000; $e++){
						$input1 = <LAST3>;
						$input2 = <LAST4>;
						chomp($input1);
						chomp($input2);
						my @lbox3 = split(/\t/,$input1);
						my @lbox4 = split(/\t/,$input2);
						$lbox1[3] = $lbox4[3] if(($lbox1[1] + 300) > $lbox4[1] && $lbox1[3] < $lbox4[3]);
						$endP = $lbox4[2] if(($lbox1[1] + 300) > $lbox4[1]);
						print OUTD "$lbox1[0]\t$firstP\t$lbox1[2]\t$lbox1[3]\n" if((($lbox1[1] != $lbox4[1]) && ($lbox1[2] != $lbox4[2]))  or ($lbox4[1] eq ""));
						print OUTDR "$lbox1[0]\tType_II_SV\t$firstP\-$lbox1[2]\t$lbox1[3]\n" if((($lbox1[1] != $lbox4[1]) && ($lbox1[2] != $lbox4[2])) or ($lbox4[1] eq ""));
						last if((($lbox1[1] != $lbox4[1]) && ($lbox1[2] != $lbox4[2])) or ($lbox4[1] eq ""));
					}#for
				}else{
					print OUTD "$input1\n";
					print OUTDR "$lbox1[0]\tType_II_SV\t$lbox1[1]\-$lbox1[2]\t$lbox1[3]\n";
				}
			}
		}
	}#while1end
}#&extract end


####################################################################################################################################
