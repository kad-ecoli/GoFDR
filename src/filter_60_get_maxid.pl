#!/usr/bin/perl -w

$file = pop @ARGV;
open B,"tmp_data/psibls_aln/$file";
open O1,">tmp_data/low_psibls_aln/$file";
open O2,">tmp_data/up_go_list/$file";
$tag = 0;
while(<B>){
	chomp;
	if(/^>(\S+)\s(\S+)\s(\S+)\s(\S+)/){
		$p = $1;
		$map = $3;
		$all = $4;
		$ratio = $map/$all;
		if($ratio <= 0.6){
			$id{$p} = $ratio;
			$tag = 1;
			print O1 $_."\n";
		}else{
			$up{$p}=1;
			$tag = 0;
		}
	}else{
		if($tag == 1){
			print O1 $_."\n";
		}
	}
}
close B;
close O1;

open A,"tmp_data/assigned_GO_terms_psibls/$file";
open O4,">tmp_data/low_assigned_GO_terms_psibls/$file";
while(<A>){
	chomp;
	if(/(\S+)\s(\S+)/){
		$p = $1;
		$g = $2;
		if($up{$p}){
			$up_go{$g} = 1;
		}
		if($id{$p}){
			$g2p{$g}{$p}=$id{$p};
			$x{$_}=1;
		}
	}
}
close A;

for $k(sort keys %x){
	print O4 $k."\n";
}
close O4;

for $k(sort keys %up_go){
	print O2 $k."\n";
}
close O2;

open O3,">tmp_data/max_id/$file";
for $k(sort keys %g2p){
	$i = 0;
	for $j(sort {$g2p{$k}{$b}<=>$g2p{$k}{$a}} keys %{$g2p{$k}}){
		last if($i == 1);
		printf O3 "%s\t%.4f\n",$k,$g2p{$k}{$j};
		$i++;
	}
}
close O3;
