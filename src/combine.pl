#!/usr/bin/perl -w

$query = pop @ARGV;
open A,"database/annotations/categories";
while(<A>){
	if(/(\S+)\s(\S+)/){
		$type{$1} = $2;
	}
}
close A;

for $cate("BP","MF","CC"){
	open B,"out_data/$cate/$query\_updated.txt";
	while(<B>){
		chomp;
		if(/(\S+)\s(\S+)\s(\S+)/){
			$p = $1;
			$g = $2;
			$s = $3;
		}
		$update{$cate}{$p}{$g}=$s;
	}
	close B;
}
open B,"tmp_data/up_go_list/$query";
while(<B>){
	chomp;
	$cate = $type{$_};
	$update{$cate}{$query}{$_}=1;
}
close B;

for $cate("BP","MF","CC"){
	open O,">out_data/$cate/$query\_final.txt";
	for $k(keys %{$update{$cate}{$query}}){
		printf O "%s\t%s\t%.3f\n",$query,$k,$update{$cate}{$query}{$k};
	}
	close O;
}
