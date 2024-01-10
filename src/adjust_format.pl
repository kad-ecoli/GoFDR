#!/usr/bin/perl -w
$file = pop @ARGV;

open A,"tmp_data/max_id/$file";
while(<A>){
	chomp;
	($g,$max_id) = split /\t/;
	$max{$g} = $max_id;
}
close A;

for $type("BP","MF","CC"){
	open A,"tmp_data/results_psibls_$type/$file";
	open O,">tmp_data/format_results_psibls_$type/$file";
	while(<A>){
		chomp;
		if(/(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)/){
			$g = $1;
			$r = $3;
			$f = $4;
			$m = $max{$g};
			printf O "%s\t%s\t%.3f\t%.3f\t%.3f\n",$g,$file,$r,$f,$m;
		}
	}
	close A;
	close O;
}
