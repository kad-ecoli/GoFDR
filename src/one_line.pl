#!/usr/bin/perl -w

$q = pop @ARGV;

open A,"in_data/sequences/$q";
open O,">in_data/one_line_sequences/$q";
while(<A>){
	chomp;
	if(/^>/){
		print O $_."\n";
	}else{
		print O $_;
	}
}
print O "\n";
close A;
close O;


