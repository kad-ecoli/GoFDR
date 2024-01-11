#!/usr/bin/perl -w
use strict;
use File::Basename;
use Cwd 'abs_path';
my $src=dirname(abs_path(__FILE__));
my $rootdir=dirname("$src");
chdir($rootdir);

&mkdir("out_data/MF");
&mkdir("out_data/BP");
&mkdir("out_data/CC");
&mkdir("in_data/one_line_sequences");
&mkdir("tmp_data/alignments_psibls");
&mkdir("tmp_data/assigned_GO_terms_psibls");
&mkdir("tmp_data/format_results_psibls_BP");
&mkdir("tmp_data/format_results_psibls_CC");
&mkdir("tmp_data/format_results_psibls_MF");
&mkdir("tmp_data/low_assigned_GO_terms_psibls");
&mkdir("tmp_data/low_psibls_aln");
&mkdir("tmp_data/max_id");
&mkdir("tmp_data/psibls_aln");
&mkdir("tmp_data/results_psibls");
&mkdir("tmp_data/results_psibls_BP");
&mkdir("tmp_data/results_psibls_CC");
&mkdir("tmp_data/results_psibls_MF");
&mkdir("tmp_data/up_go_list");


my $uniprot_db = "uniref90.fasta";

system qq(src/read_list_psibls in_data/query_list $uniprot_db);
system qq(src/get_assigned_GO_terms_psibls in_data/query_list);
system qq(perl src/one_line.pl Query1);
system qq(perl src/filter_60_get_maxid.pl Query1);
system qq(src/func_pred_psibls Query1);
system qq(src/divide_categories_psibls);
system qq(perl src/adjust_format.pl Query1);
system qq(src/produce_prob MF Query1);
system qq(src/produce_prob CC Query1);
system qq(src/produce_prob BP Query1);
system qq(src/correct_parent_child_psibls MF Query1);
system qq(src/correct_parent_child_psibls CC Query1);
system qq(src/correct_parent_child_psibls BP Query1);
system qq(perl src/combine.pl Query1);

exit();

sub mkdir
{
    my ($dir)=@_;
    system("mkdir -p $dir") if (!-d "$dir");
}
