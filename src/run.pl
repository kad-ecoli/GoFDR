#!/usr/bin/perl -w
$uniprot_db = "uniref90.fasta";

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
