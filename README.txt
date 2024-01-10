#################################
#                               #
#             GoFDR             # 
#                               #
#################################
#                               #
#        Authors:               #
#                               #
#          Qingtian Gong        #
#          Wei      Ning        #
#          Weidong  Tian*       #
#                               #
#                               #
#     *Corresponding Author     #
#            E-mail:            #
#   weidong.tian@fudan.edu.cn   #
#                               #
#################################

Usage:

1. Decompress function annotation by 
   "bunzip2 database/annotations/noiea_confirmed_iea_annotations.txt.bz2"

2.Set UniRef90 protein sequence database in "database/uniprot_sequence/uniref90.fasta"; (Because the size of the latest sequence database is very huge, we didn't include it in this source code package. You may download it from http://www.uniprot.org/downloads) 

3.Input sequence of query protein in "in_data/sequences" directory, named with its ID such as "Query1";

4.List all IDs of query proteins in "in_data/query_list" file;

5.Run "perl src/run.pl"

6.Check "out_data/[MF,BP,CC]/#ID#_final.txt" for prediction result.


################################

Contents:

---./
  |
  |---database
  |  |
  |  |---annotations #include noIEA+confirmed IEA annotations, GO specificity, GO categories, child-ancestor relationship
  |  |
  |  |---uniprot_sequence #NEED SET#for customer to set sequence database for PSI-BLAST by self
  |  |
  |  |---prob_tables #raw_score2probablity conversion tables generated from training sequences
  |
  |---in_data
  |  |
  |  |---sequences #NEED SET#store each query protein sequence in FASTA format with its ID as filename
  |  |
  |  |---query_list #NEED SET#a list of query proteins for prediction
  |  |
  |  |---one_line_sequences #intermediate dir
  |
  |---tmp_data #intermediate dir
  |
  |---out_data
  |  |
  |  |---BP #Biological Process
  |  |
  |  |---MF #Molecular Function
  |  |
  |  |---CC #Cellular Component
  |     |
  |     |-Query1_final.txt #IMPORTANT# final prediction result
  |
  |---src #source code dir


################################

Compiler Version:

gcc,g++  (GCC 5.1.0)
perl     (Perl v5.16.3)
