CC=g++
CFLAGS=-static
PROGRAM=correct_parent_child_psibls divide_categories_psibls func_pred_psibls get_assigned_GO_terms_psibls get_blast_aln_psibls read_list_psibls produce_prob

all: ${PROGRAM}

correct_parent_child_psibls: correct_parent_child_psibls.cpp
	${CC} ${CFLAGS} correct_parent_child_psibls.cpp -o correct_parent_child_psibls

divide_categories_psibls: divide_categories_psibls.cpp
	${CC} ${CFLAGS} divide_categories_psibls.cpp -o divide_categories_psibls

func_pred_psibls: func_pred_psibls.cpp
	${CC} ${CFLAGS} func_pred_psibls.cpp -o func_pred_psibls

get_assigned_GO_terms_psibls: get_assigned_GO_terms_psibls.cpp
	${CC} ${CFLAGS} get_assigned_GO_terms_psibls.cpp -o get_assigned_GO_terms_psibls

get_blast_aln_psibls: get_blast_aln_psibls.cpp
	${CC} ${CFLAGS} get_blast_aln_psibls.cpp -o get_blast_aln_psibls

read_list_psibls: read_list_psibls.cpp
	${CC} ${CFLAGS} read_list_psibls.cpp -o read_list_psibls

produce_prob: produce_prob.c
	${CC} ${CFLAGS}	produce_prob.c -o produce_prob

clean:
	rm ${PROGRAM}
