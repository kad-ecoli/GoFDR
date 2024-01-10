//usage src/produce_prob MF Query1
#include <stdio.h> 
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define DELIMITER "\t\n"
#define INPUT_DIR ""
#define MAX_ATTR_NAME 400
#define MAX_ATTR_BRIEF 400
#define MAX_ATTR_DESCRIP 400
#define MAX_LINE_LENGTH 2000

typedef struct _Cutoffs *CutoffsPtr;    
typedef struct _Cutoffs {
	double left;
	double right;
	double *score_range; 
	double *prob_table; 
}CutoffsStr;

typedef struct _ScoreTables *ScoreTablesPtr;    
typedef struct _ScoreTables {
	int num_scores;
	double *score_range; 
	double *prob_table; 
}ScoreTablesStr;

typedef struct _PredScore *PredScorePtr;    
typedef struct _PredScore {
	char status; 
	double score; 
	double prob;
	double go_freq;
	double go_spec; 
	char marker; 
	int go_id; 
	int gene_id; 
	char freq_th;
	char spec_th;
	double evalue;
	double max_id;
	double avr_id; 
	char go[400];
	char gene[400]; 
} PredScoreStr;

typedef struct _Vertice *VerticePtr;              // pointer to functional attribute structure 
typedef struct _Vertice {
	int num_trues;
	int num_candidate_gos; 
	int *candidate_gos; 
	char marker; 
	PredScorePtr *pred;
	char name[400]; 
}
VerticeStr;

VerticePtr Training_genes; 

typedef struct _VertexAttr *VertexAttrPtr;              // pointer to functional attribute structure 
typedef struct _VertexAttr {
	int num_trues;
	int num_homologs; 
	int *homologs; 
	double go_spec; 
	PredScorePtr *pred;
	char name[400]; 
}
VertexAttrStr;
VertexAttrPtr GO_attrs;

int wc (char any_file[]); 
int GetScores(char attr_file[], PredScorePtr *predictions);
double compute_auc (int num_data_points, PredScorePtr *predictions, int num_true_points, char *out_file, char roc_option, char write_option); 
int cmp( const void *a ,const void *b); 
int reverse_cmp( const void *a ,const void *b); 
int get_max_ID(char any_file[]); 
int ReadAttrs(int num_go_attrs, char attr_file[]); 
int ReadGenes(int num_genes, char attr_file[]); 
double ReadTable (char *score_file, int avr_id_range, int freq_range, double ***prob_table);
void tostring(char str[], int num); 
void ftoa(float n, char *res, int afterpoint);
void reverse(char *str, int len);
int intToStr(int x, char str[], int d);

int main(int argc, char *argv[]) {
	char raw_score_file[400] = INPUT_DIR;
	char go_attr_file[400] = INPUT_DIR;
	char gene_attr_file[400] = INPUT_DIR;
	char output_table[400] = INPUT_DIR; 
	char score_table[400] = INPUT_DIR; 
	char output_file[400] = INPUT_DIR; 
	char blast_option[400] = INPUT_DIR; 
	char GO_branch[400];
	char query_id[400];
	int i, j, l, m, k, num_data_points,  TP, FP, TN, FN; 
	int num_true_points; 
	double area, prec, recall, fp, tp, old_fp, old_tp, tmp_tp, roc_auc_01; 
	double roc_auc, roc_auc_2; 
	char *target_attrs_flag; 
	FILE *OUT; 
	double score, z, max_auc, max_z; 
	char roc_option, write_option; 
	double avr_t, avr_f; 
	int NT, NF; 
	char abbr[400]; 
	int num_go_attrs, num_genes; 
	PredScorePtr *predictions, pred_s; 
	PredScorePtr *tmp_predictions, tmp_pred; 
	double *score_range; 
	double *freq_range; 
	double *spec_range; 
	double freq_th_left, spec_th_left; 
	double freq_th_right, spec_th_right; 
	char tmp[400]; 
	char fold; 
	char table_abbr[400] = INPUT_DIR; 
	double ***prob_table; 
	CutoffsPtr spec_cutoffs, freq_cutoffs, prob_cutoffs; 
	ScoreTablesPtr **adjustment_table; 
	char option; 
	
	strcpy(GO_branch, argv[1]);
	strcpy(query_id, argv[2]);

	strcat(raw_score_file, "tmp_data/format_results_psibls_"); 
	strcat(raw_score_file, GO_branch); 
	strcat(raw_score_file, "/"); 
	strcat(raw_score_file, query_id); 

	strcat(blast_option, "database/prob_tables/"); 
	
	strcat(output_file, "out_data/"); 
	strcat(output_file, GO_branch); 
	strcat(output_file, "/"); 
	strcat(output_file, query_id); 
	strcat(output_file, "_adjusted.txt"); 

	printf ("%s\n%s\n", raw_score_file, output_file); 

	score_range = (double *) malloc (100 * sizeof(double)) - 1; 
	freq_range = (double *) malloc (10 * sizeof(double)) - 1; 
	spec_range = (double *) malloc (10 * sizeof(double)) - 1; 
	for (i = 1; i <= 13; i ++) { 
		score_range[i] = 3.0 - 0.5 * (i - 1); 
	}
	score_range[14] = -9; 

	freq_range[1] = 0; 
	freq_range[2] = 0.05; 
	freq_range[3] = 0.1; 
	freq_range[4] = 0.2; 
	freq_range[5] = 0.3; 
	freq_range[6] = 0.5; 
	freq_range[7] = 0.7; 
	freq_range[8] = 1.0;

	for (i = 1; i <= 6; i ++) { 
		spec_range[i] = 0 + 0.1 * i; 
	}

	prob_table = (double ***) calloc(10, sizeof(double **)) - 1; 
	for (i = 1; i <= 10; i ++) { 
		prob_table[i] = (double **) calloc(10, sizeof(double *)) - 1; 
		for (j = 1; j <= 10; j ++) { 
			prob_table[i][j] = (double *) calloc (100, sizeof(double)) - 1;
		}
	}

	num_data_points = wc(raw_score_file); 
	printf ("num_data_points %d\n", num_data_points); 
	predictions = (PredScorePtr *) malloc(num_data_points * sizeof (PredScorePtr)); 
	for (i = 0; i < num_data_points; i ++) { 
		predictions[i] = (PredScorePtr) calloc(1, sizeof (PredScoreStr));
	}
	num_true_points = GetScores(raw_score_file, predictions); 
	printf ("%d\n", num_true_points); 
	printf ("%d\n", num_data_points); 
	NT = 0; 
	for (i = 1; i <= 6; i ++) { 
		if (i == 1) { 
			spec_th_left = -20; 
		} else {  
			spec_th_left = spec_range[i-1]; 
		}
		spec_th_right = spec_range[i];
		for (j = 1; j <= 7; j ++) { 
			freq_th_left = freq_range[j]; 
			if (j == 1) { 
				freq_th_left = -20;
			}
			if (j == 7) { 
				freq_th_right = 1; 
			} else { 
				freq_th_right = freq_range[j + 1]; 
			}
			m = 0; 
			for (l = 0; l < num_data_points; l ++) { 
				tmp_pred = predictions[l]; 
				if (tmp_pred->max_id > spec_th_left && tmp_pred->max_id <= spec_th_right) { 
					if (tmp_pred->go_freq > freq_th_left && tmp_pred->go_freq <= freq_th_right) { 
						tmp_pred->spec_th = i;
						tmp_pred->freq_th = j; 
						m ++; 
						NT ++; 
					}
				}
			}
			printf ("max_id %f-%f, freq %f-%f, size %d i = %d j = %d\n", spec_th_left, spec_th_right, freq_th_left, freq_th_right, m, i,j); 
		}
	}
	printf ("NT %d %d \n", NT, num_data_points); 


	NT = 0; 

	if ((OUT= fopen(output_file, "w")) == NULL) {
		printf ("cannot open file to write\n"); 
		return 0;
	}

	for (i = 1; i <= 6; i ++) { 
		for (j = 1; j <= 7; j ++) { 
			strcpy(score_table, blast_option); 
			strcat (score_table, GO_branch); 
			strcat (score_table, "_table_"); 
			tostring(tmp, i); 
			strcat(score_table, tmp); 
			strcat(score_table, "_"); 
			tostring(tmp, j); 
			strcat(score_table, tmp); 
			strcat(score_table, ".txt"); 
			printf ("%s\n", score_table); 

			ReadTable(score_table, i, j, prob_table); 
		}
	}
	for (l = 0; l < num_data_points; l ++) { 
		tmp_pred = predictions[l]; 
		i = tmp_pred->spec_th; j = tmp_pred->freq_th; 
		if (tmp_pred->score < score_range[14]) { 
			fprintf(OUT, "%s\t%s\t%.3f\n", tmp_pred->gene, tmp_pred->go, 0.0);
			continue;
		}
		for (m = 1; m <= 14; m ++) { 
			if (tmp_pred->score >= score_range[m]) { 
				fprintf(OUT, "%s\t%s\t%.3f\n", tmp_pred->gene, tmp_pred->go, prob_table[i][j][m]);
				break; 
			}
		}
	}
	fclose(OUT);

	return 1; 
}
int cmp( const void *a ,const void *b) 
{
	return (*(PredScorePtr *)a)->score < (*(PredScorePtr *)b)->score ? 1 : -1; 
}

int reverse_cmp( const void *a ,const void *b) 
{
	return (*(PredScorePtr *)a)->score < (*(PredScorePtr *)b)->score ? -1 : 1; 
}

int wc (char any_file[]) {
	FILE *fp;
	int line_number;
	char line[400];  // length 2000 

	line_number=0;
	if ((fp = fopen(any_file, "r")) == NULL) {
		printf("\nFailure to open file %s in read mode\n", any_file);
		fflush(NULL);
		return(-1);
	} 
	while (fgets(line, sizeof(line), fp) != NULL) {
		line_number++;
	} 
	fclose(fp);
	return (line_number);
} 

int GetScores(char score_file[], PredScorePtr *predictions) {
	FILE *fp;
	int  attr_ID, num_entries, num_target_vertex_attrs;
	char line[400];  // length 2000 
	char *token, flag; 
	int i, num_data_points, all_index, n, true_index, false_index, num_true; 
	double score; 
	char mark[40], abbr[40]; 
	double prob;
	double go_freq, go_spec;
	char marker;
	int go_id;
	int gene_id;
	double evalue, max_id, avr_id; 
	PredScorePtr tmp_pred;

	if ((fp = fopen(score_file, "r")) == NULL) {
		printf("\nFailure in opening file %s in read mode!\n", score_file);
		fflush(NULL);
		return(-1);
	}
	i = 0; num_true = 0; n = 0; 
	while(fgets(line, sizeof(line), fp) != NULL) {
		n ++; 
		if ((token = strtok(line, DELIMITER)) == NULL) {
			continue;
		}
		strncpy(predictions[i]->go,token, MAX_ATTR_NAME);  
		if ((token = strtok(NULL, DELIMITER)) == NULL) {
			continue;
		}
		strncpy(predictions[i]->gene,token, MAX_ATTR_NAME);  
		if ((token = strtok(NULL, DELIMITER)) == NULL) {
			continue;
		}
		predictions[i]->score = atof(token); 
		if ((token = strtok(NULL, DELIMITER)) == NULL) {
			continue;
		}
		predictions[i]->go_freq = atof(token); 
		if ((token = strtok(NULL, DELIMITER)) == NULL) {
			continue;
		}
		predictions[i]->max_id = atof(token); 
		i ++; 
	}
	fclose(fp);

	return num_true; 
}

int ReadAttrs(int num_go_attrs, char attr_file[]) { 
	FILE *fp;
	int current, parent, parent_index, attr_ID,  max_attr_ID, num_parents;
	char line[MAX_LINE_LENGTH];  // length 2000 
	char abbr[MAX_ATTR_BRIEF]; 
	char *token, test, type;
	int i, num_types; 
	int  num_types_index; 
	int  max_attr_brief, num_entries; 

	if (! num_go_attrs) {
		printf("No go attributes to read!\n");
		fflush(NULL);
		return(0);
	}
	GO_attrs =  (VertexAttrPtr) malloc((num_go_attrs + 1)  * sizeof(VertexAttrStr)); 
	for (attr_ID = 1; attr_ID <= num_go_attrs; attr_ID ++) {
		(GO_attrs + attr_ID)->num_trues = 0; 
		(GO_attrs + attr_ID)->num_homologs = 0; 
		(GO_attrs + attr_ID)->go_spec = 0; 
		strcpy((GO_attrs+attr_ID)->name, ""); 
	}
	if ((fp = fopen(attr_file, "r")) == NULL) {
		printf("\nFailure in opening file %s in read mode!\n", attr_file);
		fflush(NULL);
		return(-1);
	}
	while(fgets(line, sizeof(line), fp) != NULL) {
		if ((token = strtok(line, DELIMITER)) != NULL) {
			attr_ID = atoi(token); 
		}
		if ((token = strtok(NULL, DELIMITER)) != NULL) {
			strncpy((GO_attrs + attr_ID)->name, token, MAX_ATTR_NAME);    // MAX_ATTR_NAME = 400 
		}
	}
	return attr_ID;

}

int get_max_ID (char any_file[]) { 
	FILE *fp;
	char line[MAX_LINE_LENGTH];
	int max_ID, ID; 
	char *token;
	if ((fp = fopen(any_file, "r")) == NULL) { 
		fprintf(stderr, "\nFailure to open file %s in read mode\n", any_file);
		fflush(NULL); 
		return(-1);
	}
	max_ID = 0; 
	while (fgets(line, sizeof(line), fp) != NULL) {
		if ((token = strtok(line, DELIMITER)) == NULL) {
			continue;
		}
		ID = atoi(token);        // index of vertex (genes)
		if (ID > max_ID) { 
			max_ID = ID; 
		}
	}
	return max_ID; 
}


int ReadGenes(int num_genes, char attr_file[]) { 
	FILE *fp;
	int current, parent, parent_index, attr_ID,  max_attr_ID, num_parents;
	char line[MAX_LINE_LENGTH];  // length 2000 
	char abbr[MAX_ATTR_BRIEF]; 
	char *token, test, type;
	int i, num_types; 
	int  num_types_index; 
	int  max_attr_brief, num_entries; 
	int gene_ID; 
	char marker; 

	if (! num_genes) {
		printf("No go attributes to read!\n");
		fflush(NULL);
		return(0);
	}
	Training_genes =  (VerticePtr) malloc((num_genes + 1)  * sizeof(VerticeStr)); 
	for (gene_ID = 1; gene_ID <= num_genes; gene_ID ++) {
		(Training_genes + gene_ID)->num_trues = 0; 
		(Training_genes + gene_ID)->num_candidate_gos = 0; 
		strcpy((Training_genes+gene_ID)->name, ""); 
	}
	if ((fp = fopen(attr_file, "r")) == NULL) {
		printf("\nFailure in opening file %s in read mode!\n", attr_file);
		fflush(NULL);
		return(-1);
	}
	while(fgets(line, sizeof(line), fp) != NULL) {
		if ((token = strtok(line, DELIMITER)) != NULL) {
			gene_ID = atoi(token); 
		}
		if ((token = strtok(NULL, DELIMITER)) != NULL) {
			marker = atoi(token); 
		}
		(Training_genes + gene_ID)->marker = marker; 
		if ((token = strtok(NULL, DELIMITER)) != NULL) {
			strncpy((Training_genes + gene_ID)->name, token, MAX_ATTR_NAME);    // MAX_ATTR_NAME = 400 
		}
	}
	return gene_ID;

}

double ReadTable (char *score_file, int avr_id_range, int freq_range, double ***prob_table) { 
	int FP, i, j, TP, AP; 
	double score, old_fp, old_tp, tmp_tp, roc_auc, roc_auc_01, fp, tp, recall, prec, area; 
	double old_prec, old_recall, fmax, th; 
	int test; 
	FILE *IN; 
	FILE *OUT; 
	PredScorePtr prediction; 
	int *true_hits;
	int *all_hits;
	double *prob, f, th1, th2; 
	char line[400], abbr[400];
	char *token; 

	printf ("id %d freq %d\n", avr_id_range, freq_range); 

	if ((IN = fopen(score_file, "r")) == NULL) {
		printf ("cannot open file to read\n"); 
		return 0;
	}
	strcpy(abbr, "raw"); 
	i = 0; 
	while(fgets(line, sizeof(line), IN) != NULL) {
		if ((token = strtok(line, DELIMITER)) == NULL) {
			continue;
		}
		if ((strncmp(token, abbr, 3)) == 0) {continue;}
		if ((token = strtok(NULL, DELIMITER)) == NULL) {
			continue;
		}
		i ++; 
		prob_table[avr_id_range][freq_range][i] = atof(token); 
	}
	fclose(IN); 

	return i; 
}


void tostring(char str[], int num)
{
	int i, rem, len = 0, n;

	n = num;
	while (n != 0)
	{
		len++;
		n /= 10;
	}
	for (i = 0; i < len; i++)
	{
		rem = num % 10;
		num = num / 10;
		str[len - (i + 1)] = rem + '0';
	}
	str[len] = '\0';
}

void ftoa(float n, char *res, int afterpoint)
{
	// Extract integer part
	int ipart = (int)n;

	// Extract floating part
	float fpart = n - (float)ipart;

	// convert integer part to string
	int i = intToStr(ipart, res, 0);

	// check for display option after point
	if (afterpoint != 0)
	{
		res[i] = '.';  // add dot

		// Get the value of fraction part upto given no.
		// of points after dot. The third parameter is needed
		// to handle cases like 233.007
		fpart = fpart * pow(10, afterpoint);

		intToStr((int)fpart, res + i + 1, afterpoint);
	}
}

int intToStr(int x, char str[], int d)
{
	int i = 0;
	while (x)
	{
		str[i++] = (x%10) + '0';
		x = x/10;
	}

	// If number of digits required is more, then
	// add 0s at the beginning
	while (i < d)
		str[i++] = '0';

	reverse(str, i);
	str[i] = '\0';
	return i;
}


void reverse(char *str, int len)
{
	int i=0, j=len-1, temp;
	while (i<j)
	{
		temp = str[i];
		str[i] = str[j];
		str[j] = temp;
		i++; j--;
	}
}
