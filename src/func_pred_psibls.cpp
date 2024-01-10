#include <fstream>
#include <iostream>
#include <string>
#include <ctime>
#include <map>
#include <vector>
#include <set>
#include <cmath>
#include <cstdlib>
using namespace std;

struct Seq
{
	string sequence;
};

struct GO
{
	set<int> sub_type;
	vector<float> entropy_relative;
	vector<float> z_score;
	vector<int> site_list;
	float sumlog;
	float avelog;
	int n_homo;
	int n_hetero;
	int act_size;
	string gene_string;
	int model_go_i;
	float freq;
};

string aa_list="ACDEFGHILKMNPQRSTVWY-";
map<string,int> gene2index,go2index;
map<int,string> index2go;
vector<GO> all_go;
vector<string> sequences;
set<int> sequences_annotated;
string query_seq;
map<char,float> aa_freq;
set<string> answers;
map<string,float> go_specificity;
vector<map<char,int> > hetero_distribution,annotated_distribution;
string type;
float id_cutoff1,id_cutoff2;
string id_cutoff1_str,id_cutoff2_str;

float cal_id(string seq1,string seq2)
{
	int len=0,overlap=0,i;
	if(seq1.size() != seq2.size()) return 0.0;
	for(i=0;i<seq1.size();i++)
	{
		if(seq1[i]!='-' || seq2[i]!='-')
		{
			len++;
			if(seq1[i]==seq2[i]) overlap++;
		}
	}
	return float(overlap)/seq1.size();
}

void read_data(string gene_name)
{
	int i,n=0;
	float r;
	string str,title,seq,name;
	str="in_data/one_line_sequences/"+gene_name;
	ifstream ifseq(str.c_str());
	getline(ifseq,title);
	getline(ifseq,query_seq);
	ifseq.close();
	str="tmp_data/low_psibls_aln/"+gene_name;
	ifstream fin_aln(str.c_str());
	if(!fin_aln) exit(0);
	int gene_index=0;
	while(getline(fin_aln,title))
	{
		getline(fin_aln,seq);
		name=title.substr(1,6);
		if(gene2index.count(name)==0)
		{
			sequences.push_back(seq);
			gene2index[name]=gene_index;
			gene_index++;
		}
	}
	fin_aln.close();

	int go_index=0;
	string gene,go,c;
	str="tmp_data/low_assigned_GO_terms_psibls/"+gene_name;
	ifstream fin_data(str.c_str());
	while(fin_data>>gene>>go)
	{
			if(go2index.count(go)==0)
			{
				GO new_go;
				new_go.gene_string="";
				all_go.push_back(new_go);
				go2index[go]=go_index;
				index2go[go_index]=go;
				go_index++;
			}
			all_go[go2index[go]].sub_type.insert(gene2index[gene]);
			all_go[go2index[go]].gene_string+=gene;
			sequences_annotated.insert(gene2index[gene]);
	}
	fin_data.close();
}

void read_general_data()
{
	string go;
	int a,b;
	float spc;
	ifstream fin("database/annotations/GO_specificity_log10.txt");
	while(fin>>go>>a>>b>>spc)
	{
		if(go2index.count(go)>0) go_specificity[go]=spc;
	}
	fin.close();
}

void select_sequences()
{
	float r;
	int i,j;
	int n;
	int all_aa=0;
	for(vector<string>::iterator it_i=sequences.begin();it_i!=sequences.end();it_i++)
	{
		for(j=0;j<query_seq.size();j++)
		{
			char c;
			c=(*it_i)[j];
			if(c>='A' && c<='Z')
			{
				aa_freq[c]++;
				all_aa++;
			}
		}
	}
	for(i=0;i<20;i++) aa_freq[aa_list[i]]=(aa_freq[aa_list[i]]+0.05)/(all_aa+1);
	
	vector<int> v_unique;
	for(i=0;i<all_go.size();i++)
	{
		bool same=false;
		for(vector<int>::iterator it_k=v_unique.begin();it_k!=v_unique.end();it_k++)
		{
			if(all_go[i].gene_string==all_go[*it_k].gene_string)
			{
				all_go[i].model_go_i=(*it_k);
				same=true;
				break;
			}
		}
		if(!same)
		{
			v_unique.push_back(i);
			all_go[i].model_go_i=-1;
		}
	}
}

void get_hetero_distribution()
{
	for(int i=0;i<query_seq.size();i++)
	{
		map<char,int> aa_n_hetero,aa_n_annotated;
		for(set<int>::iterator it_j=sequences_annotated.begin();it_j!=sequences_annotated.end();it_j++) aa_n_annotated[sequences[*it_j][i]]+=1;
		annotated_distribution.push_back(aa_n_annotated);
	}
}

void cal_entropy_relative()
{
	if(query_seq.size()<=10)
	{
		for(int igo=0;igo<all_go.size();igo++)
		{
			if(all_go[igo].sub_type.size()==0 || all_go[igo].model_go_i>=0) continue;
			for(int i=0;i<query_seq.size();i++) all_go[igo].site_list.push_back(i);
		}
		return;
	}
	for(int igo=0;igo<all_go.size();igo++)
	{
		if(all_go[igo].sub_type.size()==0 || all_go[igo].model_go_i>=0) continue;
		float s=0,s2=0,sd,z;
		vector<float> single_position;
		for(int i=0;i<query_seq.size();i++)
		{
			map<char,int> aa_n;
			float h=0,h1=0,h2=0,v;
			int n=0;
			for(set<int>::iterator it_j=all_go[igo].sub_type.begin();it_j!=all_go[igo].sub_type.end();it_j++) aa_n[sequences[*it_j][i]]+=1;
			for(int j=0;j<20;j++)
			{
				float p,q,r,s;
				p=(float(aa_n[aa_list[j]])+float(aa_n['-'])/20+aa_freq[aa_list[j]])/(all_go[igo].sub_type.size()+1);
				q=(float(annotated_distribution[i][aa_list[j]])+float(annotated_distribution[i]['-'])/20-float(aa_n[aa_list[j]])-float(aa_n['-'])/20+aa_freq[aa_list[j]])/(sequences_annotated.size()-all_go[igo].sub_type.size()+1);
				if(p>0)
				{
					if(q>0) h-=p*log(p/q);
					else h-=10;
				}
			}
			v=h;
			all_go[igo].entropy_relative.push_back(v);
			s+=v;
			s2+=v*v;
		}
		s/=query_seq.size();
		s2/=query_seq.size();
		sd=sqrt(s2-s*s);
		for(int i=0;i<query_seq.size();i++)
		{
			z=(s-all_go[igo].entropy_relative[i])/sd;
			all_go[igo].z_score.push_back(z);
			if(z>1) all_go[igo].site_list.push_back(i);
		}
		if(all_go[igo].site_list.size()<10)
		{
			vector<float> tmp_v;
			vector<int> tmp_v_i;
			for(int i=0;i<query_seq.size();i++)
			{
				int start;
				float z=all_go[igo].z_score[i];
				if(tmp_v.size()>10)
				{
					tmp_v[10]=z;
					tmp_v_i[10]=i;
					start=10;
				}
				else
				{
					tmp_v.push_back(z);
					tmp_v_i.push_back(i);
					start=tmp_v.size()-1;
				}
				for(int j=start;j>0;j--)
				{
					if(tmp_v[j]>tmp_v[j-1])
					{
						float tmpf;
						int tmpi;
						tmpf=tmp_v[j];
						tmp_v[j]=tmp_v[j-1];
						tmp_v[j-1]=tmpf;
						tmpi=tmp_v_i[j];
						tmp_v_i[j]=tmp_v_i[j-1];
						tmp_v_i[j-1]=tmpi;
					}
					else break;
				}
			}
			all_go[igo].site_list.clear();
			for(int i=0;i<10;i++) all_go[igo].site_list.push_back(tmp_v_i[i]);
		}
	}
}

void cal_freq(string filename)
{
	string str;
	str="tmp_data/results_psibls/"+filename;
	ofstream fout(str.c_str());
	for(int igo=0;igo<all_go.size();igo++)
	{
		if(all_go[igo].model_go_i>=0)
		{
			fout<<index2go[igo]<<'\t'<<all_go[all_go[igo].model_go_i].sumlog<<'\t'<<all_go[all_go[igo].model_go_i].avelog<<'\t'<<all_go[all_go[igo].model_go_i].freq<<'\t'<<all_go[all_go[igo].model_go_i].n_homo<<'\t'<<all_go[all_go[igo].model_go_i].n_hetero<<'\t'<<query_seq.size()<<'\t'<<all_go[all_go[igo].model_go_i].act_size<<'\t'<<1<<'\t'<<go_specificity[index2go[igo]]<<endl;
			continue;
		}
		int n_homo=0,n_hetero=0;
		float final_score=0;
		float si,max_si=0,ave_si=0;
		string max_gene="";
		int n_si=0;
		int i;
		map<int,map<char,float> > pos_score;
		for(vector<int>::iterator it_k=all_go[igo].site_list.begin();it_k!=all_go[igo].site_list.end();it_k++)
		{
			n_homo=0,n_hetero=0;
			float freq_homo=0,freq_hetero=0;
			map<char,int> count_homo,count_hetero;
			for(set<int>::iterator it_j=all_go[igo].sub_type.begin();it_j!=all_go[igo].sub_type.end();it_j++)
			{
					n_homo++;
					count_homo[sequences[*it_j][*it_k]]++;
			}
			freq_homo=(float(count_homo[query_seq[*it_k]])+float(count_homo['-'])/20+aa_freq[query_seq[*it_k]])/(n_homo+1);
			freq_hetero=(float(annotated_distribution[*it_k][query_seq[*it_k]])+float(annotated_distribution[*it_k]['-'])/20-count_homo[query_seq[*it_k]]-float(count_homo['-'])/20+aa_freq[query_seq[*it_k]])/(sequences_annotated.size()-n_homo+1);
			final_score+=log(freq_homo/freq_hetero);
		}
		n_hetero=sequences.size();
		float ave_score,e_value;
		ave_score=final_score/all_go[igo].site_list.size();
		all_go[igo].sumlog=final_score;
		all_go[igo].avelog=ave_score;
		all_go[igo].n_homo=n_homo;
		all_go[igo].n_hetero=n_hetero;
		all_go[igo].act_size=all_go[igo].site_list.size();
		float freq;
		freq=float(all_go[igo].sub_type.size())/sequences_annotated.size();
		all_go[igo].freq=freq;
		fout<<index2go[igo]<<'\t'<<final_score<<'\t'<<ave_score<<'\t'<<freq<<'\t'<<n_homo<<'\t'<<n_hetero<<'\t'<<query_seq.size()<<'\t'<<all_go[igo].site_list.size()<<'\t'<<1<<'\t'<<go_specificity[index2go[igo]]<<endl;
	}
	fout.close();
}

int main(int argc,char* argv[])
{
	read_data(argv[1]);
	read_general_data();
	select_sequences();
	get_hetero_distribution();
	cal_entropy_relative();
	cal_freq(argv[1]);
	return 0;
}
