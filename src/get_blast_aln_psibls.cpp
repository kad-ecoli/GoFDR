#include <fstream>
#include <iostream>
#include <string>
#include <ctime>
#include <map>
#include <vector>
#include <set>
#include <cstdlib>
using namespace std;

string read_value(string line,string start,string end)
{
	int i,j;
	string s="";
	i=line.find(start);
	j=line.find(end);
	s=line.substr(i+start.size(),j-i-start.size());
	return s;
}

string match_seq(string qseq,string hseq,int qstart,int qend,int hstart,int hend,int length)
{
	string s="";
	int i,hit_index=0;
	s.append(qstart-1,'-');
	for(i=0;i<qseq.size();i++) if(qseq[i]!='-') s+=hseq[i];
	s.append(length-qend,'-');
	return s;
}

void get_aln(string gene)
{
	string str,line;
	str="tmp_data/alignments_psibls/"+gene;
	ifstream fin(str.c_str());
	if(! fin)
	{
		cerr<<"Input Error"<<endl;
		return;
	}
	string round,ter_round;
	string qlen,name,e_value,qstart,qend,hstart,hend,id,seq,tmp_query_seq;
	fin.clear();
	fin.seekg(0);
	str="tmp_data/psibls_aln/"+gene;
	ofstream ofblast(str.c_str());
	while(getline(fin,line))
	{
		if(line.find("<BlastOutput_query-len>")==2) qlen=read_value(line,"<BlastOutput_query-len>","</BlastOutput_query-len>");
		if(line.find("<Hit_def>")==10) name=line.substr(28,6);
		if(line.find("<Hsp_evalue>")==14 && name!="") e_value=read_value(line,"<Hsp_evalue>","</Hsp_evalue>");
		if(line.find("<Hsp_query-from>")==14 && name!="") qstart=read_value(line,"<Hsp_query-from>","</Hsp_query-from>");
		if(line.find("<Hsp_query-to>")==14 && name!="") qend=read_value(line,"<Hsp_query-to>","</Hsp_query-to>");
		if(line.find("<Hsp_hit-from>")==14 && name!="") hstart=read_value(line,"<Hsp_hit-from>","</Hsp_hit-from>");
		if(line.find("<Hsp_hit-to>")==14 && name!="") hend=read_value(line,"<Hsp_hit-to>","</Hsp_hit-to>");
		if(line.find("<Hsp_identity>")==14 && name!="") id=read_value(line,"<Hsp_identity>","</Hsp_identity>");
		if(line.find("<Hsp_qseq>")==14 && name!="") tmp_query_seq=read_value(line,"<Hsp_qseq>","</Hsp_qseq>");
		if(line.find("<Hsp_hseq>")==14 && name!="")
		{
			seq=read_value(line,"<Hsp_hseq>","</Hsp_hseq>");
			int n_real=0,n_gap=0;
			string matched_seq,title;
			matched_seq=match_seq(tmp_query_seq,seq,atoi(qstart.c_str()),atoi(qend.c_str()),atoi(hstart.c_str()),atoi(hend.c_str()),atoi(qlen.c_str()));
			title=">"+name+"\t"+e_value+"\t"+id+"\t"+qlen;
			if(name.substr(0,3)!="UPI") ofblast<<title<<endl<<matched_seq<<endl;
			name="";
		}
	}
	ofblast.close();
	fin.close();
}

int main(int argc,char *argv[])
{
	get_aln(argv[1]);
	return 0;
}

