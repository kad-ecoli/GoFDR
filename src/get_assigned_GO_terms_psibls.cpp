#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <stdlib.h>
using namespace std;

struct Seq
{
	vector<string> all_pro_order;
	set<string> all_pro;
	map<string,set<string> > go_terms;
};

int main(int argc,char *argv[])
{
	string f,line,name,s;
	int i;
	map<string,Seq> list;
	set<string> all_seq;

	ifstream ifls(argv[1]);
	while(ifls>>s)
	{
		Seq tmp_seq;
		f="tmp_data/psibls_aln/"+s;
		ifstream ifaln(f.c_str());
		if(!ifaln) continue;
		while(getline(ifaln,line))
		{
			if(line[0]=='>')
			{
				name=line.substr(1,6);
				tmp_seq.all_pro_order.push_back(name);
				tmp_seq.all_pro.insert(name);
				all_seq.insert(name);
			}
		}
		ifaln.close();
		list[s]=tmp_seq;
	}
	ifls.close();


		string pro,go;
		f="database/annotations/noiea_confirmed_iea_annotations.txt";
		ifstream ifanno(f.c_str());
		while(ifanno>>pro>>go)
		{
			if(all_seq.count(pro)==0) continue;
			for(map<string,Seq>::iterator it_i=list.begin();it_i!=list.end();it_i++)
			{
				if(it_i->second.all_pro.count(pro)>0) list[it_i->first].go_terms[pro].insert(go);
			}
		}
		ifanno.close();

	for(map<string,Seq>::iterator it_k=list.begin();it_k!=list.end();it_k++)
	{
		f="tmp_data/assigned_GO_terms_psibls/";
		f+=it_k->first;
		ofstream outfile(f.c_str());
		for(vector<string>::iterator it_i=it_k->second.all_pro_order.begin();it_i!=it_k->second.all_pro_order.end();it_i++)
		{
			for(set<string>::iterator it_j=it_k->second.go_terms[*it_i].begin();it_j!=it_k->second.go_terms[*it_i].end();it_j++) outfile<<*it_i<<'\t'<<*it_j<<endl;
		}
		outfile.close();
	}

	return 0;
}

