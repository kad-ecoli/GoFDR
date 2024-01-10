#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <cstdlib>
#include <cmath>
#include <stdlib.h>
using namespace std;

int main(int argc,char *argv[])
{
	map<string,set<string> > children;
	string str,gene,go,child,type,query;
	float p;
	type=argv[1];
	query=argv[2];
	ifstream ifls("database/annotations/list_child.txt");
	while(ifls>>go)
	{
		str="database/annotations/child_node/"+go;
		ifstream ifchi(str.c_str());
		while(ifchi>>child) children[go].insert(child);
		ifchi.close();
	}
	ifls.close();
		set<string> all_gene;
		map<string,map<string,float> > result;
		str="out_data/"+type+"/"+query+"_adjusted.txt";
		ifstream if1(str.c_str());
		while(if1>>gene>>go>>p)
		{
			result[gene][go]=p;
			all_gene.insert(gene);
		}
		if1.close();
		str="out_data/"+type+"/"+query+"_updated.txt";
		ofstream of1(str.c_str());
		for(set<string>::iterator it_i=all_gene.begin();it_i!=all_gene.end();it_i++)
		{
			for(map<string,float>::iterator it_j=result[*it_i].begin();it_j!=result[*it_i].end();it_j++)
			{
				if(children.count(it_j->first)==0)
				{
					of1<<*it_i<<'\t'<<it_j->first<<'\t'<<it_j->second<<endl;
					continue;
				}
				float max=it_j->second;
				for(map<string,float>::iterator it_k=result[*it_i].begin();it_k!=result[*it_i].end();it_k++)
				{
					if(children[it_j->first].count(it_k->first)>0 && it_k->second>max)
					{
						max=it_k->second;
					}
				}
				of1<<*it_i<<'\t'<<it_j->first<<'\t'<<max<<endl;
			}
		}
		of1.close();
	return 0;
}
