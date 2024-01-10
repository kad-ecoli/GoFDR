#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <cstdlib>
using namespace std;

map<string,string> categories;

void read_GO_categories()
{
	string go,c;
	ifstream fin("database/annotations/categories");
	while(fin>>go>>c) categories[go]=c;
	fin.close();
}

int main(int argc,char *argv[])
{
	float s,b,c;
	string gene,str,type;
	string a,d,e,f,g,h,i;
	int j;
	string categories_list[3]={"BP","MF","CC"};
	read_GO_categories();
	ifstream ifls("in_data/query_list");
	while(ifls>>gene)
	{
		ofstream of1[3],of0[3];
		for(j=0;j<3;j++)
		{
			str="tmp_data/results_psibls_"+categories_list[j]+"/"+gene;
			of1[j].open(str.c_str());
		}
		str="tmp_data/results_psibls/"+gene;
		ifstream if1(str.c_str());
		while(if1>>a>>b>>s>>c>>d>>e>>f>>g>>h>>i)
		{
			for(j=0;j<3;j++)
			{
				if(categories[a]==categories_list[j]) of1[j]<<a<<'\t'<<b<<'\t'<<s<<'\t'<<c<<'\t'<<d<<'\t'<<e<<'\t'<<f<<'\t'<<g<<'\t'<<h<<'\t'<<i<<endl;
			}
		}
		if1.close();
	}
	ifls.close();
	return 0;
}

