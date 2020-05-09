#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cstring>
#include <sstream>
#include <set>

void manual()
{
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~Usage:~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
  std::cout << "createbin <input file> <file type (gff/vcf)> > <outbut file>\n";
  exit(1);
}

using namespace std;
struct ltstr
{
  bool operator()(string s1, string s2) const
  {
    int mark1 = s1.find("-");
    int mark2 = s2.find("-");
    return atoi(s1.substr(0,mark1-1).c_str()) < atoi(s2.substr(0,mark2-1).c_str());
  }
};

map <string,int,ltstr> IndexbinMap;

int main(int argc, char* argv[])
{
    if(argc != 3)
       manual();
   
    int start;
    const char* ftype;
   
    if(strcmp(argv[2],"gff") == 0)
    {
       start = 3;
       ftype = "gene";
    }
    else if(strcmp(argv[2],"vcf") == 0)
    {
       start = 1;
       ftype = "snp";
    }
   
    std::multimap <string, string> chrmap;
    std::set<string> chrset;
    ifstream infile(argv[1]);
    if(!infile.is_open())
    {
       cout<<"could not open the feature file\n";
       exit(1);
    }
    char linee[45000],rec[45000];
    while(infile)
    {
        infile.getline(linee,45000); 
        strcpy(rec,linee);
        if(infile)
        {
           if(rec[0] != '#')
           {
               char *pch;
               pch = strtok(linee,"\t ");
               vector<const char* > vecc;
               while(pch != NULL )
               {
                   vecc.push_back(pch);
                   pch = strtok(NULL,"\t ");
               }

               chrset.insert(vecc[0]);
              
               if(strcmp(ftype,"gene") == 0)
               { 
                  if((strcmp(vecc[2], "gene") == 0))
	          { 
                     chrmap.insert(std::pair<string,string>(vecc[0], vecc[start]));
                  } 
               } 
               else 
               { 
                   chrmap.insert(std::pair<string,string>(vecc[0], vecc[start]));
               } 
           }
	}
    }

    std::multimap<string, string>::iterator itmap;

    for (std::set<string>::iterator it= chrset.begin(); it!= chrset.end(); ++it)
    {
         string chrnum =  *it;
         std::pair <std::multimap<string, string>::iterator, std::multimap<string, string>::iterator> ret;
         ret = chrmap.equal_range(chrnum);
         for (std::multimap<string,string>::iterator ite=ret.first; ite!=ret.second; ++ite)
         { 
               const char* lline = (ite->second).c_str();

               size_t offset = (size_t)atol((ite->second).c_str());
               
               int i = offset/500000;
               size_t min = i*500000;
               size_t max = i*500000 + 500000;
               
               std::ostringstream minval, maxval;
               minval << min;
               maxval << max;

               string bin;
               bin += minval.str();
               bin += "\t";
               bin += maxval.str();

               IndexbinMap[bin]++;
          }
 
          for(map<string, int >::iterator ii= IndexbinMap.begin(); ii!= IndexbinMap.end(); ++ii)
              cout <<chrnum<<"\t"<< (*ii).first << "\t" << IndexbinMap[(*ii).first]<<endl;
 
          IndexbinMap.clear(); 
     } 
   
     return 0;
}
