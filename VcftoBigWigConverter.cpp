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
  std::cout << "VcftoBigWigConverter <chrlen file> <input file> <file type (gff/vcf)> <outbut file>\n";
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

map<string, string > lenmap;

void read_chrlenggth(const char* chrlen_file)
{
    ifstream lenfile(chrlen_file);
    if(!lenfile.is_open())
    {
       cout<<"could not open the chrlen file\n";
       exit(1);
    }
    char chrlenline[500];
    while(lenfile)
    {  
        lenfile.getline(chrlenline,500);
        vector<const char* > vec;
        char *pch;
        pch = strtok(chrlenline,"\t ");
        if(lenfile)
        {
           while(pch != NULL )
           {
              vec.push_back(pch);
              pch = strtok(NULL,"\t ");
           }
           lenmap[vec[0]] = vec[1];     
        }
    } 
}


void convert_bedgraph_to_bigwig(string bedgraph_file, string chrlength_file, string output_file)
{
    string command = "./bedGraphToBigWig " + bedgraph_file  + " " + chrlength_file + " " +  output_file;
    system(command.c_str()); 
}

map <string,int,ltstr> IndexbinMap;

int main(int argc, char* argv[])
{
    if(argc != 5)
       manual();
   
    read_chrlenggth(argv[1]);
    
    int start;
    const char* ftype;
   
    if(strcmp(argv[3],"gff") == 0)
    {
       start = 3;
       ftype = "gene";
    }
    else if(strcmp(argv[3],"vcf") == 0)
    {
       start = 1;
       ftype = "snp";
    }
   
    std::multimap <string, string> chrmap;
    std::set<string> chrset;
    ifstream infile(argv[2]);
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
 
          map<string, int >::iterator end_iter = IndexbinMap.end();
          
          string rbegin = IndexbinMap.rbegin()->first;
          long int value = IndexbinMap.rbegin()->second;
          char last_element[500];
          strcpy(last_element, rbegin.c_str());
          
          IndexbinMap.erase (last_element);

          vector<const char*> tokenmap;
          char *token;
          
          token = strtok(last_element,"\t ");
       
          while(token != NULL )
          {
              tokenmap.push_back(token);
              token = strtok(NULL,"\t ");
          }
                  
          string min_length = tokenmap[0];
          string contig_length = tokenmap[1];

          string key;
          key += min_length;
          key += "\t"; 
          key += lenmap[chrnum];
          
          IndexbinMap.insert(std::pair<string, int>(key, value));
   
          ofstream outfile(argv[4]);
  
          for(map<string, int >::iterator ii= IndexbinMap.begin(); ii!= IndexbinMap.end(); ++ii)
              outfile <<chrnum<<"\t"<< (*ii).first << "\t" << IndexbinMap[(*ii).first]<<endl;
          IndexbinMap.clear(); 

          outfile.close();
     } 
   
     convert_bedgraph_to_bigwig(argv[4] , argv[1], string(argv[2]) + ".bw");
     return 0;
}
