import json
import sys
import os
import gzip
import subprocess
import re
import gzip
import logging
from collections import Counter

class prepare_data:

   def __init__(self):
        self.contig_length = {}
        pass

   def convert_assembly_to_refseq_json(self, assembly_json_file):
       refseq = []

       flen =  open("chrlength_file", "w")

       with open(assembly_json_file) as json_file:
            data = json.load(json_file)
            for key in data['contigs']:

                self.contig_length[data['contigs'][key]["contig_id"]] = data['contigs'][key]["length"]
                flen.write(data['contigs'][key]["contig_id"] + "\t" + str(data['contigs'][key]["length"])+"\n")

                refseq.append(
                        {"end": data['contigs'][key]["length"],
                         "length":data['contigs'][key]["length"],
                         "name":data['contigs'][key]["contig_id"],
                         "seqChunkSize":20000,
                         "start":0
                        } 
                      )

       flen.close()

       output = json.dumps(refseq)
 
       if not os.path.isdir("seq"):
          os.mkdir("seq")

       with open("seq/refSeq.json", "w") as outfile: 
            outfile.write(output) 
    
       return "refseq.json"
    
   def run_cmd(self, cmd):
       try:
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
            stdout, stderr = process.communicate()
            if stdout:
               print("ret> ", process.returncode)
               print("OK> output ", stdout)
            if stderr:
               print("ret> ", process.returncode)
               print("Error> error ", stderr.strip())

       except OSError as e:
           print("OSError > ", e.errno)
           print("OSError > ", e.strerror)
           print("OSError > ", e.filename)


   def convert_bedgrpah_to_bigwig(self, bedgraph_file, chr_length_file, output_bigwig_file):
       sort_cmd = "sort -k1,1 -k2,2n " + bedgraph_file + "> sorted_" + bedgraph_file
       self.run_cmd(sort_cmd)
       cmd = "./bedGraphToBigWig " + "sorted_" + bedgraph_file + " " + chr_length_file + " " + output_bigwig_file
       logging.info("Generating bigwig ..\n" +cmd + "\n")
       self.run_cmd(cmd)
       return output_bigwig_file


   def prepare_vcf(self, vcf_file):
    
       zip_cmd = "bgzip " + vcf_file
       self.run_cmd(zip_cmd) 

       index_cmd = "tabix -p vcf " + vcf_file + ".gz"
       self.run_cmd(index_cmd) 

       return {"vcf_file_path" : vcf_file+".gz", "index_file_path" : vcf_file + "gz.tbi"}

   def prepare_gff(self, gff_file):

       sorted_gff_cmd = "sort -k1,1 -k4,4n " + gff_file + " > " + "sorted_" + gff_file
       self.run_cmd(sorted_gff_cmd)

       zip_cmd = "bgzip " + "sorted_" + gff_file 
       self.run_cmd(zip_cmd) 

       index_gff_cmd = "tabix -p gff " + "sorted_" + gff_file + ".gz"
       self.run_cmd(index_gff_cmd)

       return {"gff_file_path" : "sorted_" + gff_file + ".gz", "index_file_path" : "sorted_" + gff_file + ".gz.tbi"}        

   def prepare_ref (self, assembly_file):
   
       index_cmd = "samtools faidx " + assembly_file 
       self.run_cmd(index_cmd)
       return {"assembly_file" : assembly_file , "assembly_index_file" : assembly_file + ".fai"}

   def update_bed_graph(self, bedgraph_file, chrlength_file):
        xdict = {}
        bedgraph_dict = {}

        try: 
           with open(bedgraph_file, "r") as f:
              for x in f:
                 x=x.rstrip()
                 data = x.split("\t")
                 xdict[data[0]] = data[2]
                 bedgraph_dict[data[0] + "-" + data[2]] = x

        except IOError:
                print("Unable to read from"+ bedgraph_file)

        ydict = {}

        try:
           with open(chrlength_file, "r") as fchr:
              for y in fchr:
                 y = y.rstrip()
                 length = y.split("\t")
                 key1 = length[0]
                 if(length[0] in xdict):
                   key2 = xdict[length[0]]
                   key = str(key1) + "-"+ str(key2) 
                   ydict[key] = length[1]
        except IOError:
                print("Unable to read from"+ chrlength_file)

         
        try: 
           with open("updated_" + bedgraph_file, "w") as fupdate:
              for key in bedgraph_dict:
                 rec = (bedgraph_dict[key]).split("\t")
                 if key in ydict:
                    fupdate.write(rec[0]+"\t"+rec[1]+"\t"+ydict[key]+"\t"+rec[3]+"\n") 
                 else:
                    fupdate.write(rec[0]+"\t"+rec[1]+"\t"+rec[2]+"\t"+rec[3]+"\n")

        except IOError:
                print("Unable to write to updated_"+ bedgraph_file, + " file on disk.") 
   

   def  parse_vcf_data(self, vcf_filepath, chrlength_file, binsize):

        reader = gzip.open(vcf_filepath, "rt")

        version = ""
        genotypes = []
        contigs = {}
        chromosomes = []
        totalvars = 0
        counts=Counter()

        counter = 0
        var_length = {}

        logging.info("Generating bedgraph file\n")
        for record in reader:
            record = record.rstrip()
            if record[0]=="#":
               if(record.startswith("##fileformat=")):
                  version=(record.split("=")[1]).replace("VCFv", "")

               if(record.startswith("#CHROM")):
                  header = record.split("\t")
                  for hd in range(9,len(header)): 
                      genotypes.append(header[hd])
               continue
           
            totalvars = totalvars +1
            counter=counter+1
            d=record.split("\t")
            bf = int(int(d[1])/binsize)
            bin_id = str(d[0]) + "\t" + str(bf)

            if d[0] not in var_length:
               var_length[d[0]] = 1 
            else:
               var_length[d[0]] = var_length[d[0]] + 1

            if d[0] not in chromosomes:
                chromosomes.append(d[0])
            counts[bin_id] +=1

        for chrm in chromosomes:
            contigs[chrm] = {
                    'contig_id': chrm,
                    'totalvariants': var_length[chrm],
                    'length': self.contig_length[chrm]
                }     
        
        #wiriting updated bedgraph to file
        bedgraph_file = vcf_filepath.replace("gz","bedgraph")

        try:
           with open(bedgraph_file, "w") as fout:
              for j,k in counts.items():

                 chromosome, bin_num = j.split("\t")
        
                 bin_start=int(bin_num)*binsize
                 bin_end = bin_start + binsize
                 fout.write (chromosome  + "\t" + str(bin_start) + "\t" + str(bin_end) + "\t" + str(k) +"\n")

        except IOError:
                print("Unable to write "+ bedgraph_file, + " file on disk.")

        self.update_bed_graph(bedgraph_file, chrlength_file) 
        


        vcf_info = {
            'version': version,
            'contigs': contigs,
            'total_variants': counter,
            'genotype_ids': genotypes,
            'chromosome_ids': chromosomes,
            'file_ref': vcf_filepath
        }

        return vcf_info     

   def prepare_tracklist_json(self, json_template, assembly_file, vcf_file, gff_file, assembly_json):
       
       refSeq_json = self.convert_assembly_to_refseq_json(assembly_json)
       assembly_info  = self.prepare_ref(assembly_file)
       vcf_info = self.prepare_vcf(vcf_file)
       gff_info = self.prepare_gff(gff_file)
       self.parse_vcf_data(vcf_info["vcf_file_path"], "chrlength_file", 500000)
       bedgraph_file = "updated_" + vcf_file + ".bedgraph"
       output_bigwig_file =  vcf_file + ".bw"
       bigwig_file = self.convert_bedgrpah_to_bigwig(bedgraph_file, "chrlength_file", output_bigwig_file)
        
       with open(json_template) as tjson:
            template_data = json.load(tjson)
            template_data["tracks"][0]["urlTemplate"] = assembly_info["assembly_file"]
            template_data["tracks"][0]["faiUrlTemplate"] = assembly_info["assembly_index_file"]
            template_data["tracks"][1]["urlTemplate"] = vcf_info["vcf_file_path"]
            template_data["tracks"][1]["faiUrlTemplate"] = vcf_info["index_file_path"]
            template_data["tracks"][2]["urlTemplate"] = gff_info["gff_file_path"]
            template_data["tracks"][2]["faiUrlTemplate"] = gff_info["index_file_path"]
            template_data["tracks"][3]["urlTemplate"] = bigwig_file

       template_output = json.dumps(template_data)

       with open("trackList.json", "w") as outfile:
            outfile.write(template_output)
    

if __name__ == "__main__":
 
  pd = prepare_data()
  pd.prepare_tracklist_json("jbrowse_template.json", "Ptrichocarpa_v3.1.assembly.fa", "data.vcf", "Ptrichocarpa_444_v3.1.gene.gff3", "poplar_assembly.json") 
    
