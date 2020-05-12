import re
import sys
import vcf
import json
import subprocess
import gzip
import logging
from collections import Counter

class create_bigwig:
    def __init__(self):
        self.contig_length = {}
        pass

    def _run_cmd(self, cmd):
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

    def _convert_bedgrpah_to_bigwig(self, bedgraph_file, chr_length_file, output_bigwig_file):
        sort_cmd = "sort -k1,1 -k2,2n " + bedgraph_file + "> sorted_" + bedgraph_file
        self._run_cmd(sort_cmd)
        cmd = "./bedGraphToBigWig " + "sorted_" + bedgraph_file + " " + chr_length_file + " " + output_bigwig_file
        logging.info("Generating bigwig ..\n" +cmd + "\n")
        self._run_cmd(cmd)

    def _read_refseq_json(self, jsonfile, chrlength_file):
        with open(jsonfile) as json_file:
            data = json.load(json_file)

            try:
               with open(chrlength_file, "w") as flen:
                    for i in range(0, len(data)):
                        self.contig_length[data[i]["name"]] = data[i]["end"]
                        flen.write(data[i]["name"] + "\t" + str(data[i]["end"])+"\n")

            except IOError:
                print("Unable to write "+ chrlength_file, + " file on disk.")        

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

    def _index_gff(self, gff_file):
        zipcmd = "bgzip " + gff_file
        logging.info("zipping and indexing gff")
        self._run_cmd(zipcmd)
        indexcmd = "tabix -p gff " + gff_file + ".gz"
        "tabix -p gff <gff filename>"
        self._run_cmd(indexcmd)
        gff_info = {
            'file_path' : gff_file+ ".gz",
            'index_file_path' :  gff_file + ".gz.tbi"
        }
        return gff_info

    def _parse_vcf_data(self, params, binszie, chrlength_file):
        #vcf_filepath = self._stage_input(params)

        vcf_filepath = "athaliana.vcf.gz"        # file is validated by this point, can assume vcf_filepath is valid
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


if __name__ == "__main__":
   cw = create_bigwig()
   params= {"file":"xyz"}
   jsonfile = "refSeqs.json"
   binsize = 500000
   chrlength_file = "chr_length.tsv"
   cw._read_refseq_json(jsonfile, chrlength_file)
   vcf_info = cw._parse_vcf_data(params, binsize, chrlength_file)
   print(vcf_info)
   gff_file = "Sorted_Ptrichocarpa_444_v3.1.gene.gff3"
   cw._convert_bedgrpah_to_bigwig("updated_athaliana.vcf.bedgraph", "chr_length.tsv", "Ptrichocarpa.vcf.bw" )
   index_info = cw._index_gff(gff_file)

