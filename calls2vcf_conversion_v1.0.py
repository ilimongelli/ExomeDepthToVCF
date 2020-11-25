# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 13:00:47 2020

@author: Federica De Paoli (federica.depaoli02@universitadipavia.it)
"""

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 12:01:22 2020
Modified on Tue Apr 28 09:00:00 2020
@author: Federica De Paoli (federica.depaoli02@universitadipavia.it)
"""
# Conversion of exomedepth and convading calls to vcf file

import argparse
import subprocess
import pandas as pd
import datetime
from argparse import RawTextHelpFormatter


def manage_bed_calls(bed,exome_depth):
    dict_bed={}
    with open(bed, 'r') as bedFile:
        for line in bedFile:
            row=line.rstrip().split("\t")
            if not row[0] in dict_bed.keys():
                dict_bed[row[0]]={}
                prev_stop=row[2]
            else:
                dict_bed[row[0]][prev_stop]=row[1]
                prev_stop=row[2]
    bedFile.close()
    with open(exome_depth,'r') as infile:
        count=0
        prev_start=-1
        prev_target_stop=-1
        prev_stop=-1
        prev_chr=""
        prev_normal_RC=0
        prev_cnv_RC=0
        chromosomes=[]
        starts=[]
        ends=[]
        types=[]
        quals=[]
        readratios=[]
        dp=[]
        for line in infile:
            if count==0:
                if(len(line.rstrip().split(",")))>1:
                    sep=","
                elif(len(line.rstrip().split("\t")))>1:
                    sep="\t"
                else:
                    sep=","
                    
                count+=1
                pass
            else:
                row=line.rstrip().split(sep)
                row[4]=int(row[4])-1
                if row[4]==prev_start and row[6]==prev_chr:
                    row[4]=dict_bed[row[6]][prev_stop]
                    prev_stop=row[5]
                    row[0]=prev_target_stop+1
                    normal_RC=int(row[9])-prev_normal_RC
                    cnv_RC=int(row[10])-prev_cnv_RC
                    prev_normal_RC=int(row[9])
                    prev_cnv_RC=int(row[10])
                    prev_target_stop=int(row[1])
                    row[10]=cnv_RC
                    row[11]=float(cnv_RC)/normal_RC
                else:
                    prev_start=int(row[4])
                    prev_stop=row[5]
                    prev_chr=row[6]
                    prev_normal_RC=int(row[9])
                    prev_cnv_RC=int(row[10])
                    prev_target_stop=int(row[1])
                chromosomes.append(str(row[6]))
                starts.append(int(row[4]))
                ends.append(int(row[5]))
                if(row[2]=="deletion"):
                    types.append("DEL")
                else:
                    types.append("DUP")
                quals.append(float(row[8]))
                readratios.append(float(row[11]))
                dp.append(int(row[10]))
            all_calls=pd.DataFrame(data={
                    'chromosome':chromosomes, 'start':starts, 'end': ends, 'type':types,  'qual':quals, 'readratio':readratios, 'dp':dp
                    })
    return all_calls
                    
def bed2vcf(all_calls,fasta,output,sample,upper, lower,qualthreshold, bedtools):
    if all_calls.empty:
         with open(output,'w') as output_file:
            output_file.write(header(upper,lower,qualthreshold,sample,fasta))
    else:
        vcf_fields=pd.DataFrame(columns=["CHR","POS","ID","REF","ALT", "QUAL","FILTER","INFO","FORMAT"])    
        vcf_fields["CHR"]=all_calls["chromosome"]   
        vcf_fields["POS"]=all_calls["start"]+1
        vcf_fields["ID"]="."
        vcf_fields["ALT"]="<"+all_calls["type"]+">"
        vcf_fields["QUAL"]="."
        vcf_fields["FILTER"]=""
        
        for index, row in all_calls.iterrows():
            vcf_fields.loc[index,'INFO']=";".join(["IMPRECISE",
                      "=".join(["SVTYPE",row["type"]]),
                      "=".join(["END",str(row["end"])]), 
                      "=".join(["SVLEN",str(row["end"]-row["start"])]),
                      "=".join(["READRATIO",str(row["readratio"])]),
                      "=".join(["BF",str(row["qual"])]),
                      "=".join(["DP",str(row["dp"])])
                      ])
            if float(row['qual'])<qualthreshold:
                vcf_fields.loc[index,"FILTER"]="low_bayes_factor"
            if row['readratio']>lower and row['readratio']<upper:
                if vcf_fields.loc[index,"FILTER"]!="":
                    vcf_fields.loc[index,"FILTER"]=";".join([vcf_fields.loc[index,"FILTER"],"reads_ratio_close_to_1"])
                else:
                    vcf_fields.loc[index,"FILTER"]="reads_ratio_close_to_1"
            if not vcf_fields.loc[index,"FILTER"]:
                vcf_fields.loc[index,"FILTER"]="PASS"
                
        bed_calls=vcf_fields.copy(deep=True)[["CHR","POS"]]
        bed_calls.loc[:,"END"]=bed_calls.loc[:,"POS"]+1
        bed_calls.loc[:,"END"] = pd.to_numeric(bed_calls.loc[:,"END"],downcast='integer')
        bed_calls.loc[:,"POS"] = pd.to_numeric(bed_calls.loc[:,"POS"],downcast='integer')
        #bed_calls["CHR"] = pd.to_numeric(bed_calls["CHR"],downcast='integer')
        bed_calls.to_csv ('tmp.bed', index = None, header=False,sep="\t")
        #infer reference sequence for REF field 
        bedtools_command=bedtools + " getfasta -fi " + fasta + " -bed tmp.bed -fo tmp.seq.txt "
        subprocess.call(bedtools_command, shell=True)
        with open('tmp.seq.txt','r') as seq_file:
            count=0
            for line in seq_file:
                if line.startswith(">"):
                    pass
                else:
                    vcf_fields.loc[count,"REF"]=line.rstrip()
                    count+=1
    
        seq_file.close()    
        with open(output,'w') as output_file:
            output_file.write(header(upper,lower,qualthreshold,sample,fasta))
            for index,row in vcf_fields.iterrows():
                output_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (row["CHR"],row["POS"],row["ID"],row["REF"],row["ALT"],row["QUAL"],row["FILTER"],row["INFO"],".","."))
        rm_command1="rm tmp.bed"
        rm_command2="rm tmp.seq.txt"
        subprocess.call(rm_command1,shell=True)
        subprocess.call(rm_command2,shell=True)
    
def header(upper,lower,qualthreshold,sample,fasta):
    TS_NOW = datetime.datetime.now()
    
    VCF_HEADER = """##fileformat=VCFv4.2
##fileDate={filedate}
##reference={reference_fasta}
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##FILTER=<ID=reads_ratio_close_to_1,Description="Filtered due to read ratio between {lower} and {upper}">
##FILTER=<ID=low_bayes_factor,Description="Filtered due to Bayes Factor < {qualthreshold}">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=READRATIO,Number=.,Type=Float,Description="Reads ratio between individual sample and baseline">
##INFO=<ID=BF,Number=.,Type=Float,Description="Bayes Factor for CNV call">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##FORMAT=<ID=.,Number=1,Type=String,Description="No format available">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n""".format( filedate=TS_NOW.strftime( "%Y%m%d" ), lower=lower, upper=upper, qualthreshold=qualthreshold, sample=sample, reference_fasta=fasta )
    
    return VCF_HEADER
#-------------------------------------------------------------------------------------------
#     MAIN
#-------------------------------------------------------------------------------------------
if __name__ == '__main__':

    p = argparse.ArgumentParser(description=__doc__,formatter_class=RawTextHelpFormatter)
    
    p.add_argument('-e', '--exomedepth', type=str, required=True, help='Path to ExomeDepth input file')
    p.add_argument('-b', '--bed', type=str, required=True, help='Path to Bed CNV file')
    p.add_argument('-o', '--output', type=str, required=True, help='Path to vcf output file')
    p.add_argument('-f', '--fasta', type=str, required=True, help='Path to reference genome fasta file')
    p.add_argument('-ut','--upper_threshold', type=float, required=False, help='Read ratio threshold for duplications', default=1.2)
    p.add_argument('-lt','--lower_threshold', type=float, required=False, help='Read ratio threshold for deletions', default=0.8)
    p.add_argument('-qt','--qual_threshold', type=float, required=False, help='Quality threshold for Bayes Factor', default=6)
    p.add_argument('-s','--sample', type=str, required=True, help='Sample name')
    p.add_argument('-t','--bedtools', type=str, required=True, help='Path to bedtools binary')
    global args
    args = p.parse_args()
    
    all_calls=manage_bed_calls(args.bed,args.exomedepth)
    
    
    # Create VCF file from BED
    bed2vcf(all_calls,args.fasta,args.output, args.sample, args.upper_threshold, args.lower_threshold, args.qual_threshold, args.bedtools)
    