#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  8 13:38:52 2017

@author: yuhan
"""



import sys,  getopt

def parameter_used():
    print("detail information is in the README.txt.")
    print("parameter_used: python %s -option <argument>" %sys.argv[0])
    print("   -U     <STRING>    unpaired-fastq")
    print("   -1     <STRING>    paired-fastq-R1")
    print("   -2     <STRING>    paired-fastq-R2")
    print("   -m     <STRING>    database mask type: hybrid-mask (default), repeat-mask, homology-mask")
    print("   -o     <STRING>    prefix of output file")
    print("   -p     <INT>       number of threaded, default is 1")
    print("   -cov    <INT>    minimal coverage threshold to determine HPV present, default is 150 bp.")

    print("   -h/--help                     ")

try:
    opts, args = getopt.getopt( sys.argv[1:], "U:1:2:m:o:p:cov:h", ["help" ] )
except getopt.GetoptError:
    print("miss parameters, option error!")
    parameter_used()
    sys.exit(2)


#default settings
nthreades="1"
database_type="hybrid-mask"
min_cov=str(150)

## retrieve settings
for opt, val in opts:
    if opt in ( "-h", "--help" ):
        parameter_used()
        sys.exit(1)
    else:
        if opt in ( "-U", ):
            unpair = val
        if opt in ( "-1", ):
            R1_pair = val
        if opt in ( "-2", ):
            R2_pair = val
        if opt in ( "-m", ):
            database_type = val
        if opt in ( "-o", ):
            outprefix = val
        if opt in ( "-p", ):
            nthreades = str(val)
        if opt in ( "-cov", ):
            min_cov = str(val)


try: 
    R1_pair, R2_pair, outprefix, database_type
    sequence_type='paired alignments'
except:
    try: 
        unpair,outprefix, database_type
        sequence_type='unpaired alignments'
    except:
      print("\ninput fastq files error!")
      parameter_used()
      
      
      
print (sequence_type)
import os

HPViewer_path = os.path.dirname(os.path.realpath(__file__))

HPViewer_path=HPViewer_path+'/'
#this path should be ended with '/'

#chosee database according to the selected mode
if database_type=='repeat-mask' or database_type=='hybrid-mask' :
  database_1=HPViewer_path+'database/repeat-mask/HPV_mask'
  HPV_length=HPViewer_path+'database/repeat-mask/HPV_mask_length.txt'
elif database_type=='homology-mask':
  database_1=HPViewer_path+'database/homology-mask/HPV_homology_mask'
  HPV_length=HPViewer_path+'database/homology-mask/HPV_homology_mask_length.txt'
else:
  print ("-m parameter has problems. -m repeat-mask or -m homology-mask")

print ("HPView mode: "+database_type)

database_2=HPViewer_path+'database/homology-mask/HPV_homology_mask'
import subprocess
from subprocess import call

call("mkdir "+outprefix,shell=True)
call("mkdir "+outprefix+'/temp',shell=True)


output_sample=outprefix.split('/')[-1]
output_location=outprefix
outprefix=outprefix+'/temp/'+output_sample


#alignment unpair and paried reads
def align_unpaired( database, unpair, outprefix, nthreades):

    try:
        subprocess.check_output(['bowtie2', '-h'])
    except OSError:
        raise RuntimeError('bowtie2 not found; check if it is installed and in $PATH\n')

    # stream output from bowtie2
    bowtie_args = ['bowtie2 -x '+ database+ ' -U '+ unpair +' -p '+nthreades+" --quiet --no-unal -S " +outprefix+".sam" ]
    for command in bowtie_args:
        call(command,shell=True)

        
def align_paired( database, R1_pair,R2_pair, outprefix, nthreades, flags=("--quiet","--no-unal")):


    # check that we have access to bowtie2
    try:
        subprocess.check_output(['bowtie2', '-h'])
    except OSError:
        raise RuntimeError('bowtie2 not found; check if it is installed and in $PATH\n')

    # stream output from bowtie2
    bowtie_args = ['bowtie2 -x '+ database+ ' -1 '+ R1_pair+' -2 '+R2_pair+' -p '+nthreades+" --quiet --no-unal -S " +outprefix+".sam" ]
    for command in bowtie_args:
        call(command,shell=True)

if sequence_type=='paired':
    align_paired(database_1, R1_pair,R2_pair,outprefix,nthreades)
else:
    align_unpaired(database_1, unpair, outprefix,nthreades)



#samtools to select potential present HPV which at least is covered for 150 bp (min_cov)
def run_samtools (outprefix,min_cov):
    
    try:
        subprocess.check_output("samtools view -? ",shell=True)
    except OSError:
        raise RuntimeError('samtools not found')

    samtools_arg=['samtools view -bS '+outprefix+'.sam'+' | samtools sort  -o '+outprefix+'.bam',\
    'samtools depth '+outprefix+'.bam'+ " | cut -f1 | uniq -c | rev| cut -d ' ' -f 1,2  | rev | awk '{if ($1 > "+min_cov+") {print $2}}'  > " + outprefix+'_L1_id.txt', \
    'samtools view '+outprefix+'.bam'+" | cut -f3 |uniq -c | rev | cut -d ' ' -f 1,2   |  rev | grep   -Ff "+ outprefix+"_L1_id.txt  >" + outprefix+'_summary_L1.txt']
    for command in samtools_arg:
        call(command,shell=True)

run_samtools(outprefix,min_cov)



def summary_1(FILE):
  HPV_homo=dict()
  HPV_homo_id=list()
  homo_matirx=HPViewer_path+'database/homology_distance_matrix.csv'
  with open (homo_matirx) as f:
        HPV_homo_id = f.readline().rstrip().split(',')
        del HPV_homo_id[0]
        for line in f:
            line=line.rstrip().split(',')
            HPV_homo[line[0]]=list()
            for l in line[1:len(line)]:
              HPV_homo[line[0]].append(float(l))
      

  def hpv_neighbor(hpv):
    min_HPV_list=list()
    if  hpv not in HPV_homo:
       return min_HPV_list 
    else:
      for i in range(len(HPV_homo[hpv])):
        if HPV_homo[hpv][i]< 0.35 and HPV_homo_id[i]!=hpv:
            #min_HPV_list.append((HPV_homo_id[i],HPV_homo[hpv][i]))
            min_HPV_list.append(HPV_homo_id[i])
      
      return min_HPV_list
  

  OUTPUT=FILE.replace("_summary_L1.txt","_summary_L2.txt")
  with open (FILE) as f2:
    high_hpv=list()
    low_hpv=list()
    for line in f2:
        line=line.rstrip()
        hpv_count=int(line.split(' ')[0])
        hpv_id=line.split('|')[-2].replace("REF.1","").replace("REF.2","")
        if hpv_count > 99:
          high_hpv.append(hpv_id)

        else:
            low_hpv.append(hpv_id)
    suspected_pool=list()
    confirm_hpv=list()
    suspected_hpv=list()
    confirm_hpv.extend(high_hpv)
    for high in high_hpv:
        suspected_pool.extend(hpv_neighbor(high))
    for low in low_hpv:
        if low in suspected_pool:
            suspected_hpv.append(low)
        else:
            confirm_hpv.append(low)
          
  with open (OUTPUT,'w') as g:
    for c in confirm_hpv:
        newline="confirm"+' '+c+'\n'
        g.write(newline)
    for s in suspected_hpv:
        newline="suspected"+' '+s+'\n'
        g.write(newline)
        
        

def summary_2(F3):
  F2=F3.replace("_hybrid_summary_L3.txt","_summary_L2.txt")  
  F4=F3.replace("_hybrid_summary_L3.txt","_final_summary_L4.txt")
  with open (F2) as f2, open (F3) as f3, open (F4,'w') as g4:
    HPV_f2=dict()
    for line2 in f2:
        line2=line2.rstrip().split(' ')
        HPV_f2[line2[1]]=line2[0]
    for line3 in f3:
        h3=line3.rstrip().split('|')[-2].replace("REF.1","").replace("REF.2","")
        HPV_f2[h3]='confirm'
    for h in HPV_f2:
        if HPV_f2[h] == 'confirm':
            g4.write(h+'\n')
        
        

def HPVname(old):
    if 'REF' in old:
      new=old.split('|')[-2].replace("REF.1","").replace("REF.2","")
    else:
      new=old
    return new

summary_1_input=outprefix+'_summary_L1.txt'
summary_1(summary_1_input)




def hybrid_database (outprefix,min_cov):
    with open (outprefix+'_summary_L2.txt') as hf:
        suspect=0
        for line_hf in hf:
            if 'suspected' in line_hf:
                suspect=suspect+1
    if suspect == 0:
      hybrid_arg0=["cut -d ' ' -f2 " +outprefix+"_summary_L2.txt > "+outprefix+"_final_summary_L4.txt"]
      for command in hybrid_arg0:
         call(command,shell=True)
    else:
        hybrid_arg1=["bedtools bamtofastq -i "+outprefix+".bam -fq "+outprefix+"_HPV.fastq"]
        for command in hybrid_arg1:
          call(command,shell=True)
        hybrid_bowtie=["bowtie2 -x "+database_2+ " --no-unal --quiet -U "+outprefix+"_HPV.fastq -S "+outprefix+"_hybrid.sam -p "+nthreades]
        for command in hybrid_bowtie:
          call(command,shell=True)
     
        hybrid_samtools=['samtools view -bS '+outprefix+'_hybrid.sam'+' | samtools sort  -o '+outprefix+'_hybrid.bam ',\
        "samtools depth "+outprefix+"_hybrid.bam | cut -f1 | uniq -c | rev| cut -d ' ' -f 1,2  | rev | awk -F ' '  '{if ($1 > "+min_cov+") {print $2}}' > "  + outprefix+"_hybrid_summary_L3.txt"]
        for command in hybrid_samtools:
          call(command,shell=True)

        summary_2_input=outprefix+'_hybrid_summary_L3.txt'
        summary_2(summary_2_input)
    


    
    
def quant_HPV(outprefix,min_cov):
    if database_type=="hybrid-mask":
      hybrid_database (outprefix,min_cov)
      f1=open(outprefix+'_final_summary_L4.txt') 
      f2=open(outprefix+'_summary_L1.txt')
      f3=open(HPV_length)

    else:
      f1=open(outprefix+'_L1_id.txt') 
      f2=open(outprefix+'_summary_L1.txt')
      f3=open(HPV_length)

    d2=dict()
    d3=dict()
    for line in f2:
        line=line.rstrip().split(' ')
        d2[HPVname(line[1])]=line[0]
    f2.close()
    for line in f3:
        line=line.rstrip().split(':')
        d3[HPVname(line[0])]=line[1]
    f3.close()
    g=open(outprefix+'_HPV_profile.txt','w')
    headline='HPV_type'+'\t'+'RPK(reads_per_kilobase)'+'\t'+'count_of_reads'+'\n'
    g.write(headline)
    for line in f1:
        line=line.rstrip()
        hpv=HPVname(HPVname(line))
        coverage=str("%.4f" % (1000*float(d2[hpv])/float(d3[hpv])))
        newline=hpv+'\t'+coverage+'\t'+d2[hpv]+'\n'
        g.write(newline)
    f1.close()
    g.close()


quant_HPV(outprefix,min_cov)        
call('mv  '+outprefix+'_HPV_profile.txt '+output_location,shell=True)


