#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 23:08:32 2020

@author: pronaya
@email: pronaya.prosun@gmail.com
"""

import os
import sys
import ntpath
import a_script_rw as pipeline
import time
import functools 
import re
from termcolor import colored

def secToStr(t):
    return ("%dh, %02dm, %02d.%03ds" % functools.reduce(
        lambda ll,b : divmod(ll[0],b) + ll[1:], [(t*1000,),1000,60,60]))
        
start_time = time.time()
#modpath = "/Home/ppd/works/others/B.Napus"
modpath = os.path.dirname(os.path.abspath(__file__))
SEQ_DIR = ""
REF_GENOME = ""
SNPEFF_DB = ""
NO_OF_INTERNAL_THREAD = 2
OVERWRITE_FILE = False
EXIT_ON_ERROR = True

QD_SNP = 0
FS_SNP = 0       
MQ_SNP = 0     
MQRankSum_SNP = 0
ReadPosRankSum_SNP = 0    
SOR_SNP = 0

QD_INDEL = 0       
FS_INDEL = 0          
ReadPosRankSum_INDEL = 0
SOR_INDEL = 0

try:
    def getFloatVal(line, search, defaultVal, check=False):
        matched_key = re.findall(search, line)
        if(len(matched_key)==1 and matched_key[0]==search):
            if(not check):
                print(matched_key, matched_value)
            if(matched_value[0].lstrip('-').replace('.','',1).isdigit()):
                return float(matched_value[0])
            else:
                return defaultVal        
        return "NULL"
    
    with open("pipe.inp", 'r', encoding='utf-8') as f:
        for line in f:
            matched_value = re.findall("\[([^)]*)\]", line)
        
            matched_key = re.findall("SEQUENCE_DIR", line)
            if(len(matched_key)==1 and matched_key[0]=='SEQUENCE_DIR'):
                print(matched_key, matched_value)
                SEQ_DIR = matched_value[0]
        
            matched_key = re.findall("REF_GENOME", line)
            if(len(matched_key)==1 and matched_key[0]=='REF_GENOME'):
                print(matched_key, matched_value)
                REF_GENOME = matched_value[0]
                
            matched_key = re.findall("SNPEFF_DB", line)
            if(len(matched_key)==1 and matched_key[0]=='SNPEFF_DB'):
                print(matched_key, matched_value)
                SNPEFF_DB = matched_value[0]
            
            matched_key = re.findall("NO_OF_INTERNAL_THREAD", line)
            if(len(matched_key)==1 and matched_key[0]=='NO_OF_INTERNAL_THREAD'):
                print(matched_key, matched_value)
                if(matched_value[0].isdigit()):
                    NO_OF_INTERNAL_THREAD = int(matched_value[0])
                else:
                    NO_OF_INTERNAL_THREAD = 2
                    
            matched_key = re.findall("EXIT_ON_ERROR", line)
            if(len(matched_key)==1 and matched_key[0]=='EXIT_ON_ERROR'):
                print(matched_key, matched_value)
                if(matched_value[0].lower()=="t"):
                    EXIT_ON_ERROR = True
                elif(matched_value[0].lower()=="f"):
                    EXIT_ON_ERROR = False
                    
            matched_key = re.findall("OVERWRITE_FILE", line)
            if(len(matched_key)==1 and matched_key[0]=='OVERWRITE_FILE'):
                print(matched_key, matched_value)
                if(matched_value[0].lower()=="t"):
                    OVERWRITE_FILE = True
                elif(matched_value[0].lower()=="f"):
                    OVERWRITE_FILE = False                
            
            if(getFloatVal(line=line, search="QD_SNP", defaultVal=2.0, check=True)!="NULL"):
                QD_SNP = getFloatVal(line=line, search="QD_SNP", defaultVal=2.0)
            
            if(getFloatVal(line=line, search="FS_SNP", defaultVal=60.0, check=True)!="NULL"):
                FS_SNP = getFloatVal(line=line, search="FS_SNP", defaultVal=60.0)
            
            if(getFloatVal(line=line, search="MQ_SNP", defaultVal=40.0, check=True)!="NULL"):
                MQ_SNP = getFloatVal(line=line, search="MQ_SNP", defaultVal=40.0)    
            
            if(getFloatVal(line=line, search="MQRankSum_SNP", defaultVal=-12.5, check=True)!="NULL"):
                MQRankSum_SNP = getFloatVal(line=line, search="MQRankSum_SNP", defaultVal=-12.5)
            
            if(getFloatVal(line=line, search="ReadPosRankSum_SNP", defaultVal=-8.0, check=True)!="NULL"):
                ReadPosRankSum_SNP = getFloatVal(line=line, search="ReadPosRankSum_SNP", defaultVal=-8.0)  
            
            if(getFloatVal(line=line, search="SOR_SNP", defaultVal=4.0, check=True)!="NULL"):    
                SOR_SNP = getFloatVal(line=line, search="SOR_SNP", defaultVal=4.0)
            
            if(getFloatVal(line=line, search="QD_INDEL", defaultVal=2.0, check=True)!="NULL"): 
                QD_INDEL = getFloatVal(line=line, search="QD_INDEL", defaultVal=2.0) 
            
            if(getFloatVal(line=line, search="FS_INDEL", defaultVal=200.0, check=True)!="NULL"):
                FS_INDEL = getFloatVal(line=line, search="FS_INDEL", defaultVal=200.0)
                
            if(getFloatVal(line=line, search="ReadPosRankSum_INDEL", defaultVal=-20.0, check=True)!="NULL"):    
                ReadPosRankSum_INDEL = getFloatVal(line=line, search="ReadPosRankSum_INDEL", defaultVal=0) 
            
            if(getFloatVal(line=line, search="SOR_INDEL", defaultVal=10.0, check=True)!="NULL"):
                SOR_INDEL = getFloatVal(line=line, search="SOR_INDEL", defaultVal=10.0)

except:
    print(colored("Error with pipe.inp file!!!", "red"))
    #exit(1)


VCP = pipeline.VCPipeline(srcPath=modpath, refGenomePath=REF_GENOME, noOfInternalThread=NO_OF_INTERNAL_THREAD)
VCP.setParam(snpEffDb=SNPEFF_DB, overwriteFile=OVERWRITE_FILE, exit_on_error=EXIT_ON_ERROR, qdSnp=QD_SNP, 
             fsSnp=FS_SNP, mqSnp=MQ_SNP, mqRankSumSnp=MQRankSum_SNP, readPosRankSumSnp=ReadPosRankSum_SNP, 
             sorSnp=SOR_SNP, qdIndel=QD_INDEL, fsIndel=FS_INDEL, readPosRankSumIndel=ReadPosRankSum_INDEL, 
             sorIndel=SOR_INDEL)
VCP.checkRequiredTools()


for dirInfo in os.walk(SEQ_DIR):
    seq_time = time.time()
    path, dirs, files = dirInfo
    matchingFiles = [s for s in files if ((".gz" in s) or (".fq" in s))]
    if(len(matchingFiles)>=2):
        selFiles = [s for s in matchingFiles if s.endswith(".fq" )]
        if(len(selFiles)<2):
            selFiles = [s for s in matchingFiles if s.endswith(".gz" )]
        if(len(selFiles)>=2):
            print("Accessing sequence directory:", ntpath.basename(path))
            selFiles.sort()
            SEQ1 = os.path.join(path, selFiles[0])
            SEQ2 = os.path.join(path, selFiles[1])
            print(SEQ1)
            print(SEQ2)
            VCP.setSeqFiles(seq1=SEQ1, seq2=SEQ2)
            VCP.printAllParam()
            VCP.executePipeline()
            print("Total time (per seqs):", secToStr(time.time()-seq_time), "\n") 
            sys.stdout.flush()

print("Total time:", secToStr(time.time()-start_time), "\n")   
