#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 21:08:05 2020

@author: pronaya
@email: pronaya.prosun@gmail.com
"""

import sys
import os
import subprocess
import threading
import time
import re
from termcolor import colored
import ntpath
from sys import exit

class VCPipeline:
    def __init__(self, srcPath="", refGenomePath="", seq1="", seq2="", snpEffDb="",
                 noOfInternalThread=40, heapSizeGb=2, overwriteFile=False, exit_on_error=True, 
                 qdSnp=2.0, fsSnp=60.0, mqSnp=40.0, mqRankSumSnp=-12.5, readPosRankSumSnp=-8.0, 
                 sorSnp=4.0, qdIndel=2.0, fsIndel=200.0, readPosRankSumIndel=-20.0, 
                 sorIndel=10.0):
        self.NO_OF_INTERNAL_THREAD = noOfInternalThread
        self.AMOUNT_OF_HEAP = heapSizeGb #in GB 
        self.OVER_WRITE_FILE = overwriteFile
        self.EXIT_ON_ERROR = exit_on_error

        ###Checking all the necessary tools for running pipeline
        #self.modpath = '/media/pronaya/My Passport/RRRRR/B.Napus'
        self.modpath = srcPath
        #self.modpath = os.path.dirname(os.path.abspath(__file__))
        self.toolspath = os.path.join(self.modpath, "tools")
        self.JAVA = os.path.join(self.toolspath, "java8jre/bin/java")
        self.PICARD_JAR = os.path.join(self.toolspath, "picard.jar")
        self.GATK_JAR = os.path.join(self.toolspath, "GenomeAnalysisTK.jar")
        self.SNPEFF_JAR = os.path.join(self.toolspath, "snpEff_v4_1i_core/snpEff/snpEff.jar")
        #self.SNPEFF_DB = os.path.join(self.toolspath, "snpEffectPredictor.bin")
        self.SNPEFF_DB = snpEffDb

        #Reference Genome
        #self.REF_GENOME = os.path.join(self.modpath, "refGenom/_Brassica_napus.AST_PRJEB5043_v1.dna.toplevel.fa")
        self.REF_GENOME = refGenomePath
        #self.SEQ1 = os.path.join(self.modpath, "SeqFiles/CRR045427-/CRR045427_f1.fq.gz")
        #self.SEQ2 = os.path.join(self.modpath, "SeqFiles/CRR045427-/CRR045427_r2.fq.gz")
        self.SEQ1 = seq1
        self.SEQ2 = seq2
        self.OUTPUT = os.path.join(ntpath.dirname(self.SEQ1), "results")
        self.ID = ntpath.basename(ntpath.dirname(self.SEQ1))
        self.OUTPUT_STAT = ntpath.dirname(ntpath.dirname(self.SEQ1))

        #Parameter for filtering SNPs
        self.QD_SNP = qdSnp
        self.FS_SNP = fsSnp
        self.MQ_SNP = mqSnp
        self.MQRankSum_SNP = mqRankSumSnp
        self.ReadPosRankSum_SNP = readPosRankSumSnp
        self.SOR_SNP = sorSnp

        #Parameter for filtering Indels
        self.QD_INDEL = qdIndel
        self.FS_INDEL = fsIndel
        self.ReadPosRankSum_INDEL = readPosRankSumIndel
        self.SOR_INDEL = sorIndel
        
        self.tool_not_installed = []
        
        #self.checkRequiredTools()
        
        #self.executePipeline()
    def printAllParam(self):
        print(colored("##########################################################", "yellow"))
        print("Printing all parameters...")
        print(" NO_OF_INTERNAL_THREAD:", self.NO_OF_INTERNAL_THREAD, "\n",
              "AMOUNT_OF_HEAP:", self.AMOUNT_OF_HEAP, "\n", 
              "OVER_WRITE_FILE: ", self.OVER_WRITE_FILE, "\n",
              "EXIT_ON_ERROR:", self.EXIT_ON_ERROR, "\n", 
              "SRC_PATH:", self.modpath, "\n", 
              "toolspath:", self.toolspath, "\n", 
              "JAVA:", self.JAVA, "\n" 
              "PICARD_JAR:", self.PICARD_JAR, "\n",
              "GATK_JAR:", self.GATK_JAR, "\n", 
              "SNPEFF_JAR:", self.SNPEFF_JAR, "\n", 
              "SNPEFF_DB:", self.SNPEFF_DB, "\n", 
              "REF_GENOME:", self.REF_GENOME, "\n", 
              "SEQ1:", self.SEQ1, "\n",
              "SEQ2:", self.SEQ2, "\n", 
              "OUTPUT:", self.OUTPUT, "\n", 
              "ID:", self.ID, "\n", 
              "OUTPUT_STAT:", self.OUTPUT_STAT, "\n", 
              "QD_SNP:", self.QD_SNP, "\n", 
              "FS_SNP:", self.FS_SNP, "\n", 
              "MQ_SNP:", self.MQ_SNP, "\n", 
              "MQRankSum_SNP:", self.MQRankSum_SNP, "\n", 
              "ReadPosRankSum_SNP:", self.ReadPosRankSum_SNP, "\n", 
              "SOR_SNP:", self.SOR_SNP, "\n", 
              "QD_INDEL:", self.QD_INDEL, "\n", 
              "FS_INDEL:", self.FS_INDEL, "\n", 
              "ReadPosRankSum_INDEL:", self.ReadPosRankSum_INDEL, "\n", 
              "SOR_INDEL:", self.SOR_INDEL, "\n")
        print(colored("##########################################################", "yellow"))
        
    def setParam(self, snpEffDb="", overwriteFile=False, exit_on_error=True, qdSnp=2.0, fsSnp=60.0, mqSnp=40.0, 
                 mqRankSumSnp=-12.5, readPosRankSumSnp=-8.0, sorSnp=4.0, qdIndel=2.0, fsIndel=200.0, 
                 readPosRankSumIndel=-20.0, sorIndel=10.0):
        self.SNPEFF_DB = snpEffDb
        self.OVER_WRITE_FILE = overwriteFile
        self.EXIT_ON_ERROR = exit_on_error
        #Parameter for filtering SNPs
        self.QD_SNP = qdSnp
        self.FS_SNP = fsSnp
        self.MQ_SNP = mqSnp
        self.MQRankSum_SNP = mqRankSumSnp
        self.ReadPosRankSum_SNP = readPosRankSumSnp
        self.SOR_SNP = sorSnp
        #Parameter for filtering Indels
        self.QD_INDEL = qdIndel
        self.FS_INDEL = fsIndel
        self.ReadPosRankSum_INDEL = readPosRankSumIndel
        self.SOR_INDEL = sorIndel
        
    def setSeqFiles(self, seq1, seq2):
        self.SEQ1 = seq1
        self.SEQ2 = seq2
        self.OUTPUT = os.path.join(ntpath.dirname(self.SEQ1), "results")
        self.ID = ntpath.basename(ntpath.dirname(self.SEQ1))
        self.OUTPUT_STAT = ntpath.dirname(ntpath.dirname(self.SEQ1))

            
    def executePipeline(self):
        self.decompressFiles(self.EXIT_ON_ERROR, self.OVER_WRITE_FILE)
        #Step 0
        self.indexRefGenome(self.EXIT_ON_ERROR, self.OVER_WRITE_FILE)
        #Step 1: Alignment – Map to Reference
        self.generateAlignedRead(self.EXIT_ON_ERROR, self.OVER_WRITE_FILE)
        #Step 2: Sort SAM file by coordinate, convert to BAM
        self.sortSAMandConvertToBAM(self.EXIT_ON_ERROR, self.OVER_WRITE_FILE)
        #Step 3: Collect Alignment & Insert Size Metrics
        self.collectAlignmentAndInsertMatrix(self.EXIT_ON_ERROR, self.OVER_WRITE_FILE)
        #Step 4: Mark Duplicates
        self.markDuplicates(self.EXIT_ON_ERROR, self.OVER_WRITE_FILE)
        #Step 5: Build BAM Index
        self.buildBamIndex(self.EXIT_ON_ERROR, self.OVER_WRITE_FILE)
        #Step 6: Create Realignment Targets
        self.createRealignmentTargets(self.EXIT_ON_ERROR, self.OVER_WRITE_FILE)
        #Step 7: Realign Indels
        self.realignIndels(self.EXIT_ON_ERROR, self.OVER_WRITE_FILE)
        #Step 8: Call Variants
        self.OUTPUT_RAW_VCF = self.variantsCalling(inp=self.OUTPUT_REALIGNED_BAM, 
                        exit_on_error=self.EXIT_ON_ERROR, overwrite=self.OVER_WRITE_FILE)
        #Step 9: Extract SNPs and Indels
        self.OUTPUT_RAW_SNPS_VCF, self.OUTPUT_RAW_INDELS_VCF = self.extractSNPsAndIndel(
                        inp=self.OUTPUT_RAW_VCF, exit_on_error=self.EXIT_ON_ERROR, 
                        overwrite=self.OVER_WRITE_FILE)
        #Step 10: Filter SNPs
        self.OUTPUT_FILTERED_SNPS_VCF = self.filterSNPs(inp=self.OUTPUT_RAW_SNPS_VCF, 
                        exit_on_error=self.EXIT_ON_ERROR, overwrite=self.OVER_WRITE_FILE)
        #Step 11: Filter Indels
        self.OUTPUT_FILTERED_INDELS_VCF = self.filterIndels(inp=self.OUTPUT_RAW_INDELS_VCF, 
                        exit_on_error=self.EXIT_ON_ERROR, overwrite=self.OVER_WRITE_FILE)
        #Step 12: Base Quality Score Recalibration (BQSR) 1
        self.scoreRecalibration(self.EXIT_ON_ERROR, self.OVER_WRITE_FILE)
        #Step 13: Base Quality Score Recalibration (BQSR) 2
        self.scoreRecalibration2(self.EXIT_ON_ERROR, self.OVER_WRITE_FILE)
        #Step 14: Analyze Covariates  
        self.analyzeCovariates(self.EXIT_ON_ERROR, self.OVER_WRITE_FILE)
        #Step 15: Apply BQSR
        self.applyBQSR(self.EXIT_ON_ERROR, self.OVER_WRITE_FILE)   
        #Step 16: Call Variants
        self.OUTPUT_RAW_RECAL_VCF = self.variantsCalling(inp=self.OUTPUT_RECAL_BAM, 
                        out="raw_variants_recal.vcf", re1="Step 16", 
                        exit_on_error=self.EXIT_ON_ERROR, overwrite=self.OVER_WRITE_FILE)
        #Step 17: Extract SNPs & Indels
        self.OUTPUT_RAW_SNPS_RECAL_VCF, self.OUTPUT_RAW_INDELS_RECAL_VCF = self.extractSNPsAndIndel(
                        inp=self.OUTPUT_RAW_RECAL_VCF, out1="raw_snps_recal.vcf", 
                        out2="raw_indels_recal.vcf", re1="Step 17", 
                        exit_on_error=self.EXIT_ON_ERROR, overwrite=self.OVER_WRITE_FILE)
        #Step 18: Filter SNPs
        self.OUTPUT_FILTERED_SNPS_FINAL_VCF = self.filterSNPs(inp=self.OUTPUT_RAW_SNPS_RECAL_VCF, 
                        out="filtered_snps_final.vcf", re1="Step 18", 
                        exit_on_error=self.EXIT_ON_ERROR, overwrite=self.OVER_WRITE_FILE)
        #Step 19: 	Filter Indels
        self.OUTPUT_FILTERED_INDELS_FINAL_VCF = self.filterIndels(
                        inp=self.OUTPUT_RAW_INDELS_RECAL_VCF, out="filtered_indels_recal.vcf",
                        re1="Step 19", exit_on_error=self.EXIT_ON_ERROR, 
                        overwrite=self.OVER_WRITE_FILE)
        #Step 20: Annotate SNPs and Predict Effects
        self.annotateSNPsAndPredictEffects(False, self.OVER_WRITE_FILE)
        #Step 21: Annotate SNPs and Predict Effects
        #self.computeCoverageStat(self.EXIT_ON_ERROR, self.OVER_WRITE_FILE)
        #Step 22
        #self.computeStat(self.EXIT_ON_ERROR, self.OVER_WRITE_FILE)
        
    def checkRequiredTools(self):
        print("Checking all the necessary tools for running pipeline...")
        self.checkJAVA()
        self.checkR()
        self.checkOtherTools(check="bwa", remark="BWA")
        self.checkOtherTools(check="samtools", remark="Samtools")
        #self.checkOtherTools(check="bedtools", remark="Bedtools")
        self.checkJARs(JAR_FILE=self.PICARD_JAR, remark="Picard")
        self.checkJARs(JAR_FILE=self.GATK_JAR, remark="GATK3")
        
        if(len(self.tool_not_installed)==1):
            print("Following tool is not available.")
            print(self.tool_not_installed)
            if(self.EXIT_ON_ERROR):
                exit(1)
        elif(len(self.tool_not_installed)>1):
            print("Following tools are not available.")
            print(self.tool_not_installed)
            if(self.EXIT_ON_ERROR):
                exit(1)
        
        print(colored("##########################################################", "yellow"))
           
    #checking java installation and its version
    def checkJAVA(self):
        command = "which java"
        try:
            res = subprocess.check_output(command, shell=True)
            res = (re.sub("\\\\n'", "", str(res))).split("/")
            if('java'==res[-1]):
                java_version = subprocess.check_output(['java', '-version'], stderr=subprocess.STDOUT)
                java_version = float(re.search('\"(\d+\.\d+).*\"', str(java_version)).groups()[0])
                
                if(java_version>=9.0):
                    print("Your system has java", java_version, 
                          "installed. Unfortunately GATK3 also needs java 8.")
                    if(os.path.isfile(self.JAVA)):
                        print("Java 8:", colored("Found", "green"), "(in tools)")
                        print("path:", self.JAVA)
                    else:
                        print("Java 8 is not available in the tools folder either.")
                        self.tool_not_installed.append("java 8")
                else:
                    self.JAVA = "java"
                    print("Java ", str(java_version), ": ", colored("Found", "green"), " (in system)", sep="")
                    print("path:", self.JAVA)
            else:
                print("Java is not available in the system.")
                self.tool_not_installed.append("java in system")
                    
        except:
            if(os.path.isfile(self.JAVA)):
                print("Java 8:", colored("Found", "green"), "(in tools)")
                print("path:", self.JAVA)
            else:
                print("Java 8 is not available in the tools folder either.")
                self.tool_not_installed.append("java 8")
        sys.stdout.flush()

    #R
    def checkR(self):
        command = "which R"
        try:
            res = subprocess.check_output(command, shell=True)
            res = (re.sub("\\\\n'", "", str(res))).split("/")
            if('R'==res[-1]):
                print("R:", colored("Installed", "green"), "(in system)")
                print("Checking R packages, ")
                pkList = ["ggplot2", "gplots", "reshape", "grid", "tools", "gsalib"]
                for pk in pkList:
                    resRPk = subprocess.check_output("Rscript --vanilla "+ 
                            "%r"%os.path.join(self.toolspath, "checkRPackage.R")+" "+pk, shell=True)
                    resRPk = (re.sub("\\\\n'", "", str(resRPk))).split(" ")
                    if(resRPk[-1]=="1"):
                        print("R:", colored('"'+pk+'" installed', "green"))
                    else:
                        print("R:", colored('"'+pk+'" not installed', "red"))
                        self.tool_not_installed.append("R "+pk)
            else:
                print("R:", colored("Not installed", "red"))
                self.tool_not_installed.append("R")
        except:
            print("R:", colored("Not installed", "red"))
            self.tool_not_installed.append("R/R package(s)")
        sys.stdout.flush()

    #BWA, Samtools and Bedtools
    def checkOtherTools(self, check, remark="Name of the tool"):
        command = "which "+check
        try:
            res = subprocess.check_output(command, shell=True)
            res = (re.sub("\\\\n'", "", str(res))).split("/")
            if(check==res[-1]):
                print(remark+":", colored("Installed", "green"), "(in system)")
            else:
                print(remark+":", colored("Not installed", "red"))
                self.tool_not_installed.append(remark)
        except:
            print(remark+":", colored("Not installed", "red"))
            self.tool_not_installed.append(remark)
        sys.stdout.flush()
    
    #Picard and GATK
    def checkJARs(self, JAR_FILE, remark="Name of the JAR"):
        try:
            if(os.path.isfile(JAR_FILE)):
                print(remark+" (jar):", colored("Found", "green"), "(in tools)")
                print("path:", JAR_FILE)
            else:
                print(remark+" (jar) is", colored("not available", "red"), "in the tools folder.")
                self.tool_not_installed.append(remark)
        except:
            print(remark+" (jar) is", colored("not available", "red"), "in the tools folder.")
            self.tool_not_installed.append(remark)
        sys.stdout.flush()
    
    #Extract files
    def decompressFiles(self, exit_on_error=True, overwrite=False):
        if not (os.path.isfile(self.REF_GENOME)):
            print("Reference genome", colored("not found", "red"))
            if(exit_on_error):
                exit(1)
        if not (os.path.isfile(self.SEQ1)):
            print("Sequence 1", colored("not found", "red"))
            if(exit_on_error):
                exit(1)
        if not (os.path.isfile(self.SEQ2)):
            print("Sequence 2", colored("not found", "red"))
            if(exit_on_error):
                exit(1)
                
        #internal function
        def checkIfCompressed(file):
            filename, file_extension = os.path.splitext(file)
            #print(filename, file_extension)
            if(file_extension==".gz"):
                if(overwrite==False and os.path.isfile(filename)):
                    print(filename, colored("Decompressed file already exists!", "yellow"))
                    file = filename
                    return file
                try:
                    print(file, colored("Decompressing...", "yellow"))
                    command = "gunzip -k -v '"+file+"'"
                    subprocess.check_output(command, shell=True)
                    print(file, colored("Decompressed!", "green"))
                    file = filename
                except:
                    print(colored("Error occured during decompression!", "red"))
                    if(exit_on_error):
                        exit(1)
            sys.stdout.flush()
            return file
        self.REF_GENOME = checkIfCompressed(self.REF_GENOME)
        self.SEQ1 = checkIfCompressed(self.SEQ1)
        self.SEQ2 = checkIfCompressed(self.SEQ2)

    #Step 0
    def indexRefGenome(self, exit_on_error=True, overwrite=False):   
        print(colored("##########################################################", "yellow"), "\n")
        print("Executing Step 0")
        try:
            print(colored("Indexing reference genome...", "yellow"))
            sys.stdout.flush()
            if(overwrite==False):
                given_file_list_ext = [".amb", ".ann", ".bwt", ".pac", ".sa"]
                head, tail = os.path.split(self.REF_GENOME)
                files_in_dir = [f for f in os.listdir(head) if os.path.isfile(os.path.join(head, f))]
                given_file_list = [tail+f for f in given_file_list_ext]
                set_dir_files = set(files_in_dir)
                set_given_file_list = set(given_file_list)
                
                if not (set_given_file_list.issubset(set_dir_files)):
                    command = "bwa index '" + self.REF_GENOME + "'"
                    subprocess.check_output(command, shell=True)
                else:
                    print("File already exists: \n", given_file_list)
                
                #checking .fai file
                if not (set([tail+".fai"]).issubset(set_dir_files)):
                    command = "samtools faidx '" + self.REF_GENOME + "'"
                    subprocess.check_output(command, shell=True)
                else:
                    print("File already exists: ", tail+".fai")
                
                #checking .dict file
                filename, file_extension = os.path.splitext(self.REF_GENOME)
                head, tail = os.path.split(filename)
                if not (set([tail+".dict"]).issubset(set_dir_files)):
                #java -XX:ParallelGCThreads=40 -jar picard.jar CreateSequenceDictionary REFERENCE=refGenom/refGenom.fa OUTPUT=refGenom/refGenom.dict
                    command = "%r"%self.JAVA+" -XX:ParallelGCThreads="+str(self.NO_OF_INTERNAL_THREAD)+ \
                              " -jar "+"%r"%self.PICARD_JAR+" CreateSequenceDictionary REFERENCE="+ \
                              "%r"%self.REF_GENOME+" OUTPUT="+"%r"%os.path.join(filename+".dict") 
                    subprocess.check_output(command, shell=True)
                else:
                    print("File already exists: ", tail+".dict")
            else:
                command = "bwa index '" + self.REF_GENOME + "'"
                subprocess.check_output(command, shell=True)
                command = "samtools faidx '" + self.REF_GENOME + "'"
                subprocess.check_output(command, shell=True)
                filename, file_extension = os.path.splitext(self.REF_GENOME)
                command = "%r"%self.JAVA+" -XX:ParallelGCThreads="+str(self.NO_OF_INTERNAL_THREAD)+ \
                          " -jar "+"%r"%self.PICARD_JAR+" CreateSequenceDictionary REFERENCE="+ \
                          "%r"%self.REF_GENOME+" OUTPUT="+"%r"%os.path.join(filename+".dict") 
                subprocess.check_output(command, shell=True)
            print(colored("Step 0 execution finished!", "green"))
            print(colored("##########################################################", "yellow"), "\n")
        except:
            print(colored("Step 0 execution halted. Error occured!", "red"))
            if(exit_on_error):
                exit(1)
        sys.stdout.flush()
    
    #Common function
    def executeCommand(self, out, command, remark1, remark2, remark3, exit_on_exception=True, overwrite=False):
        sys.stdout.flush()
        out_main = ""
        if(isinstance(out, list)):
            out_main = out[0]
        else:
            out_main = out
        try: 
            if(overwrite==False and os.path.isfile(out_main)):
                print(out_main, colored("file already exists!", "yellow"))
            else:
                print("\nGenerating", remark2, "file(s). Please wait...")
                sys.stdout.flush()
                subprocess.check_output(command, shell=True)
                print(colored(remark3+" file(s) generated!", "green"), out)
            print(colored(remark1+" execution finished!", "green"))
            print(colored("##########################################################", "yellow"), "\n")
        except:
            try:
                if(isinstance(out, list)):
                    for i in range(len(out)):
                        if os.path.exists(out[i]):
                            os.remove(out[i])
                else:
                    if os.path.exists(out):
                        os.remove(out)
            except:
                print("")         
            print(colored(remark1+" execution halted. Error occured!", "red"))
            if(exit_on_exception):
                exit(1)
        sys.stdout.flush()

    #Step 1: Alignment – Map to Reference
    #Creating output directory if needed
    def generateAlignedRead(self, exit_on_error=True, overwrite=False):
        if not os.path.exists(self.OUTPUT):
            os.mkdir(self.OUTPUT)
        print("Executing Step 1")
        print(colored("Alignment - Map to Reference...", "yellow"))
        print("Following files are needed: ", self.REF_GENOME, self.SEQ1, self.SEQ2, sep="\n")
        self.OUTPUT_SAM = os.path.join(self.OUTPUT, "aligned_reads.sam")
        print("Output:", colored(self.OUTPUT_SAM, "yellow"))
        command = "bwa mem -t "+str(self.NO_OF_INTERNAL_THREAD)+" -M -R "+ \
                  r"'@RG\tID:sample_1\tLB:sample_1\tPL:ILLUMINA\tPM:HISEQ\tSM:sample_1'"+ \
                  " "+"%r"%self.REF_GENOME+" "+"%r"%self.SEQ1+" "+"%r"%self.SEQ2+" > "+"%r"%self.OUTPUT_SAM
        self.executeCommand(self.OUTPUT_SAM, command, remark1="Step 1", remark2="SAM", 
                            remark3="SAM", exit_on_exception=exit_on_error, overwrite=overwrite)
        sys.stdout.flush()
        
    #Step 2: Sort SAM file by coordinate, convert to BAM
    def sortSAMandConvertToBAM(self, exit_on_error=True, overwrite=False): 
        print("Executing Step 2")
        print(colored("Sort SAM file by coordinate, convert to BAM...", "yellow"))
        print("Following file is needed: ", self.OUTPUT_SAM, sep="\n")
        self.OUTPUT_BAM = os.path.join(self.OUTPUT, "sorted_reads.bam")
        print("Output:", colored(self.OUTPUT_BAM, "yellow"))
        command = "java"+" -XX:ParallelGCThreads="+str(self.NO_OF_INTERNAL_THREAD)+" -jar "+ \
                  "%r"%self.PICARD_JAR+" SortSam INPUT="+"%r"%self.OUTPUT_SAM+" OUTPUT="+ \
                  "%r"%self.OUTPUT_BAM+" SORT_ORDER=coordinate"
        self.executeCommand(self.OUTPUT_BAM, command, remark1="Step 2", remark2="BAM (Sorted SAM)", 
                            remark3="BAM", exit_on_exception=exit_on_error, overwrite=overwrite)
        sys.stdout.flush()
    
    #Step 3: Collect Alignment & Insert Size Metrics
    def collectAlignmentAndInsertMatrix(self, exit_on_error=True, overwrite=False):     
        print("Executing Step 3")
        print(colored("Collect Alignment and Insert Size Metrics...", "yellow"))
        print("Following files are needed: ", self.REF_GENOME, self.OUTPUT_BAM, sep="\n")
        self.OUTPUT_ALIGN_MAT_TXT = os.path.join(self.OUTPUT, "alignment_metrics.txt")
        self.OUTPUT_INSERT_MAT = os.path.join(self.OUTPUT, "insert_metrics.txt")
        self.OUTPUT_INSERT_SZ_HIST_PDF = os.path.join(self.OUTPUT, "insert_size_histogram.pdf")
        self.OUTPUT_DEPTH_TXT = os.path.join(self.OUTPUT, "depth_out.txt")
        print("Output:", colored(self.OUTPUT_ALIGN_MAT_TXT, "yellow"), colored(self.OUTPUT_INSERT_MAT, "yellow"),
              colored(self.OUTPUT_INSERT_SZ_HIST_PDF, "yellow"), colored(self.OUTPUT_DEPTH_TXT, "yellow"), sep="\n")
        command ="java"+" -XX:ParallelGCThreads="+str(self.NO_OF_INTERNAL_THREAD)+" -jar "+ \
                 "%r"%self.PICARD_JAR+" CollectAlignmentSummaryMetrics R="+"%r"%self.REF_GENOME+ \
                 " I="+"%r"%self.OUTPUT_BAM+" O="+"%r"%self.OUTPUT_ALIGN_MAT_TXT
        self.executeCommand(self.OUTPUT_ALIGN_MAT_TXT, command, remark1="Step 3 part 1/3", 
                            remark2="alignment metrics (TXT)", remark3="Alignment metrics (TXT)", 
                            exit_on_exception=exit_on_error, overwrite=overwrite)
        command = "java"+" -XX:ParallelGCThreads="+str(self.NO_OF_INTERNAL_THREAD)+" -jar "+ \
                  "%r"%self.PICARD_JAR+" CollectInsertSizeMetrics INPUT="+"%r"%self.OUTPUT_BAM+ \
                  " OUTPUT="+"%r"%self.OUTPUT_INSERT_MAT+" HISTOGRAM_FILE="+"%r"%self.OUTPUT_INSERT_SZ_HIST_PDF
        self.executeCommand([self.OUTPUT_INSERT_MAT, self.OUTPUT_INSERT_SZ_HIST_PDF], command, 
                            remark1="Step 3 part 2/3", remark2="insert matrix (TXT) & histogram (PDF)", 
                            remark3="TXT and PDF", exit_on_exception=exit_on_error, overwrite=overwrite)
        command = "samtools depth -a "+"%r"%self.OUTPUT_BAM+" > "+"%r"%self.OUTPUT_DEPTH_TXT
        self.executeCommand(self.OUTPUT_DEPTH_TXT, command, remark1="Step 3 part 3/3", remark2="depth (TXT)", 
                       remark3="Depth (TXT)", exit_on_exception=exit_on_error, overwrite=overwrite)
        sys.stdout.flush()
    
    #Step 4: Mark Duplicates
    def markDuplicates(self, exit_on_error=True, overwrite=False): 
        print("Executing Step 4")
        print(colored("Mark Duplicates...", "yellow"))
        print("Following file is needed: ", self.OUTPUT_BAM, sep="\n")
        self.OUTPUT_DEDUP_BAM = os.path.join(self.OUTPUT, "dedup_reads.bam")
        self.OUTPUT_MAT_TXT = os.path.join(self.OUTPUT, "metrics.txt")
        print("Output:", colored(self.OUTPUT_DEDUP_BAM, "yellow"), colored(self.OUTPUT_MAT_TXT, "yellow"), sep="\n")
        command = "java"+" -XX:ParallelGCThreads="+str(self.NO_OF_INTERNAL_THREAD)+" -jar "+ \
                 "%r"%self.PICARD_JAR+" MarkDuplicates INPUT="+"%r"%self.OUTPUT_BAM+" OUTPUT="+ \
                 "%r"%self.OUTPUT_DEDUP_BAM+" METRICS_FILE="+"%r"%self.OUTPUT_MAT_TXT
        self.executeCommand([self.OUTPUT_DEDUP_BAM, self.OUTPUT_MAT_TXT], command, remark1="Step 4", 
                            remark2="BAM (duplicates marked)", remark3="DEDUP BAM", 
                            exit_on_exception=exit_on_error, overwrite=overwrite)
        sys.stdout.flush()
    
    #Step 5: Build BAM Index
    def buildBamIndex(self, exit_on_error=True, overwrite=False): 
        print("Executing Step 5")
        print(colored("Building BAM Index...", "yellow"))
        print("Following file is needed: ", self.OUTPUT_DEDUP_BAM, sep="\n")
        self.OUTPUT_DEDUP_BAI = os.path.join(self.OUTPUT, "dedup_reads.bai")
        print("Output:", colored(self.OUTPUT_DEDUP_BAI, "yellow"))
        command = "java"+" -XX:ParallelGCThreads="+str(self.NO_OF_INTERNAL_THREAD)+" -jar "+ \
                  "%r"%self.PICARD_JAR+" BuildBamIndex INPUT="+"%r"%self.OUTPUT_DEDUP_BAM
        self.executeCommand(self.OUTPUT_DEDUP_BAI, command, remark1="Step 5", remark2="BAI (BAM Index)", 
                       remark3="BAI", exit_on_exception=exit_on_error, overwrite=overwrite)
        sys.stdout.flush()
    
    #Step 6: Create Realignment Targets
    def createRealignmentTargets(self, exit_on_error=True, overwrite=False): 
        print("Executing Step 6")
        print(colored("Creating Realignment Targets...", "yellow"))
        print("Following files are needed: ", self.REF_GENOME, self.OUTPUT_DEDUP_BAM, sep="\n")
        self.OUTPUT_REALIGNMENT_LIST = os.path.join(self.OUTPUT, "realignment_targets.list")
        print("Output:", colored(self.OUTPUT_REALIGNMENT_LIST, "yellow"))
        command = "%r"%self.JAVA+" -XX:ParallelGCThreads="+str(self.NO_OF_INTERNAL_THREAD)+ \
                  " -jar "+"%r"%self.GATK_JAR+" -T RealignerTargetCreator -R "+"%r"%self.REF_GENOME+ \
                  " -I "+"%r"%self.OUTPUT_DEDUP_BAM+" -o "+"%r"%self.OUTPUT_REALIGNMENT_LIST
        self.executeCommand(self.OUTPUT_REALIGNMENT_LIST, command, remark1="Step 6", 
                            remark2="list (Realignment)", remark3="List", 
                            exit_on_exception=exit_on_error, overwrite=overwrite)
        sys.stdout.flush()

    #Step 7: Realign Indels
    def realignIndels(self, exit_on_error=True, overwrite=False): 
        print("Executing Step 7")
        print(colored("Realigning Indels...", "yellow"))
        print("Following files are needed: ", self.REF_GENOME, self.OUTPUT_DEDUP_BAM, 
              self.OUTPUT_REALIGNMENT_LIST, sep="\n")
        self.OUTPUT_REALIGNED_BAM = os.path.join(self.OUTPUT, "realigned_reads.bam")
        print("Output:", colored(self.OUTPUT_REALIGNED_BAM, "yellow"))
        command = "%r"%self.JAVA+" -XX:ParallelGCThreads="+str(self.NO_OF_INTERNAL_THREAD)+ \
                  " -jar "+"%r"%self.GATK_JAR+" -T IndelRealigner -R "+"%r"%self.REF_GENOME+ \
                  " -I "+"%r"%self.OUTPUT_DEDUP_BAM+" -targetIntervals "+ \
                  "%r"%self.OUTPUT_REALIGNMENT_LIST+" -o "+"%r"%self.OUTPUT_REALIGNED_BAM
        self.executeCommand(self.OUTPUT_REALIGNED_BAM, command, remark1="Step 7", 
                            remark2="BAM (Realigned reads)", remark3="BAM", 
                            exit_on_exception=exit_on_error, overwrite=overwrite)
        sys.stdout.flush()

    #Step 8: Call Variants
    def variantsCalling(self, inp, out="raw_variants.vcf", re1="Step 8", showCommand=False, 
                        exit_on_error=True, overwrite=False): 
        print("Executing "+re1)
        print(colored("Variants calling...", "yellow"))
        print("Following files are needed: ", self.REF_GENOME, inp, sep="\n")
        OUTPUT_RAW_VCF = os.path.join(self.OUTPUT, out)
        print("Output:", colored(OUTPUT_RAW_VCF, "yellow"))
        command = "%r"%self.JAVA+" -XX:ParallelGCThreads="+str(self.NO_OF_INTERNAL_THREAD)+ \
                  " -Xmx"+str(self.AMOUNT_OF_HEAP)+"g -jar "+"%r"%self.GATK_JAR+ \
                  " -T HaplotypeCaller -nct "+str(self.NO_OF_INTERNAL_THREAD)+" -R "+ \
                  "%r"%self.REF_GENOME+" -I "+"%r"%inp+" -o "+"%r"%OUTPUT_RAW_VCF
        if(showCommand):
            print("Command: \n",command)
        self.executeCommand(OUTPUT_RAW_VCF, command, remark1=re1, 
                            remark2="VCF (raw variants)", remark3="VCF (raw variants)", 
                            exit_on_exception=exit_on_error, overwrite=overwrite)
        sys.stdout.flush()
        return OUTPUT_RAW_VCF

    #Step 9: Extract SNPs and Indels
    def extractSNPsAndIndel(self, inp, out1="raw_snps.vcf", out2="raw_indels.vcf", re1="Step 9", 
                            showCommand=False, exit_on_error=True, overwrite=False): 
        print("Executing "+re1)
        print(colored("Extracting SNPs...", "yellow"))
        print("Following files are needed: ", self.REF_GENOME, inp, sep="\n")
        OUTPUT_RAW_SNPS_VCF = os.path.join(self.OUTPUT, out1)
        print("Output:", colored(OUTPUT_RAW_SNPS_VCF, "yellow"))
        command = "%r"%self.JAVA+" -XX:ParallelGCThreads="+str(self.NO_OF_INTERNAL_THREAD)+ \
                  " -jar "+"%r"%self.GATK_JAR+" -T SelectVariants "+" -R "+"%r"%self.REF_GENOME+ \
                  " -V "+"%r"%inp+" -selectType SNP -o "+"%r"%OUTPUT_RAW_SNPS_VCF
        if(showCommand):
            print("Command: \n",command)
        self.executeCommand(OUTPUT_RAW_SNPS_VCF, command, remark1=re1+" part 1/2", 
                            remark2="VCF (raw SNPS)", remark3="VCF (raw SNPS)", 
                            exit_on_exception=exit_on_error, overwrite=overwrite)
        
        print(colored("Extracting Indels...", "yellow"))
        print("Following files are needed: ", self.REF_GENOME, inp, sep="\n")
        OUTPUT_RAW_INDELS_VCF = os.path.join(self.OUTPUT, out2)
        print("Output:", colored(OUTPUT_RAW_INDELS_VCF, "yellow"))
        command = "%r"%self.JAVA+" -XX:ParallelGCThreads="+str(self.NO_OF_INTERNAL_THREAD)+ \
                  " -jar "+"%r"%self.GATK_JAR+" -T SelectVariants "+" -R "+"%r"%self.REF_GENOME+ \
                  " -V "+"%r"%inp+" -selectType INDEL -o "+"%r"%OUTPUT_RAW_INDELS_VCF
        if(showCommand):
            print("Command: \n",command)
        self.executeCommand(OUTPUT_RAW_INDELS_VCF, command, remark1=re1+" part 2/2", 
                            remark2="VCF (raw INDELS)", remark3="VCF (raw INDELS)", 
                            exit_on_exception=exit_on_error, overwrite=overwrite)
        sys.stdout.flush()
        return (OUTPUT_RAW_SNPS_VCF, OUTPUT_RAW_INDELS_VCF)

    #Step 10: Filter SNPs
    def filterSNPs(self, inp, out="filtered_snps.vcf", re1="Step 10", 
                   showCommand=False, exit_on_error=True, overwrite=False): 
        print("Executing "+re1)
        print(colored("Filtering SNPs...", "yellow"))
        print("Following files are needed: ", self.REF_GENOME, inp, sep="\n")
        OUTPUT_FILTERED_SNPS_VCF = os.path.join(self.OUTPUT, out)
        print("Output:", colored(OUTPUT_FILTERED_SNPS_VCF, "yellow"))
        condition="'QD < {} || FS > {} || MQ < {} || MQRankSum < {} || ReadPosRankSum < {} || SOR > {}'".format(
            self.QD_SNP, self.FS_SNP, self.MQ_SNP, self.MQRankSum_SNP, self.ReadPosRankSum_SNP, self.SOR_SNP)
        command = "%r"%self.JAVA+" -XX:ParallelGCThreads="+str(self.NO_OF_INTERNAL_THREAD)+" -jar "+ \
                  "%r"%self.GATK_JAR+" -T VariantFiltration "+" -R "+"%r"%self.REF_GENOME+" -V "+ \
                  "%r"%inp+" --filterExpression "+condition+" --filterName 'basic_snp_filter' -o "+ \
                  "%r"%OUTPUT_FILTERED_SNPS_VCF
        if(showCommand):
            print("Command: \n",command)
        self.executeCommand(OUTPUT_FILTERED_SNPS_VCF, command, remark1=re1, 
                            remark2="VCF (filtered SNPs)", remark3="VCF (filtered SNPs)", 
                            exit_on_exception=exit_on_error, overwrite=overwrite)
        sys.stdout.flush()
        return OUTPUT_FILTERED_SNPS_VCF

    #Step 11: Filter Indels
    def filterIndels(self, inp, out="filtered_indels.vcf", re1="Step 11", 
                     showCommand=False, exit_on_error=True, overwrite=False): 
        print("Executing "+re1)
        print(colored("Filtering Indels...", "yellow"))
        print("Following files are needed: ", self.REF_GENOME, inp, sep="\n")
        OUTPUT_FILTERED_INDELS_VCF = os.path.join(self.OUTPUT, out)
        print("Output:", colored(OUTPUT_FILTERED_INDELS_VCF, "yellow"))
        condition="'QD < {} || FS > {} || ReadPosRankSum < {} || SOR > {}'".format(
                    self.QD_INDEL, self.FS_INDEL, self.ReadPosRankSum_INDEL, self.SOR_INDEL)
        command = "%r"%self.JAVA+" -XX:ParallelGCThreads="+str(self.NO_OF_INTERNAL_THREAD)+ \
                  " -jar "+"%r"%self.GATK_JAR+" -T VariantFiltration "+" -R "+"%r"%self.REF_GENOME+ \
                  " -V "+"%r"%inp+" --filterExpression "+condition+ \
                  " --filterName 'basic_indel_filter' -o "+"%r"%OUTPUT_FILTERED_INDELS_VCF
        if(showCommand):
            print("Command: \n",command)
        self.executeCommand(OUTPUT_FILTERED_INDELS_VCF, command, remark1=re1, 
                            remark2="VCF (filtered Indels)", remark3="VCF (filtered Indels)",  
                            exit_on_exception=exit_on_error, overwrite=overwrite)
        sys.stdout.flush()
        return OUTPUT_FILTERED_INDELS_VCF    

    #Step 12: Base Quality Score Recalibration (BQSR) 1
    def scoreRecalibration(self, exit_on_error=True, overwrite=False): 
        print("Executing Step 12")
        print(colored("Base Quality Score Recalibration (BQSR) #1...", "yellow"))
        print("Following files are needed: ", self.REF_GENOME, self.OUTPUT_REALIGNED_BAM, 
              self.OUTPUT_FILTERED_SNPS_VCF, self.OUTPUT_FILTERED_INDELS_VCF, sep="\n")
        self.OUTPUT_RECAL_DATA = os.path.join(self.OUTPUT, "recal_data.table")
        print("Output:", colored(self.OUTPUT_RECAL_DATA, "yellow"))
        command = "%r"%self.JAVA+" -XX:ParallelGCThreads="+str(self.NO_OF_INTERNAL_THREAD)+ \
            " -Xmx"+str(self.AMOUNT_OF_HEAP)+"g -jar "+"%r"%self.GATK_JAR+" -T BaseRecalibrator "+ \
            "-nct "+str(self.NO_OF_INTERNAL_THREAD)+" -R "+"%r"%self.REF_GENOME+ \
            " -I "+"%r"%self.OUTPUT_REALIGNED_BAM+" -knownSites "+"%r"%self.OUTPUT_FILTERED_SNPS_VCF+ \
            " -knownSites "+"%r"%self.OUTPUT_FILTERED_INDELS_VCF+" -o "+"%r"%self.OUTPUT_RECAL_DATA
        self.executeCommand(self.OUTPUT_RECAL_DATA, command, remark1="Step 12", 
                            remark2="table (BQSR #1)", remark3="Table (BQSR #1)",
                            exit_on_exception=exit_on_error, overwrite=overwrite)
        sys.stdout.flush()

    #Step 13: Base Quality Score Recalibration (BQSR) 2
    def scoreRecalibration2(self, exit_on_error=True, overwrite=False): 
        print("Executing Step 13")
        print(colored("Base Quality Score Recalibration (BQSR) #2...", "yellow"))
        print("Following files are needed: ", self.REF_GENOME, self.OUTPUT_REALIGNED_BAM, 
              self.OUTPUT_RECAL_DATA, self.OUTPUT_FILTERED_SNPS_VCF, self.OUTPUT_FILTERED_INDELS_VCF, sep="\n")
        self.OUTPUT_POST_RECAL_DATA = os.path.join(self.OUTPUT, "post_recal_data.table")
        print("Output:", colored(self.OUTPUT_POST_RECAL_DATA, "yellow"))
        command = "%r"%self.JAVA+" -XX:ParallelGCThreads="+str(self.NO_OF_INTERNAL_THREAD)+" -Xmx"+ \
            str(self.AMOUNT_OF_HEAP)+"g -jar "+"%r"%self.GATK_JAR+" -T BaseRecalibrator "+"-nct "+ \
            str(self.NO_OF_INTERNAL_THREAD)+" -R "+"%r"%self.REF_GENOME+" -I "+ \
            "%r"%self.OUTPUT_REALIGNED_BAM+" -knownSites "+"%r"%self.OUTPUT_FILTERED_SNPS_VCF+ \
            " -knownSites "+"%r"%self.OUTPUT_FILTERED_INDELS_VCF+" -BQSR "+ \
            "%r"%self.OUTPUT_RECAL_DATA+" -o "+"%r"%self.OUTPUT_POST_RECAL_DATA
        self.executeCommand(self.OUTPUT_POST_RECAL_DATA, command, remark1="Step 13", 
                            remark2="table (BQSR #2)", remark3="Table (BQSR #2)", 
                            exit_on_exception=exit_on_error, overwrite=overwrite)
        sys.stdout.flush()
    
    #Step 14: Analyze Covariates  
    def analyzeCovariates(self, exit_on_error=True, overwrite=False): 
        print("Executing Step 14")
        print(colored("Analyzing Covariates...", "yellow"))
        print("Following files are needed: ", self.REF_GENOME, self.OUTPUT_RECAL_DATA, 
              self.OUTPUT_POST_RECAL_DATA, sep="\n")
        self.OUTPUT_RECAL_PLOT = os.path.join(self.OUTPUT, "recalibration_plots.pdf")
        print("Output:", colored(self.OUTPUT_RECAL_PLOT, "yellow"))
        command = "%r"%self.JAVA+" -XX:ParallelGCThreads="+str(self.NO_OF_INTERNAL_THREAD)+ \
            " -Xmx"+str(self.AMOUNT_OF_HEAP)+"g -jar "+"%r"%self.GATK_JAR+" -T AnalyzeCovariates -R "+ \
            "%r"%self.REF_GENOME+" -before "+"%r"%self.OUTPUT_RECAL_DATA+" -after "+ \
            "%r"%self.OUTPUT_POST_RECAL_DATA+" -plots "+"%r"%self.OUTPUT_RECAL_PLOT
        self.executeCommand(self.OUTPUT_RECAL_PLOT, command, remark1="Step 14", 
                       remark2="PDF (recalibration plots)", remark3="PDF (recalibration plots)", 
                       exit_on_exception=exit_on_error, overwrite=overwrite)
        sys.stdout.flush()
        
    #Step 15: Apply BQSR
    def applyBQSR(self, exit_on_error=True, overwrite=False):  
        print("Executing Step 15")
        print(colored("Applying BQSR...", "yellow"))
        print("Following files are needed: ", self.REF_GENOME, self.OUTPUT_REALIGNED_BAM, 
              self.OUTPUT_RECAL_DATA, sep="\n")
        self.OUTPUT_RECAL_BAM = os.path.join(self.OUTPUT, "recal_reads.bam")
        print("Output:", colored(self.OUTPUT_RECAL_BAM, "yellow"))
        command = "%r"%self.JAVA+" -XX:ParallelGCThreads="+str(self.NO_OF_INTERNAL_THREAD)+ \
            " -Xmx"+str(self.AMOUNT_OF_HEAP)+"g -jar "+"%r"%self.GATK_JAR+" -T PrintReads -nct "+ \
            str(self.NO_OF_INTERNAL_THREAD)+" -R "+"%r"%self.REF_GENOME+" -I "+ \
            "%r"%self.OUTPUT_REALIGNED_BAM+" -BQSR "+"%r"%self.OUTPUT_RECAL_DATA+" -o "+ \
            "%r"%self.OUTPUT_RECAL_BAM
        self.executeCommand(self.OUTPUT_RECAL_BAM, command, remark1="Step 15", 
                            remark2="BAM (recalibrated)", remark3="BAM (recalibrated)", 
                            exit_on_exception=exit_on_error, overwrite=overwrite)
        sys.stdout.flush()

    #Step 20: Annotate SNPs and Predict Effects
    def annotateSNPsAndPredictEffects(self, exit_on_error=True, overwrite=False):  
        print("Executing Step 20")
        print(colored("Annotating SNPs and Predicting Effects...", "yellow"))
        print("Following files are needed: ", self.OUTPUT_FILTERED_SNPS_FINAL_VCF, sep="\n")
        self.OUTPUT_FILTERED_SNPS_FINAL_ANN_VCF = os.path.join(self.OUTPUT, "filtered_snps_final.ann.vcf")
        print("Output:", colored(self.OUTPUT_FILTERED_SNPS_FINAL_ANN_VCF, "yellow"))
        command = "java"+" -XX:ParallelGCThreads="+str(self.NO_OF_INTERNAL_THREAD)+" -Xmx"+ \
            str(self.AMOUNT_OF_HEAP)+"g -jar "+"%r"%self.SNPEFF_JAR+" -v "+"%r"%self.SNPEFF_DB+ \
            " "+"%r"%self.OUTPUT_FILTERED_SNPS_FINAL_VCF+" > "+"%r"%self.OUTPUT_FILTERED_SNPS_FINAL_ANN_VCF
        self.executeCommand(self.OUTPUT_FILTERED_SNPS_FINAL_ANN_VCF, command, remark1="Step 20",
                remark2="filtered snps, summary and genes", remark3="filtered snps, summary and genes", 
                exit_on_exception=exit_on_error, overwrite=overwrite)
        sys.stdout.flush()
        
    #Step 21: Annotate SNPs and Predict Effects
    def computeCoverageStat(self, exit_on_error=True, overwrite=False):
        print("Executing Step 21")
        print(colored("Computing Coverage Statistics...", "yellow"))
        print("Following files are needed: ", self.OUTPUT_RECAL_BAM, sep="\n")
        self.OUTPUT_GENOME_COV = os.path.join(self.OUTPUT, "genomecov.bedgraph")
        print("Output:", colored(self.OUTPUT_GENOME_COV, "yellow"))
        command = "bedtools genomecov -bga -ibam "+"%r"%self.OUTPUT_RECAL_BAM+" > "+ \
            "%r"%self.OUTPUT_GENOME_COV
        self.executeCommand(self.OUTPUT_GENOME_COV, command, remark1="Step 21", 
                remark2="Coverage Statistics (bedgraph)", remark3="Coverage Statistics (bedgraph)", 
                exit_on_exception=exit_on_error, overwrite=overwrite)
        sys.stdout.flush()
    
    #Step 22: Compile Statistics
    def computeStat(self, exit_on_error=True, overwrite=False):
        print("Executing Step 22")
        print(colored("Compiling Statistics...", "yellow"))
        print("Following files are needed: ", self.OUTPUT_ALIGN_MAT_TXT, self.OUTPUT_INSERT_MAT, 
              self.OUTPUT_RAW_SNPS_VCF, self.OUTPUT_FILTERED_SNPS_VCF, self.OUTPUT_RAW_SNPS_RECAL_VCF, 
              self.OUTPUT_FILTERED_SNPS_FINAL_VCF, self.OUTPUT_DEPTH_TXT, sep="\n")
        self.OUTPUT_REPORT_CSV = os.path.join(self.OUTPUT, "report.csv")
        print("Output:", colored(self.OUTPUT_REPORT_CSV, "yellow"))
        command = "sh "+"%r"%os.path.join(self.toolspath, "parse_metrics.sh")+" "+"%r"%(self.OUTPUT+"/")+" "+ \
            "%r"%(self.OUTPUT_STAT+"/")+" "+"%r"%self.ID
        self.executeCommand(self.OUTPUT_GENOME_COV, command, remark1="Step 21", 
                remark2="Coverage Statistics (bedgraph)", remark3="Coverage Statistics (bedgraph)", 
                exit_on_exception=exit_on_error, overwrite=overwrite)
        sys.stdout.flush()
        
    



















