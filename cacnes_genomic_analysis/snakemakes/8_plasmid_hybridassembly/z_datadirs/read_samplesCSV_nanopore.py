# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 20:59:39 2020

@author: raval
"""

import os
from pathlib import Path
import sys
import glob
import subprocess

def read_samplesCSV(spls):
    # reads in samples.csv file, format: Batch, Batich_Nanopore, Sample,Alignments,ProviderName,Patient
    hdr_check = ['Batch_Illumina', 'Batch_Nanopore', 'Sample', 'Alignments', 'ProviderName', 'Barcode', 'NumNpFiles','Patient']
    switch = "on"
    file = open(spls, 'r')
    list_path_illumina = []
    list_path_nanopore = []
    list_splID = []
    list_providerNames = []
    list_refG = []
    list_barcode = []
    list_numnpfiles = []
    list_patient = []
    for line in file:
        line = line.strip('\n').split(',')
        # Test Header. Note: Even when header wrong code continues (w/ warning), but first line not read.
        if switch == "on":
            if (line == hdr_check):
                print("Passed CSV header check")
            else:
                Warning("CSV did NOT pass header check! Code continues, but first line ignored")
            switch = "off"
            continue
        # build lists
        list_path_illumina.append(line[0])
        list_path_nanopore.append(line[1])
        list_splID.append(line[2])
        list_refG.append(line[3])
        list_providerNames.append(line[4])
        list_barcode.append(line[5])
        list_numnpfiles.append(line[6])
        list_patient.append(line[7])
    return [list_path_illumina,list_path_nanopore,list_splID,list_refG,list_providerNames,list_barcode,list_numnpfiles,list_patient] # 2019.02.08, Arolyn: removed set(list_patient) provides only unique subject IDs


def split_samplesCSV(PATH_Illumina_ls,PATH_Nanopore_ls,SAMPLE_ls,REF_Genome_ls,PROVIDER_ls,BARCODE_ls,NumNpFiles_ls,PATIENTID_ls_all):
    # Added by Arolyn, 2019.02.11
    # takes info extracted form samples.csv; saves each line of samples.csv as sample_info.csv in data/{sampleID}
    for i, sample in enumerate(SAMPLE_ls):
        # get info for this sample
        sample_info_csv_text = PATH_Illumina_ls[i] + ',' + PATH_Nanopore_ls[i] + ',' + SAMPLE_ls[i] + ',' + REF_Genome_ls[i] + ',' + PROVIDER_ls[i] + ',' + BARCODE_ls[i] + ',' + NumNpFiles_ls[i] + ',' + PATIENTID_ls_all[i]
        #print( sample )
        #print( sample_info_csv_text )
        # make data directory for this sample if it doesn't already exist
        if not(os.path.isdir('data/' + sample)):
            os.makedirs('data/' + sample, exist_ok=True)
        # check to see if this mini csv with sample info already exists
        if os.path.isfile('data/' + sample + '/sample_info.csv'):
            # if so, read file
            old_file = open('data/' + sample + '/sample_info.csv','r')
            old_info = old_file.readline()
            old_file.close()
            # check to see if the existing file is consistent with samples.csv
            if not(old_info == sample_info_csv_text):
                # if not, remove the old file and save sample info in a new file
                #print('Information file must be updated.')
                os.remove('data/' + sample + '/sample_info.csv')
                f = open('data/' + sample + '/sample_info.csv','w')
                f.write(sample_info_csv_text)
                f.close()
            #else:
            #print('Information file already updated.')
        else: # if mini csv with sample info does not already exist
            # save sample info in mini csv
            #print('Information file must be created.')
            f = open('data/' + sample + '/sample_info.csv','w')
            f.write(sample_info_csv_text)
            f.close()
            
def findfastqfile(dr,ID,filename):
    fwd=[]
    rev=[]
    potentialhits_forward=glob.glob(dr + '/' + filename +'/*1.fastq.gz')
    potentialhits_reverse=glob.glob(dr + '/' + filename +'/*2.fastq.gz')
    if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
        fwd=potentialhits_forward[0]
        rev=potentialhits_reverse[0]
    elif len(potentialhits_forward)==0 and len(potentialhits_reverse)==0:
        potentialhits_forward=glob.glob(dr + '/' + filename +'*1.fastq.gz')
        potentialhits_reverse=glob.glob(dr + '/' + filename +'*2.fastq.gz')
        if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
            fwd=potentialhits_forward[0]
            rev=potentialhits_reverse[0]
        elif len(potentialhits_forward)==0 and len(potentialhits_reverse)==0:
            potentialhits_forward=glob.glob(dr + '/' + filename +'*1.fastq')
            potentialhits_reverse=glob.glob(dr + '/' + filename +'*2.fastq')
            if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
                subprocess.run("gzip " + potentialhits_forward[0], shell=True)
                subprocess.run("gzip " + potentialhits_reverse[0], shell=True)
                fwd=potentialhits_forward[0]+'.gz'
                rev=potentialhits_reverse[0]+'.gz'
            else:
                foldername=glob.glob(dr + '/' + filename + '*')
                if foldername and os.path.isdir(foldername[0]):
                    foldername=foldername[0]
                    potentialhits_forward=glob.glob(foldername + '/*' + filename + '*1*.fastq.gz')
                    potentialhits_reverse=glob.glob(foldername + '/*' + filename + '*2*.fastq.gz')
                    if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
                        fwd=potentialhits_forward[0]
                        rev=potentialhits_reverse[0]
                    elif len(potentialhits_forward)==0 and len(potentialhits_reverse)==0:
                        print(foldername + '/*' + filename + '*2*.fastq.gz')
                        potentialhits_forward=glob.glob(foldername +  '/*' + filename + '*1*.fastq')
                        potentialhits_reverse=glob.glob(foldername + '/*' + filename + '*2*.fastq')
                        if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
                            subprocess.run("gzip " + potentialhits_forward[0], shell=True)
                            subprocess.run("gzip " + potentialhits_reverse[0], shell=True)
                            fwd=potentialhits_forward[0]+'.gz'
                            rev=potentialhits_reverse[0]+'.gz'
    if not(fwd) or not(rev):
        raise ValueError('Either no file or more than 1 file found in ' + dr + 'for ' + ID)
    ##zip fastq files if they aren't already zipped
    subprocess.run("gzip " + fwd, shell=True)
    subprocess.run("gzip " + rev, shell=True)
    return [fwd, rev]
            
def makelink(path,sample,providername):
    #When sample is run on a single lane
    #Provider name can be either a COMPLETE directory name or a file name in batch(called path in this fx)
    [fwd_file, rev_file]=findfastqfile(path,sample, providername)
    subprocess.run('ln -s -T ' + fwd_file + ' data/' + sample + '/R1.fq.gz', shell=True)
    subprocess.run('ln -s -T ' + rev_file + ' data/' + sample + '/R2.fq.gz', shell=True)

def cp_append_files(paths,sample,providername):
    #When sample is run on multiple lanes with same barcode
    fwd_list=''
    rev_list=''
    for path in paths:
        #Provider name can be either a COMPLETE directory name or a file name in batch(called path in this fx)
        [fwd_file, rev_file]=findfastqfile(path,sample, providername)
        fwd_list=fwd_list+ ' ' +fwd_file
        rev_list=rev_list+ ' ' +rev_file
        print(rev_list)
        print(fwd_list)
    subprocess.run("zcat " + fwd_list + ' | gzip > data/' +  sample + '/R1.fq.gz', shell=True)
    subprocess.run("zcat " + rev_list + ' | gzip > data/' +  sample + '/R2.fq.gz', shell=True)

def cp_append_nanopore_files(path_nanopore,sample,barcode,numnpfiles):
    fwd_list = ''
    for i in range(0,int(numnpfiles)):
        fwd_file = path_nanopore + barcode + str(i) + '.fastq'
        subprocess.run("gzip " + fwd_file, shell=True)
        fwd_list = fwd_list=fwd_list+ ' ' +fwd_file    
    subprocess.run("zcat " + fwd_list + ' | gzip > data/' +  sample + '/np.fq.gz', shell=True)

if __name__ == "__main__":
    if len(sys.argv)<1:
        print('Usage: Reads in a csv file containing information of samples to be processed.')
        print('Copies data links to data/ folder and returns for further processing.')
        print('Example: python read_move_samplesCSV.py samples.csv')
    else:
        file=sys.argv[1]
        [PATH_ls, SAMPLE_ls, REF_Genome_ls, PROVIDER_ls, PATIENTID_ls] = read_samplesCSV(file)
        for i, path in enumerate(PATH_ls):
            print(i)
            if not(os.path.isfile(SAMPLE_ls[i])):
                os.makedirs('data/' + SAMPLE_ls[i], exist_ok=True)
            paths = path.split(' ')
            if len(paths)>1:
                cp_append_files(paths, SAMPLE_ls[i], PROVIDER_ls[i])
            else:
                makelink(path, SAMPLE_ls[i], PROVIDER_ls[i])
