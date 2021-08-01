###
#   File name  : ReadFastq_RmReps.py
#   Author     : Hyunjin Kim
#   Date       : Apr 20, 2021
#   Email      : firadazer@gmail.com
#   Purpose    : Read a FASTQ file, remove any replicates, then write out as FASTQ
#
#   Instruction
#               1. import ReadFastq_RmReps.py
#               2. Run the function ReadFastq_RmReps.start()
#               3. The results will be generated under the same directory as the input
###

### set parameters
verbose=True
input_fastq_path = ["Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/Allison_HLA/data/2-1077609/1834022_JCC223_SJ-FLC-B2_GEX_S2_L001_R2_001_HLA.fastq",
                    "Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/Allison_HLA/data/2-1077609/1834022_JCC223_SJ-FLC-B2_GEX_S2_L002_R2_001_HLA.fastq",
                    "Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/Allison_HLA/data/2-1131806/1894886_JCC209_AK_FLC_AscitesTIL_Rxn1_Gex_S1_L001_R2_001_HLA.fastq",
                    "Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/Allison_HLA/data/2-1131806/1894886_JCC209_AK_FLC_AscitesTIL_Rxn1_Gex_S1_L002_R2_001_HLA.fastq",
                    "Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/Allison_HLA/data/2-1244730/1977314_JCC223_AK-1_P9R_GEX_S11_L002_R2_001_HLA.fastq",
                    "Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/Allison_HLA/data/2-1244731/1977315_JCC223_AK-2_B3unstim_GEX_S12_L003_R2_001_HLA.fastq"]
input_fastq_path = ["Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/Allison_HLA/hisat-genotype/AK_FLC_AscitesTIL_Rxn1_Gex_S1_L001_R2_001.fastq",
                    "Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/Allison_HLA/hisat-genotype/AK_FLC_AscitesTIL_Rxn1_Gex_S1_L002_R2_001.fastq",
                    "Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/Allison_HLA/hisat-genotype/SJ-FLC-B2_GEX_S2_L001_R2_001.fastq",
                    "Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/Allison_HLA/hisat-genotype/SJ-FLC-B2_GEX_S2_L002_R2_001.fastq",
                    "Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/Allison_HLA/hisat-genotype/AK-2_B3unstim_GEX_S12_L003_R2_001.fastq",
                    "Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/Allison_HLA/hisat-genotype/AK-1_P9R_GEX_S11_L002_R2_001.fastq"]
input_fastq_path = ["Z:/ResearchHome/ProjectSpace/thomagrp/JCC282_Hematopoiesis/common/HK/2156502_JCC223_AK-FLC02_GEX_S4_L001_R2_001_HLA.fastq"]
input_fastq_path = ["Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/1912291_JCC212_SJCAR19-11_Wk-1_PB_Gex_S5_L001_R2_001_HLA.fastq",
                    "Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/1967963_JCC212_SJCAR19-11_GMP_GEX_S5_L001_R2_001_HLA.fastq",
                    "Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/1929015_JCC212_SJCAR19-11_Wk1_PB_Gex_S6_L001_R2_001_HLA.fastq",
                    "Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/1929016_JCC212_SJCAR19-11_Wk2_PB_Gex_S7_L001_R2_001_HLA.fastq",
                    "Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/1932046_JCC212_SJCAR19-11_Wk3_PB_Gex_S1_L001_R2_001_HLA.fastq"]
sample_num=300000

### import modules
import timeit
from Bio import SeqIO
import pandas as pd
from Bio.SeqRecord import SeqRecord
import os
from random import sample
import random

### a function starting this script
def start():
    print("ReadFastq_RmReps.py")

    start_time = timeit.default_timer()
    for fp in input_fastq_path:
        read_and_write(fp)
    print("Execution Time: ", timeit.default_timer() - start_time)

### read fastq and write without duplicates
def read_and_write(input_path):
    if verbose:
        print("read_and_write()")

    # read the fastq file
    records = list(SeqIO.parse(input_path, "fastq"))

    # make a dataframe with the records
    seqs = [y for y in [x.seq for x in records]]
    ids = [y for y in [x.id for x in records]]
    names = [y for y in [x.name for x in records]]
    dscrs = [y for y in [x.description for x in records]]
    dbxrefs = [y for y in [x.dbxrefs for x in records]]
    lannos = [y for y in [x.letter_annotations for x in records]]
    df = pd.DataFrame(
        {'Seq': seqs,
         'id': ids,
         'name': names,
         'description': dscrs,
         'dbxrefs': dbxrefs,
         'letter_annotations': lannos
        }
    )

    # remove duplicated sequences
    df = df.drop_duplicates(subset=['Seq'])

    # change the df to list
    new_records = df.values.tolist()
    new_records = [SeqRecord(y[0],
                             id=y[1],
                             name=y[2],
                             description=y[3],
                             dbxrefs=y[4],
                             letter_annotations=y[5]) for y in [x for x in new_records]]

    ### down sampling using the sample_num
    if len(new_records) > sample_num:
        random.seed(1234)
        new_records = sample(new_records, sample_num)

    # write out the filtered sequences
    SeqIO.write(new_records, os.path.splitext(input_path)[0] + '_RepRemoved_' + str(int(sample_num/1000)) + 'k_Sampled.fastq', "fastq")

start()
