###
#   File name  : ReadFastq_RandomReads.py
#   Author     : Hyunjin Kim
#   Date       : Apr 12, 2023
#   Email      : firadazer@gmail.com
#   Purpose    : Read a FASTQ file, select the first 1:n reads (down-sizing), then write out as new FASTQ
#
#   * THIS CODE DOES NOT SELECT READS RANDOMLY. IT'S 1:N
#
#   Instruction
#               1. import ReadFastq_RandomReads.py
#               2. Run the function ReadFastq_RandomReads.start()
#               3. The results will be generated under the same directory as the input
###

### set parameters
input_fastq_path = ["C:/Users/hyunjin.kim2/Documents/RProjects/SimpleTasks/data/PID5248/1754485-G-1_S1_L001_R1_001.fastq",
                    "C:/Users/hyunjin.kim2/Documents/RProjects/SimpleTasks/data/PID5248/1754485-G-1_S1_L001_R2_001.fastq"]
sampling_pcnt=75
verbose=True

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
    print("ReadFastq_RandomReads.py")

    start_time = timeit.default_timer()
    for fp in input_fastq_path:
        read_and_write(fp)
    print("Execution Time: ", timeit.default_timer() - start_time)

### read fastq and write without duplicates
def read_and_write(input_path):
    if verbose:
        print("read_and_write()")

    ### read the fastq file
    records = list(SeqIO.parse(input_path, "fastq"))

    ### make a dataframe with the records
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

    ### change the df to list
    new_records = df.values.tolist()
    new_records = [SeqRecord(y[0],
                             id=y[1],
                             name=y[2],
                             description=y[3],
                             dbxrefs=y[4],
                             letter_annotations=y[5]) for y in [x for x in new_records]]

    ### down sampling using the sample_num
    sample_num = round(len(new_records) * sampling_pcnt/100)
    # random.seed(1234)
    # new_records = sample(new_records, sample_num)
    new_records = [new_records[i] for i in range(sample_num)]

    # write out the filtered sequences
    SeqIO.write(new_records, os.path.splitext(input_path)[0] + '_' + str(sampling_pcnt) + 'Pcnt_Sampled.fastq', "fastq")

start()





