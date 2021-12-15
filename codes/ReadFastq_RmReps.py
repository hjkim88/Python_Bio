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
input_fastq_path = ["Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/JCC212_SJCAR19-02_GMP-redo_short.thomagrp_192790_10x-1.2-1239786.GEX_S13_L003_R2_001.fastq",
                    "Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/JCC212_SJCAR19-11_GMP-redo_short.thomagrp_192790_10x-1.2-1239779.GEX_S6_L001_R2_001.fastq"]
input_fastq_path = ["Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/JCC212_SJCAR19-02_PB_Wk2_short.thomagrp_154684_10x-1.1626554.Wk2_S1_L001_R2_001.fastq"]
input_fastq_path = ["Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/JCC212_SJCAR19-10_Wk1_PB_short.thomagrp_182465_10x-1.2-1145865.PB_Gex_S3_L001_R2_001.fastq",
                    "Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/JCC212_SJCAR19-03_PB_Wk1_short.thomagrp_156295_10x-1.2-826754.Wk1_GEX_2_S6_L001_R2_001.fastq",
                    "Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/JCC212_SJCAR19-05_9mo_PB_short.thomagrp_183630_10x-1.2-1161888.PB_Gex_S2_L001_R2_001.fastq",
                    "Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/JCC212_SJCAR19_04_PreTransB_short.thomagrp_164732_10x-2.2-945704.GEX_S10_L001_R2_001.fastq"]
input_fastq_path = ["Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/JCC212_SJCAR19-06_Wk-1_short.thomagrp_167553_premade10x-1.2-985110.GEX_S5_L001_R2_001.fastq"]
input_fastq_path = ["Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-00_GMP18015_short.thomagrp_156667_10x-1.1646454.GEX_S2_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-01_GMP18018_short.thomagrp_156667_10x-1.1646455.GEX_S3_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-02_BM_Wk4_short.thomagrp_154684_10x-1.1626557.Wk4_S4_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-02_GMP-redo_short.thomagrp_192790_10x-1.2-1239786.GEX_S13_L003_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-02_GMP18024_short.thomagrp_156295_10x-1.2-826757.GEX_1_S9_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-02_GMP18024_short.thomagrp_156295_10x-1.2-826758.GEX_2_S10_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-02_GMP18024_short.thomagrp_156295_10x-1.2-826759.GEX_3_S11_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-02_GMP18024_short.thomagrp_156295_10x-1.2-826760.GEX_4_S12_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-02_PB_Wk1_short.thomagrp_152456_10x-1.2-826758.1_S13_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-02_PB_Wk1_short.thomagrp_152456_10x-1.2-826759.2_S14_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-02_PB_Wk1_short.thomagrp_152456_10x-1.2-826760.3_S15_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-02_PB_Wk1_short.thomagrp_152456_10x-1.2-826761.4_S16_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-02_PB_Wk1_short.thomagrp_152456_10x-2.2-826749.1_S1_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-02_PB_Wk1_short.thomagrp_152456_10x-2.2-826750.2_S2_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-02_PB_Wk1_short.thomagrp_152456_10x-2.2-826751.3_S3_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-02_PB_Wk1_short.thomagrp_152456_10x-2.2-826752.4_S4_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-02_PB_Wk2_short.thomagrp_154684_10x-1.1626554.Wk2_S1_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-02_PB_Wk3_short.thomagrp_154684_10x-1.1626555.Wk3_S2_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-02_PB_Wk4_short.thomagrp_154684_10x-1.1626556.Wk4_S3_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-03_BM_Wk4_short.thomagrp_157960_10x-1.1660485.Wk4_GEX_S2_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-03_GMP-redo_short.thomagrp_192790_10x-1.2-1239787.GEX_S14_L003_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-03_GMP19002_short.thomagrp_157960_10x-1.1660486.GEX_S3_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-03_PB_Wk-1_short.thomagrp_155832_10x-1.2-852140.Wk-1_S1_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-03_PB_Wk-1_short.thomagrp_155832_10x-1.2-852141.Wk-1_S2_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-03_PB_Wk-1_short.thomagrp_155832_10x-1.2-852142.Wk-1_S3_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-03_PB_Wk-1_short.thomagrp_155832_10x-1.2-852143.Wk-1_S4_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-03_PB_Wk1_short.thomagrp_156295_10x-1.2-826753.Wk1_GEX_1_S5_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-03_PB_Wk1_short.thomagrp_156295_10x-1.2-826754.Wk1_GEX_2_S6_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-03_PB_Wk1_short.thomagrp_156295_10x-1.2-826755.Wk1_GEX_3_S7_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-03_PB_Wk1_short.thomagrp_156295_10x-1.2-826756.Wk1_GEX_4_S8_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-03_PB_Wk2_short.thomagrp_156667_10x-1.1646453.Wk2_GEX_S1_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-03_PB_Wk3_short.thomagrp_157069_10x-1.1651508.Wk3_GEX_S1_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-03_PB_Wk4_short.thomagrp_157960_10x-1.1660484.Wk4_GEX_S1_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-03_PB_Wk6_short.thomagrp_159035_10x-1.1673745.Wk6_GEX_S1_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-04_BM_Wk4_short.thomagrp_160509_10x-1.2-900344.Wk4_GEX_S3_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-04_GMP-redo_short.thomagrp_192790_10x-1.2-1239788.GEX_S15_L003_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-04_PB_Wk-1_short.thomagrp_157960_10x-1.1660487.GEX_S4_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-04_PB_Wk1_short.thomagrp_159035_10x-1.1673746.Wk1_GEX_S2_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-04_PB_Wk1b_short.thomagrp_159035_10x-1.1673747.Wk1b_GEX_S3_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-04_PB_Wk4_short.thomagrp_160509_10x-1.2-900343.Wk4_GEX_S2_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-04_Wk3_short.thomagrp_160509_10x-1.2-900342.GEX_S1_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-05_6mo_PB_short.thomagrp_176316_10x-1.2-1088444.PB_Gex_S11_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-05_9mo_PB_short.thomagrp_183630_10x-1.2-1161888.PB_Gex_S2_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-05_BM_3mo_short.thomagrp_168186_premade10x-1.2-1003613.3mo_GEX_S5_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-05_GMP-redo_short.thomagrp_192790_10x-1.2-1239789.GEX_S16_L003_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-05_GMP19016_short.thomagrp_167302_premade10x-1.2-985107.GEX_S2_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-05_PB_3mo_short.thomagrp_168186_premade10x-1.2-1003612.3mo_GEX_S4_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-05_Wk8_short.thomagrp_167302_premade10x-1.2-985106.GEX_S1_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-06_3mo_BM_short.thomagrp_176316_10x-1.2-1088449.BM_Gex_S16_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-06_3mo_PB_short.thomagrp_176316_10x-1.2-1088448.PB_Gex_S15_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-06_DisEval_BM_short.thomagrp_182465_10x-1.2-1145864.BM_Gex_S2_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-06_DisEval_PB_short.thomagrp_182465_10x-1.2-1145863.PB_Gex_S1_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-06_GMP-redo_short.thomagrp_192790_10x-1.2-1239780.GEX_S7_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-06_GMP19028_short.thomagrp_168186_premade10x-1.2-1003611.GEX_S3_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-06_Wk-1_short.thomagrp_167553_premade10x-1.2-985110.GEX_S5_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-06_Wk1_PB_short.thomagrp_169630_10x-1.2-1022634.PB_GEX_S1_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-06_Wk2_PB_short.thomagrp_169630_10x-1.2-1022635.PB_GEX_S2_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-06_Wk3_PB_short.thomagrp_169630_10x-1.2-1022639.PB_GEX_S6_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-06_Wk4_BM_short.thomagrp_176316_10x-1.2-1088435.BM_Gex_S2_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-06_Wk4_PB_short.thomagrp_176316_10x-1.2-1088434.PB_Gex_S1_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-06_Wk8_PB_short.thomagrp_176316_10x-1.2-1088437.PB_Gex_S4_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-07_GMP-redo_short.thomagrp_192790_10x-1.2-1239781.GEX_S8_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-07_GMP19047_short.thomagrp_176316_10x-1.2-1088441.Gex_S8_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-07_Wk-1_PB_short.thomagrp_176316_10x-1.2-1088436.PB_Gex_S3_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-07_Wk1_PB_short.thomagrp_176316_10x-1.2-1088438.PB_Gex_S5_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-07_Wk2_PB_short.thomagrp_176316_10x-1.2-1088439.PB_Gex_S6_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-07_Wk3_PB_short.thomagrp_176316_10x-1.2-1088440.PB_Gex_S7_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-07_Wk4_BM_short.thomagrp_176316_10x-1.2-1088443.BM_Gex_S10_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-07_Wk4_PB_short.thomagrp_176316_10x-1.2-1088442.PB_Gex_S9_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-07_Wk8_BM_short.thomagrp_176316_10x-1.2-1088447.BM_Gex_S14_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-07_Wk8_PB_short.thomagrp_176316_10x-1.2-1088446.PB_Gex_S13_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-08_GMP-redo_short.thomagrp_192790_10x-1.2-1239782.GEX_S9_L003_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-08_GMP19064_short.thomagrp_186310_10x-1.2-1188377.Gex_S2_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-08_Wk-1_PB_short.thomagrp_176316_10x-1.2-1088445.PB_Gex_S12_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-08_Wk1_PB_short.thomagrp_179433_10x-1.2-1113620.PB_Gex_S1_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-08_Wk2_PB_short.thomagrp_179433_10x-1.2-1113621.PB_Gex_S2_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-08_Wk3_PB_short.thomagrp_179433_10x-1.2-1113624.PB_Gex_S5_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-08_Wk4_BM_short.thomagrp_179433_10x-1.2-1113626.BM_Gex_S7_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-08_Wk4_PB_short.thomagrp_179433_10x-1.2-1113625.PB_Gex_S6_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-08_Wk8_PB_short.thomagrp_179433_10x-1.2-1113628.PB_Gex_S9_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-09_GMP-redo_short.thomagrp_192790_10x-1.2-1239777.GEX_S4_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-09_GMP_short.thomagrp_192790_10x-1.2-1239783.PB_GEX_S10_L003_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-09_Wk-1Run1_PB_short.thomagrp_179433_10x-1.2-1113622.Gex_S3_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-09_Wk-1Run2_PB_short.thomagrp_179433_10x-1.2-1113623.Gex_S4_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-09_Wk1_PB_short.thomagrp_179433_10x-1.2-1113627.PB_Gex_S8_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-09_Wk2_PB_short.thomagrp_179433_10x-1.2-1113629.PB_Gex_S10_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-09_Wk3_PB_short.thomagrp_181333_10x-1.2-1131743.PB_Gex_S1_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-09_Wk4_PB_short.thomagrp_181333_10x-1.2-1131744.PB_Gex_S2_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-10_GMP-redo_short.thomagrp_192790_10x-1.2-1239778.GEX_S5_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-10_GMP_short.thomagrp_192790_10x-1.2-1239784.GEX_S11_L003_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-10_Wk-1_PB_short.thomagrp_181333_10x-1.2-1131745.PB_Gex_S3_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-10_Wk1_PB_short.thomagrp_182465_10x-1.2-1145865.PB_Gex_S3_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-10_Wk2_PB_short.thomagrp_182465_10x-1.2-1145866.PB_Gex_S4_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-10_Wk3_PB_short.thomagrp_183630_10x-1.2-1161887.PB_Gex_S1_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-10_Wk4_BM_short.thomagrp_183630_10x-1.2-1161890.BM_Gex_S4_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-10_Wk4_PB_short.thomagrp_183630_10x-1.2-1161889.PB_Gex_S3_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-10_Wk4_PB_short.thomagrp_200167_10x-1.2-1301843.Gex_S1_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-10_Wk8_PB_short.thomagrp_185924_10x-1.2-1185163.PB_Gex_S5_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-11_GMP-redo_short.thomagrp_192790_10x-1.2-1239779.GEX_S6_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-11_GMP_short.thomagrp_192790_10x-1.2-1239785.GEX_S12_L003_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-11_Wk-1_PB_short.thomagrp_183630_10x-1.2-1161891.PB_Gex_S5_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-11_Wk1_PB_short.thomagrp_185924_10x-1.2-1185164.PB_Gex_S6_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-11_Wk2_PB_short.thomagrp_185924_10x-1.2-1185165.PB_Gex_S7_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-11_Wk3_PB_short.thomagrp_186310_10x-1.2-1188376.PB_Gex_S1_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-11_Wk4_PB_short.thomagrp_186791_10x-1.2-1197360.PB_Gex_S5_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-12_GMP_short.thomagrp_194001_10x-1.2-1244727.GEX_S8_L002_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-12_GMP_short.thomagrp_200168_10x-1.2-1301845.S3_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-12_Wk-1_PB_short.thomagrp_187419_10x-1.2-1206137.PB_Gex_S1_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-13_GMP_short.thomagrp_194001_10x-1.2-1244728.GEX_S9_L002_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-13_GMP_short.thomagrp_200168_10x-1.2-1301846.S4_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-13_Wk-1_PB_short.thomagrp_192790_10x-1.2-1239776.GEX_S3_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-13_Wk3_PB_short.thomagrp_197185_10x-1.2-1272743.PB_GEX_S11_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-13_Wk4_BM_short.thomagrp_197185_10x-1.2-1272746.BM_GEX_S14_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-13_Wk4_PB_short.thomagrp_197185_10x-1.2-1272745.PB_GEX_S13_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-14_Wk1_PB_short.thomagrp_197185_10x-1.2-1272744.PB_GEX_S12_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-14_Wk2_PB_short.thomagrp_197185_10x-1.2-1272747.PB_GEX_S15_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-14_Wk3_PB_short.thomagrp_197185_10x-1.2-1272749.PB_GEX_S17_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-14_Wk4_BM_short.thomagrp_197185_10x-1.2-1272803.BM_GEX_S21_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-14_Wk4_PB_short.thomagrp_197185_10x-1.2-1272802.PB_GEX_S20_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-15_Wk-1_PB_short.thomagrp_200168_10x-1.2-1301844.GEX_S2_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-15_Wk1_PB_short.thomagrp_197185_10x-1.2-1272748.PB_GEX_S16_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-15_Wk2_PB_short.thomagrp_197185_10x-1.2-1272750.PB_GEX_S18_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-15_Wk3_PB_short.thomagrp_197185_10x-1.2-1272801.PB_GEX_S19_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-15_Wk4_BM_short.thomagrp_197185_10x-1.2-1272805.BM_GEX_S23_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-15_Wk4_PB_short.thomagrp_197185_10x-1.2-1272804.PB_GEX_S22_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-6_PreTrans_short.thomagrp_169630_10x-1.2-1022638.GEX_S5_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_00_PreTransB_short.thomagrp_164732_10x-1.2-945705.GEX_S2_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_00_PreTransB_short.thomagrp_164732_10x-4.2-945705.GEX_S3_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_03_PreTransB_short.thomagrp_164732_10x-2.2-945703.GEX_S9_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_04_CD45neg_Wk2_short.thomagrp_160895_10x-1.2-906024.CD45neg_S4_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_04_GMP19003_short.thomagrp_159627_10x-1.1678894.GEX_S2_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_04_PreTransB_short.thomagrp_164732_10x-2.2-945704.GEX_S10_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_04_Wk2_short.thomagrp_159627_10x-1.1678893.GEX_S1_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_04_Wk8_short.thomagrp_163075_10x-1.2-930448.GEX_S1_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_05_BM_Wk4_short.thomagrp_164732_10x-2.2-945702.Wk4_GEX_S8_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_05_PB_Wk0_short.thomagrp_161275_10x-1.2-914748.GEX_S1_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_05_PB_Wk4_short.thomagrp_164732_10x-2.2-945701.Wk4_GEX_S7_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_05_PreTransA_short.thomagrp_164732_10x-1.2-945710.GEX_S7_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_05_PreTransB_short.thomagrp_164732_10x-1.2-945706.GEX_S3_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_05_Wk1_short.thomagrp_164732_10x-2.2-945649.GEX_S5_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_05_Wk2_short.thomagrp_163075_10x-1.2-930449.GEX_S2_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_05_Wk2_short.thomagrp_164732_10x-2.2-945650.GEX_S6_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_05_Wk3_short.thomagrp_163075_10x-1.2-930450.GEX_S3_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_Donor30_PreTransB_short.thomagrp_164732_10x-1.2-945707.GEX_S4_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_Donor32_PreTrans_short.thomagrp_162288_10x-1.2-923746.GEX_S2_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_Donor32_PreTransB_short.thomagrp_164732_10x-1.2-945708.GEX_S5_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_Donor32_PreTransB_short.thomagrp_164732_10x-4.2-945708.GEX_S4_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_Donor33_PreTrans_short.thomagrp_162288_10x-1.2-923747.GEX_S3_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_Donor33_PreTransB_short.thomagrp_164732_10x-1.2-945709.GEX_S6_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_GMPdonor30_short.thomagrp_155832_10x-1.2-852144.L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_GMPdonor30_short.thomagrp_155832_10x-1.2-852145.L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_GMPdonor30_short.thomagrp_155832_10x-1.2-852146.L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_GMPdonor30_short.thomagrp_155832_10x-1.2-852147.L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_GMPdonor32_short.thomagrp_155832_10x-1.2-852148.L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_GMPdonor32_short.thomagrp_155832_10x-1.2-852149.L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_GMPdonor32_short.thomagrp_155832_10x-1.2-852150.L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_GMPdonor32_short.thomagrp_155832_10x-1.2-852151.L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_GMPdonor32_short.thomagrp_155832_10x-2.2-852142.L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_GMPdonor33_short.thomagrp_155832_10x-1.2-852152.L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_GMPdonor33_short.thomagrp_155832_10x-1.2-852153.L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_GMPdonor33_short.thomagrp_155832_10x-1.2-852154.L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_GMPdonor33_short.thomagrp_155832_10x-1.2-852155.L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_GMPdonor33_short.thomagrp_155832_10x-2.2-852143.L001_R2_001.fastq"]
input_fastq_path = ["Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-02_PreTrans_short.thomagrp_161275_10x-1.2-914749.GEX_S2_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-05_GMP-redo_short.thomagrp_192790_10x-1.2-1239789.GEX_S16_L003_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-05relapse-2ndInfusionWk1_PB_short.thomagrp_194001_10x-1.2-1244723.GEX_S4_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-05relapse-2ndInfusionWk2_PB_short.thomagrp_194001_10x-1.2-1244724.GEX_S5_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-07_Wk8_BM_short.thomagrp_176316_10x-1.2-1088447.BM_Gex_S14_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-11_Wk4_BM_short.thomagrp_186791_10x-1.2-1197361.BM_Gex_S6_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-12_Wk1_PB_short.thomagrp_194001_10x-1.2-1244729.PB_GEX_S10_L002_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-13_Wk2_PB_short.thomagrp_194001_10x-1.2-1244725.PB_GEX_S6_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-13_Wk3_PB_short.thomagrp_197185_10x-1.2-1272743.PB_GEX_S11_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-15_Wk-1_PB_short.thomagrp_194001_10x-1.2-1244726.PB_GEX_S7_L002_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_05_PB_Wk-1_short.thomagrp_160895_10x-1.2-906023.Wk-1_GEX_S3_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19_05_Wk1_short.thomagrp_162288_10x-1.2-923745.Wk1_GEX_S1_L001_R2_001.fastq"]
input_fastq_path = ["Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/HLA_Run/JCC212_SJCAR19-02_PB_Wk1_short.thomagrp_152456_10x-1.2-826758.1_S13_L001_R2_001.fastq"]
input_fastq_path = ["Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/JCC212_SJCAR19-07_Wk1_PB_short.thomagrp_176316_10x-1.2-1088438.PB_Gex_S5_L001_R2_001.fastq"]
input_fastq_path = ["Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/CATCHAML/2227092_JCC319_CATCHAML-05_GMP_GEX_S2_L001_R2_001.fastq",
"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/CATCHAML/2194808_JCC319_CATCHAML-05_Wk-1_GEX_S1_L001_R2_001.fastq"]
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
