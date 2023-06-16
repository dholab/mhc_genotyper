# -*- coding: utf-8 -*-
import numpy as np
import logging
import tempfile
import os
import shutil
import subprocess
import time
import xlsxwriter
import pandas as pd
from pathlib import Path
from datetime import datetime
from argparse import ArgumentParser, ArgumentTypeError


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ArgumentTypeError('Boolean value expected.')

# project column name on sample sheet


# i have the bbmap directory symbolic path  saved to /usr/local/bin


def parse_process_arrays_args(parser: ArgumentParser):
    """Parses the python script arguments from bash and makes sure files/inputs are valid"""
    parser.add_argument('-r', '--run_id',
                        type=str,
                        help='Required: id for your run',
                        required=True)
    parser.add_argument('-i', '--input_dir',
                        type=str,
                        help='Required: directory where your R1, R2 fastq (paired-end read) fastq files exist. This does not search Subdirectories',
                        required=True)
    parser.add_argument('-g', '--genotyper_root_dir',
                        type=str,
                        help='Required: Directory full path of this repository',
                        required=True)
    parser.add_argument('--output_dir',
                        type=str,
                        help='Optional: directory where your output_dir',
                        default=None,
                        required=False)
    parser.add_argument('--experiment',
                        type=str,
                        help='Optional: Experiment Number if different from run_id, or will set as same',
                        default=None,
                        required=False)
    parser.add_argument('--BBDUK_PATH',
                        type=str,
                        help='Optional: path where bbduk.sh is if not symbolically linked to bbduk.sh',
                        default='bbduk.sh',
                        required=False)
    parser.add_argument('--STATS_PATH',
                        type=str,
                        help='Optional: path where bbduk.sh is if not symbolically linked to default stats.sh',
                        default='stats.sh',
                        required=False)

    parser.add_argument('--BBMAP_PATH',
                        type=str,
                        help='Optional: path where bbduk.sh is if not symbolically linked to bbmap.sh',
                        default='bbmap.sh',
                        required=False)
    parser.add_argument('--BBMERGE_PATH',
                        type=str,
                        help='Optional: path where bbduk.sh is if not symbolically linked to bbmerge.sh',
                        default='bbmerge.sh',
                        required=False)
    parser.add_argument('--sample_sheet_path',
                        type=str,
                        help='Use full path if there is multiple .csv in your input_dir or if you have you samplesheet not in your input_dir',
                        default=None,
                        required=False)

    parser.add_argument('--VSEARCH_PATH',
                        type=str,
                        help='Optional: path used to trigger vsearch',
                        default='vsearch',
                        required=False)

    parser.add_argument('--project_column',
                        type=str,
                        help='Optional: if on your Sample_Project the column is not Sample_Project',
                        default='Sample_Project',
                        required=False)
    parser.add_argument('--species_column',
                        type=str,
                        help='Optional: if on your sample sheet the Species column is not Species',
                        default='Species',
                        required=False)



def get_process_arrays_args():
    """	Inputs arguments from bash
    Gets the arguments, checks requirements, returns a dictionary of arguments
    Return: args - Arguments as a dictionary
    """
    parser = ArgumentParser()
    parse_process_arrays_args(parser)
    args = parser.parse_args()
    return args


args = get_process_arrays_args()

EXPERIMENT = args.experiment
RUN_ID = args.run_id
if EXPERIMENT is None:
    EXPERIMENT = RUN_ID
input_dir = args.input_dir
output_dir = args.output_dir
if output_dir is None:
    output_dir = os.path.join(input_dir, 'out')
genotyper_root_dir=args.genotyper_root_dir
BBDUK_PATH = args.BBDUK_PATH
STATS_PATH = args.STATS_PATH
BBMAP_PATH = args.BBMAP_PATH
BBMERGE_PATH = args.BBMERGE_PATH
sample_path = args.sample_sheet_path


vsearch_path = args.VSEARCH_PATH
USEARCH_PATH = vsearch_path
species_column = args.species_column
project_column = args.project_column
project_list = []  # optional must be perfectly spelled and is case sensitve as on your sample sheet


if os.path.exists(input_dir):
    print("Make sure fastq.gz and samplesheet (.csv) have been successfully copied: ")
    print(input_dir)
else:
    print("creat a input folder and copy  fastq.gz and samplesheet (.csv) to it in the directory:")
    print(input_dir)
    raise ('INPUT FOLDER DOES NOT EXIST!')
file_list = os.listdir(input_dir)
fastq_list = [x for x in file_list if (x.endswith('.fastq.gz') or x.endswith('.fq.gz')) and not x.startswith('._')]
if sample_path is None:
    sample_list = [x for x in file_list if x.endswith('.csv') and not x.startswith('._')]
    if len(sample_list) > 1:
        raise ('Too many .csv files in input older cannot determine samplesheet!')
    elif len(sample_list) == 1:
        sample_path = os.path.join(input_dir, sample_list[0])
        print('Sample sheet path: {0}'.format(sample_path))
    else:
        print('NO SAMPLE SHEET FOUND')
        raise ('INPUT FOLDER DOES NOT EXIST!')
else:
    if not os.path.exists(sample_path):
        raise ('Sample path declared does not exist, this needs to be a full path!')

"""# PIPELINE BEGINS: DO EDIT BELOW THIS SECTION!

## Output files:

### Dependencies

+ Jupyter Notebook/Jupyter Lab
+ Python 3 (tested on anaconda distribution of Python 3.6.4)
+ Access to dholk.primate.wisc.edu
+ pigz (in PATH)
+ bbmap (in PATH)
+ bbmerge (in PATH)
+ bbduk (in PATH)
+  DO NOT USE: USEARCH v10 (discontinued)
+ vsearch 
    + ! wget https://github.com/torognes/vsearch/releases/download/v2.21.1/vsearch-2.21.1-linux-x86_64.tar.gz
+ Pandas (tested from anaconda distribution)

## Configure file paths"""
ref_dict = {
    'REF': {
        'MCM': os.path.join(genotyper_root_dir, 'ref/MCM_MHC-all_mRNA-MiSeq_singles-RENAME_20Jun16.fasta'),
        'MAMU_04_21': os.path.join(genotyper_root_dir, 'ref/25533_ipd_miSeq_deduplicated_6Apr21.fasta'),
        'MANE': os.path.join(genotyper_root_dir, 'ref/Mane_MiSeq-IPD_17.06.01_2.2.0.0_plus_SP_RW.fasta'),
        'MAMU': os.path.join(genotyper_root_dir, 'ref/26128_ipd-mhc-mamu-2021-07-09.miseq.RWv4.fasta')
    },
    'PCR_PRIMERS': os.path.join(genotyper_root_dir, 'ref/SBT195_MHCII_primers_2Sep13.fasta')
}

"""# Results and cohort sorting intermediate input files structure
- The input files will be automatically moved by cohort. This is because the pivot tables are contructed by cohort and this will be one less step.
    - they will go to <input_dir>/input_fq_by_cohort/*miseqnumber*/*cohort_name*

# split up main dir into different cohorts
- open the sample sheet
- look split by same cohorts
- mkdir for each cohort
- mv files to each by cohort
- use cyber duck's webdav to download the files from illumina
"""

# create nested dictionaries of dictionaries to support haplotyping against multiple haplotype definitions

### Indian rhesus ###
indian_rhesus = {'PREFIX': 'Mamu'}

# Indian rhesus MHC updated by Roger 24 August 2021

indian_rhesus['MHC_A_HAPLOTYPES'] = {
    'A001.01': ['A1_001'],
    'A002.01': ['A1_002_01'],
    'A003.01': ['A1_003'],
    'A004.01': ['A1_004'],
    'A006.01': ['A1_006'],
    'A007.01': ['A1_007'],
    'A008.01': ['A1_008'],
    'A011.01': ['A1_011'],
    'A012.01': ['A1_012'],
    'A016.01': ['A1_016'],
    'A018.01': ['A1_018'],
    'A018.02': ['A1_018', 'A2_01'],
    'A019.01': ['A1_019'],
    'A019.02': ['A1_019_11', 'A1_003'],
    'A022.01': ['A1_022'],
    'A023.01': ['A1_023'],
    'A025.01': ['A1_025'],
    'A026.01': ['A1_026'],
    'A028.01': ['A1_028g'],
    'A055.01': ['A1_055'],
    'A074.01': ['A1_074'],
    'A110-A111.01': ['A1_110_A1_111'],
    'A224.01': ['A2_24', 'A1_003'],
}

indian_rhesus['MHC_B_HAPLOTYPES'] = {
    'B001.01': ['B_001', 'B_007', 'B_030'],
    'B001.03': ['B_001_02', 'B_094', 'B_095'],
    'B002.01': ['B_002'],
    'B008.01': ['B_008', 'B_006'],
    'B012.01': ['B_012', 'B_030', 'B_082'],
    'B012.02': ['B_012', 'B_022', 'B_030'],
    'B012.03': ['B_012', 'B_022', 'B_030', 'B_031g'],
    'B015.01': ['B_015g2', 'B_005g'],
    'B015.02': ['B_015g2', 'B_068g1'],
    'B015.03': ['B_015g2', 'B_031g', 'B_068g1'],
    'B017.01': ['B_017', 'B_029'],
    'B017.02': ['B_017', 'B_065', 'B_083'],
    'B017.04': ['B_017', 'B_065', 'B_068', 'B_083'],
    'B024.01': ['B_024', 'B_019'],
    'B028.01': ['B_028', 'B_021'],
    'B043.01': ['B_043', 'B_030'],
    'B043.02': ['B_043', 'B_030', 'B_031_03', 'B_073'],
    'B043.03': ['B_043', 'B_030', 'B_073'],
    'B045.01': ['B_045', 'B_037'],
    'B047.01': ['B_047'],
    'B048.01': ['B_048', 'B_041'],
    'B055.01': ['B_055', 'B_052', 'B_058'],
    'B056.01': ['B_056', 'B_067'],
    'B056.02': ['B_056', 'B_066', 'B_068'],
    'B069.01': ['B_069', 'B_065'],
    'B069.02': ['B_069', 'B_068', 'B_075'],
    'B071.01': ['B_047_B_071', 'B_006'],
    'B080.01': ['B_080', 'B_081'],
    'B091.01': ['B_091', 'B_068'],
    'B093.01': ['B_093'],
    'B106.01': ['B_106', 'B_033'],
}

indian_rhesus['MHC_DRB_HAPLOTYPES'] = {
    'DR01.01': ['DRB1_04_06_01', 'DRB5_03_01'],
    'DR01.03': ['DRB1_04_11', 'DRB5_03_09'],
    'DR01.04': ['DRB1_04_06_01', 'DRB5_03_09'],
    'DR02.01': ['DRB3_04_03', 'DRB_W003_05'],
    'DR03.01': ['DRB1_03_03_01', 'DRB1_10_07'],
    'DR03.02': ['DRB1_03_12', 'DRB1_10_07'],
    'DR03.03': ['DRB1_03_17', 'DRB1_10_08'],
    'DR03.04': ['DRB1_03_18', 'DRB1_10_03'],
    'DR03.05': ['DRB1_03_06', 'DRB1_10_07'],
    'DR03.06': ['DRB1_03_06', 'DRB1_10_03'],
    'DR03.07': ['DRB1_03_19', 'DRB1_10_03'],
    'DR03.08': ['DRB1_03_03_01', 'DRB1_10_03'],
    'DR03.09': ['DRB1_03_20', 'DRB1_10_02_02'],
    'DR04.01': ['DRB1_03_09', 'DRB_W002_01'],
    'DR04.02': ['DRB1_03_18', 'DRB_W002_01'],
    'DR04.03': ['DRB1_03_09', 'DRB_W002_03'],
    'DR05.01': ['DRB1_04_03', 'DRB_W005_01'],
    'DR05.02': ['DRB1_04_03', 'DRB_W005_02'],
    'DR06.01': ['DRB_W003_03', 'DRB_W004_01'],
    'DR08.01': ['DRB_W028_01', 'DRB3_04_09', 'DRB5_03_07'],
    'DR09.01': ['DRB1_04_04', 'DRB_W007_02_01', 'DRB_W003_07'],
    'DR09.02': ['DRB1_04_08', 'DRB_W007_01'],
    'DR10.01': ['DRB1_07_01', 'DRB3_04_05', 'DRB5_03_03'],
    'DR10.02': ['DRB1_07_01', 'DRB3_04_05', 'DRB5_03_01'],
    'DR11.01': ['DRB_W025_01'],
    'DR11.02': ['DRB_W205_w_01'],
    'DR11.03': ['DRB_W025_05', 'DRB1_07_04'],
    'DR13.01': ['DRB1_03_18', 'DRB_W006_03', 'DRB_W006_04'],
    'DR13.02': ['DRB1_03_18', 'DRB_W006_11', 'DRB_W006_04'],
    'DR14.01': ['DRB3_04_10', 'DRB_W004_02', 'DRB_W027_01'],
    'DR14.02': ['DRB3_04_10', 'DRB_W004_02', 'DRB_W027_02'],
    'DR15.01/02': ['DRB_W006_06', 'DRB_W021_04', 'DRB_W026g'],
    'DR15.03': ['DRB_W006_06', 'DRB_W021_04', 'DRB_W002_01'],
    'DR16.01': ['DRB1_03_10', 'DRB_W001_01', 'DRB_W006_02', 'DRB_W006_09_01'],
    'DR18.01': ['DRB4_01_02', 'DRB5_03_06'],
    'DR28.01': ['DRB1_07g', 'DRB4_01_04', 'DRB_W102_01'],
    'DR29.01': ['DRB1_10_11', 'DRB_W001_05'],
    'DR30.01': ['DRB1_07_05', 'DRB_W002_03']
}

indian_rhesus['MHC_DQA_HAPLOTYPES'] = {
    '01_02': ['DQA1_01_02'],
    '01_07': ['DQA1_01_07'],
    '01_09': ['DQA1_01_09'],
    '01g1': ['DQA1_01g1'],
    '01g2': ['DQA1_01g2'],
    '01g3': ['DQA1_01g3'],
    '01g4': ['DQA1_01g4'],
    '05_01': ['DQA1_05_01'],
    '05_02': ['DQA1_05_02'],
    '05_03': ['DQA1_05_03'],
    '05_04': ['DQA1_05_04'],
    '05_05': ['DQA1_05_05'],
    '05_06': ['DQA1_05_06'],
    '05_07': ['DQA1_05_07'],
    '23_01': ['DQA1_23_01'],
    '23_02': ['DQA1_23_02'],
    '23_03': ['DQA1_23_03'],
    '24_02': ['DQA1_24_02'],
    '24_04': ['DQA1_24_04'],
    '24_08': ['DQA1_24_08'],
    '24g1': ['DQA1_24g1'],
    '24g2': ['DQA1_24g2'],
    '26_01': ['DQA1_26_01'],
    '26g1': ['DQA1_26g1'],
    '26g2': ['DQA1_26g2']
}

indian_rhesus['MHC_DQB_HAPLOTYPES'] = {
    '06_01': ['DQB1_06_01'],
    '06_07': ['DQB1_06_07'],
    '06_08': ['DQB1_06_08'],
    '06_09': ['DQB1_06_09'],
    '06_10': ['DQB1_06_10'],
    '06_13_01': ['DQB1_06_13_01'],
    '06g1': ['DQB1_06g1'],
    '06g2': ['DQB1_06g2'],
    '06g3': ['DQB1_06g3'],
    '06g4': ['DQB1_06g4'],
    '15_02': ['DQB1_15_02'],
    '15g1': ['DQB1_15g1'],
    '15g2': ['DQB1_15g2'],
    '16_01': ['DQB1_16_01'],
    '16_03': ['DQB1_16_03'],
    '17_03': ['DQB1_17_03'],
    '17g1': ['DQB1_17g1'],
    '17g2': ['DQB1_17g2'],
    '17g3': ['DQB1_17g3'],
    '18_08': ['DQB1_18_08'],
    '18_10': ['DQB1_18_10'],
    '18_12': ['DQB1_18_12'],
    '18_17': ['DQB1_18_17'],
    '18_20': ['DQB1_18_20'],
    '18_24': ['DQB1_18_24'],
    '18g3': ['DQB1_18g3'],
    '18g4': ['DQB1_18g4'],
    '18g5': ['DQB1_18g5'],
    '24_01': ['DQB1_24_01'],
    '27g': ['DQB1_27g']
}

indian_rhesus['MHC_DPA_HAPLOTYPES'] = {
    '02_03': ['DPA1_02_03'],
    '02_08': ['DPA1_02_08'],
    '02_13': ['DPA1_02_13'],
    '02_14': ['DPA1_02_14'],
    '02_15': ['DPA1_02_15'],
    '02_16': ['DPA1_02_16'],
    '02_20': ['DPA1_02_20'],
    '02g1': ['DPA1_02g1'],
    '02g2': ['DPA1_02g2'],
    '02g3': ['DPA1_02g3'],
    '02g4': ['DPA1_02g4'],
    '04_01': ['DPA1_04_01'],
    '04_04': ['DPA1_04_04'],
    '04g': ['DPA1_04g'],
    '06g': ['DPA1_06g'],
    '07_01': ['DPA1_07_01'],
    '07_04': ['DPA1_07_04'],
    '07_09': ['DPA1_07_09'],
    '07g1': ['DPA1_07g1'],
    '07g2': ['DPA1_07g2'],
    '07g3': ['DPA1_07g3'],
    '08g': ['DPA1_08g'],
    '09_01': ['DPA1_09_01'],
    '10_01': ['DPA1_10_01'],
    '11_01': ['DPA1_11_01']
}

indian_rhesus['MHC_DPB_HAPLOTYPES'] = {
    '01g1': ['DPB1_01g1'],
    '01g2': ['DPB1_01g2'],
    '01g3': ['DPB1_01g3'],
    '01g4': ['DPB1_01g4'],
    '01g5': ['DPB1_01g5'],
    '02_02': ['DPB1_02_02'],
    '02g': ['DPB1_02g'],
    '03g': ['DPB1_03g'],
    '04_01': ['DPB1_04_01'],
    '05_01': ['DPB1_05_01'],
    '05_02': ['DPB1_05_02'],
    '06_04': ['DPB1_06_04'],
    '06g': ['DPB1_06g'],
    '07g1': ['DPB1_07g1'],
    '07g2': ['DPB1_07g2'],
    '08_01': ['DPB1_08_01'],
    '08_02': ['DPB1_08_02'],
    '15_03': ['DPB1_15_03'],
    '15g': ['DPB1_15g'],
    '16_01': ['DPB1_16_01'],
    '17_01': ['DPB1_17_01'],
    '18_01': ['DPB1_18_01'],
    '19_02': ['DPB1_19_02'],
    '19_06': ['DPB1_19_06'],
    '19g1': ['DPB1_19g1'],
    '19g2': ['DPB1_19g2'],
    '21_01': ['DPB1_21_01'],
    '21_02': ['DPB1_21_02'],
    '21_03': ['DPB1_21_03'],
    '23_01': ['DPB1_23_01'],
    '23_02': ['DPB1_23_02'],
    '24_01': ['DPB1_24_01']
}

### Mauritian cynomolgus macaques ###

mcm = {'PREFIX': 'Mafa'}

# MCM MHC updated by Roger 29 May 2018

mcm['MHC_A_HAPLOTYPES'] = {
    'M1A': ['05_M1M2M3_A1_063g', '07_M1M2_70_156bp', '11_M1_E_02g3|E_02_nov_09,_E_02_nov_10', '04_M1_AG_05_3mis_156bp'],
    'M2A': ['05_M1M2M3_A1_063g', '07_M1M2_70_156bp', '02_M2_G_02_06_156bp'],
    'M3A': ['05_M1M2M3_A1_063g', '07_M3_70_156bp'],
    'M4A': ['05_M4_A1_031_01'],
    'M5A': ['05_M5_A1_033_01'],
    'M6A': ['05_M6_A1_032_01', '05_M6_A1_047_01'],
    'M7A': ['05_M7_A1_060_05']
}

mcm['MHC_B_HAPLOTYPES'] = {
    'M1B': ['12_M1_B_134_02', '12_M1_B_152_01N'],
    'M2B': ['12_M2_B_019_03', '12_M2_B_150_01_01'],
    'M3B': ['12_M3_B_165_01', '12_M3_B_075_01'],
    'M4B': ['12_M4_B_088_01', '12_M4_B_127_nov_01'],
    'M5B': ['12_M5_B_167_01N', '12_M5_B_051_04'],
    'M6B': ['12_M6_B17_01_g103c', '12_M6_B_095_01'],
    'M7B': ['12_M7_B_072_02', '12_M7_B_166_01']
}

mcm['MHC_DRB_HAPLOTYPES'] = {
    'M1DR': ['13_M1_DRB_W21_01', '13_M1_DRB_W5_01'],
    'M2DR': ['13_M2_DRB1_10_01', '13_M2_DRB_W4_02'],
    'M3DR': ['13_M3_DRB1_10_02', '13_M3_DRB_W49_01_01'],
    'M4DR': ['13_M4_DRB4_01_01'],
    'M5DR': ['13_M5_DRB4_01_02'],
    'M6DR': ['13_M6_DRB1_04_02_01', '13_M6_DRB_W4_01'],
    'M7DR': ['13_M7_DRB_W1_03', '13_M7_DRB_W36_05']
}

mcm['MHC_DQA_HAPLOTYPES'] = {
    'M1DQ': ['14_M1_DQB1_18_01_01'],
    'M2DQ': ['14_M2_DQA1_01_04'],
    'M3DQ': ['14_M3_DQB1_16_01', '14_M3_DQA1_05_03_01'],
    'M4DQ': ['14_M4_DQB1_06_08', '14_M4_DQA1_01_07_01'],
    'M5DQ': ['14_M5_DQA1_01_06', '14_M5_DQB1_06_11'],
    'M6DQ': ['14_M6_DQA1_01_08_01'],
    'M7DQ': ['14_M7_DQA1_23_01', '14_M7_DQB1_18_14']
}

mcm['MHC_DQB_HAPLOTYPES'] = {
    'M1DQ': ['14_M1_DQB1_18_01_01'],
    'M2DQ': ['14_M2_DQA1_01_04'],
    'M3DQ': ['14_M3_DQB1_16_01', '14_M3_DQA1_05_03_01'],
    'M4DQ': ['14_M4_DQB1_06_08', '14_M4_DQA1_01_07_01'],
    'M5DQ': ['14_M5_DQA1_01_06', '14_M5_DQB1_06_11'],
    'M6DQ': ['14_M6_DQA1_01_08_01'],
    'M7DQ': ['14_M7_DQA1_23_01', '14_M7_DQB1_18_14']
}

mcm['MHC_DPA_HAPLOTYPES'] = {
    'M1DP': ['15_M1_DPA1_07_02', '15_M1_DPB1_19_03'],
    'M2DP': ['15_M2_DPA1_07_01', '15_M2_DPB1_20_01'],
    'M3DP': ['15_M3_DPB1_09_02'],
    'M4M7DP': ['15_M4M7_DPB1_03_03'],
    'M5M6DP': ['15_M5M6_DPB1_04_01']
}

mcm['MHC_DPB_HAPLOTYPES'] = {
    'M1DP': ['15_M1_DPA1_07_02', '15_M1_DPB1_19_03'],
    'M2DP': ['15_M2_DPA1_07_01', '15_M2_DPB1_20_01'],
    'M3DP': ['15_M3_DPB1_09_02'],
    'M4M7DP': ['15_M4M7_DPB1_03_03'],
    'M5M6DP': ['15_M5M6_DPB1_04_01']
}

haplotype_dict = {'MAMU': indian_rhesus, 'MCM': mcm, 'MANE': indian_rhesus}

"""## Genotype miSeq data against reference FASTA

This is a new implementation of the MHC genotyping pipeline. Considerations:

Current throughput is about 360 samples per hour (10 seconds per sample). 
- Export to Excel similar to current format
- Jupyter Notebook for portability and reproducible data analysis. This is really important so we can distribute users' data and the full analysis of their results. 

One possibly controversial decision in this algorithm is that I selectively include identical sequences that are found  as a fraction of total reads. This runs the risk of losing some sequences that could potentially be informative. When making this decision, I thought a lot about lossless compression of music. There is a lot of discussion about whether lossy compression of music files (e.g., 320kb MP3 is distinguishable from lossless FLAC/ALAC (https://www.npr.org/sections/therecord/2015/06/02/411473508/how-well-can-you-hear-audio-quality). I think there is a parallel in MHC genotyping -- do we really need to know all MHC if they are present in very low abundance of cDNA? Could we improve genotyping by simply reporting those sequences that comprise a significant fraction of reads (set to 0.1% of total reads by default)? I would need to be convinced that this really helps.

## Dependencies

+ Jupyter Notebook/Jupyter Lab
+ Python 3 (tested on anaconda distribution of Python 3.6.4)
+ Access to dholk.primate.wisc.edu
+ pigz (in PATH)
+ bbmap (in PATH)
+ bbmerge (in PATH)
+ bbduk (in PATH)
+ USEARCH v10 (attempts automatic installation if not available)
+ Pandas (tested from anaconda distribution)
"""

# generic functions


log = logging.getLogger(__name__)


def print_status(status):
    '''print timestamped status update'''
    print('--[' + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '] ' + status + '--')
    log.info(status)


def create_temp_folder():
    '''make temporary folder at specified location'''
    TMP_DIR = tempfile.mkdtemp(dir=os.path.join(output_dir, 'TMP'))
    return TMP_DIR


def close_temp_folder(tmp_dir):
    '''destroy temporary folder after it is no longer used'''
    os.removedirs(tmp_dir)


def create_output_folder(cwd):
    '''create timestamped output folder at specified location'''

    # fetch current time
    CURRENT_TIME = datetime.now().strftime("%Y%m%d%H%M%S")

    # path to output folder
    OUTPUT_FOLDER = cwd + '/' + CURRENT_TIME

    # create folder if it doesn't already exist
    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)

    # print output folder name
    print_status('Output folder: ' + OUTPUT_FOLDER)

    return OUTPUT_FOLDER


def run_command(cmd_list, stdout_file=None, stderr_file=None):
    '''run command with subprocess.call
    if stdout or stderr arguments are passed, save to specified file
    '''

    import subprocess

    print_status(' '.join(cmd_list))  # print status

    # if neither stdout or stderr specified
    if stdout_file is None and stderr_file is None:
        print(cmd_list)
        subprocess.call(cmd_list)

    # if only stdout is specified
    elif stdout_file is not None and stderr_file is None:
        with open(stdout_file, 'w') as so:
            subprocess.call(cmd_list, stdout=so)

    # if only stderr is specified
    elif stdout_file is None and stderr_file is not None:
        with open(stderr_file, 'w') as se:
            subprocess.call(cmd_list, stderr=se)

    # if both stdout and stderr are specified
    elif stdout_file is not None and stderr_file is not None:
        with open(stdout_file, 'w') as so:
            with open(stderr_file, 'w') as se:
                subprocess.call(cmd_list, stdout=so, stderr=se)

    else:
        pass


def test_executable(cmd):
    '''check that a particular command can be run as an executable'''

    import shutil

    assert shutil.which(cmd) is not None, 'Executable ' + cmd + ' cannot be run'


def get_notebook_path(out_dir):
    '''get name of  20835-genotyping.ipynb file in current working directory
    copy to output folder
    '''

    import os
    import shutil

    cwd = os.getcwd()  # get working directory
    notebook_path = cwd + '/26887-miseq-genotyping.ipynb'

    # copy to output folder
    shutil.copy2(notebook_path, out_dir + '/' + EXPERIMENT + '.ipynb')


def file_size(f):
    '''return file size'''
    import os

    return os.stat(f).st_size


"""## Create data structure for paired-end reads

Make dictionary containing R1/R2 miSeq read pairs, alternative sample identifiers, and miSeq run IDs and uses sample names as the dictionary key. This should only cause problems if multiple samples with the same name are run in the same workflow invokation, which seems unlikely. The system requires gzip-FASTQ files and will abort if FASTQ files are not gzip-compressed. This is to protect users from themselves - having uncompressed FASTQ files in the filesystem is a quick way to fill hard drives. There _is_ an uncompressed, merged FASTQ file that gets generated in the temporary intermediate files; such files are necessary for USEARCH functionality.

As a convenience, animal identifiers and run information is downloaded directly from the dholk LabKey server. Data from every genotyping run should be in this system so I have built this workflow to deliberately break if the data isn't in LabKey. Just today (2018-05-25 --dho) I had to spend time scouring BaseSpace for miSeq files that weren't properly archived in our LabKey system. So it does not seem unreasonable to me to enforce correct usage of the miSeq LabKey system.
"""


def is_gz_file(filepath):
    '''test if file is gzip-compressed
    return True if gzip-compressed
    source: https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed
    '''

    import binascii

    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'


def get_read_dict(df_samples_i):
    READS = {}
    sample_list = list(df_samples_i['Sample_ID'].unique())
    for sample_i in sample_list:

        df_samples_j = df_samples_i[df_samples_i['Sample_ID'] == sample_i]

        CLIENT_ID = list(df_samples_j['Sample_Name'])[0]
        RUN_ID = list(df_samples_j['Run'])[0]
        read_list = list(df_samples_j['FILEPATH'])

        for filepath_i in read_list:
            if 'R1' == filepath_i.split('_')[-2]:
                R1 = filepath_i
                break
        for filepath_i in read_list:
            if 'R2' == filepath_i.split('_')[-2]:
                R2 = filepath_i
                break
        comments = ""

        if "Comments" in df_samples_j.columns:
            df_samples_j['Comments'] = df_samples_j['Comments'].astype(str)
            comments = ';'.join(df_samples_j['Comments'].unique())
            if comments == 'nan':
                comments = ''
        READS[sample_i] = [R1, R2, CLIENT_ID, RUN_ID, comments]
    return READS


def decompress_fastq(reads):
    '''decompress FASTQ file with pigz and return path to decompressed file'''

    # ensure bbduk can be run from path
    test_executable('pigz')

    # command
    decompress_fastq_cmd = ['pigz',
                            '-d',
                            '-k',
                            reads]

    # run command
    run_command(decompress_fastq_cmd)

    # return path to decompressed file
    return os.path.splitext(reads)[0]


def vsearch_unique(reads, out_dir):
    # vsearch_path
    # test_executable(vsearch_path)

    READS_BN = [x for x in map(str.strip, (os.path.basename(reads)).split('.')) if x][0]

    vsearch_unique_cmd = [vsearch_path,
                          '--fastx_uniques',
                          reads,
                          '--sizeout',
                          '--relabel',
                          'Uniq',
                          '--fastaout',
                          out_dir + '/' + READS_BN + '.unique.fasta', ]

    # run command
    run_command(vsearch_unique_cmd)

    # test that output file exists before exiting function
    assert os.path.exists(
        out_dir + '/' + READS_BN + '.unique.fasta') == 1, out_dir + '/' + READS_BN + '.unique.fasta' + ' does not exist'

    # return unique sequence FASTQ and total number of merged reads
    return out_dir + '/' + READS_BN + '.unique.fasta'


def usearch_unique(reads, out_dir):
    '''expect merged gzip-compressed FASTQ files
    run clumpify'''

    import os

    # ensure bbduk can be run from path
    test_executable(USEARCH_PATH)

    # make sure FASTQ files exists
    assert os.path.exists(reads) == 1, 'FASTQ file ' + reads + ' does not exist'

    # get basename for reads
    READS_BN = [x for x in map(str.strip, (os.path.basename(reads)).split('.')) if x][0]

    # decompress FASTQ for use with USEARCH
    READS_DECOMPRESSED = decompress_fastq(reads)
    assert os.path.exists(READS_DECOMPRESSED) == 1, READS_DECOMPRESSED + ' does not exist'

    # USEARCH unique command

    usearch_unique_cmd = [USEARCH_PATH,
                          '-fastx_uniques',
                          READS_DECOMPRESSED,
                          '-sizeout',
                          '-relabel',
                          'Uniq',
                          '-fastaout',
                          out_dir + '/' + READS_BN + '.unique.fasta', ]

    # run command
    run_command(usearch_unique_cmd)

    # test that output file exists before exiting function
    assert os.path.exists(
        out_dir + '/' + READS_BN + '.unique.fasta') == 1, out_dir + '/' + READS_BN + '.unique.fasta' + ' does not exist'

    # return unique sequence FASTQ and total number of merged reads
    return (out_dir + '/' + READS_BN + '.unique.fasta')


def vsearch_denoise(reads, read_ct, out_dir):
    '''remove sequencing artifacts and chimeras by running UNOISE from the USEARCH package
    Expects decompressed FASTA unique sequences created from fastx_uniques command
    Preserves output sequences over a calculated minimum abundance
    '''

    import os

    # ensure bbduk can be run from path
    # test_executable(vsearch_path)

    # make sure FASTA file exists
    assert os.path.exists(reads) == 1, 'Unique FASTA file ' + reads + ' does not exist'

    # get basename for reads
    READS_BN = [x for x in map(str.strip, (os.path.basename(reads)).split('.')) if x][0]

    # calculate UNOISE minuniquesize threshold
    MIN_READS = str(calculate_threshold(min_freq=0.0002, reads=reads, read_ct=read_ct))
    print("MINREADS: {0}".format(MIN_READS))
    # if MIN_READS < 2 (as can happen for samples with low coverage, set MIN_READS = 2)

    if int(MIN_READS) < 2:
        MIN_READS = '2'

    # USEARCH unoise command
    # include only sequences greater than min_reads threshold

    vsearch_unoise_cmd = [vsearch_path,
                          '--cluster_unoise',
                          reads,
                          '--minsize',
                          MIN_READS,
                          '--unoise_alpha',
                          '2',
                          '--centroids',
                          os.path.join(out_dir, READS_BN + '.unoise.fasta')]

    # run command
    run_command(vsearch_unoise_cmd)

    # extract ZOTU sequences from unique FASTA
    vsearch_chimera_cmd = [vsearch_path,
                           '--uchime_denovo',
                           os.path.join(out_dir, READS_BN + '.unoise.fasta'),
                           '--abskew',
                           '16',
                           '--nonchimeras',
                           os.path.join(out_dir, READS_BN + '.zotu_descriptive.fasta')]

    # run command
    run_command(vsearch_chimera_cmd, stdout_file=os.path.join(out_dir, 'stdout.txt'),
                stderr_file=os.path.join(out_dir, 'stderr.txt'))

    # test that output file exists before exiting function
    assert os.path.exists(os.path.join(out_dir,
                                       READS_BN + '.zotu_descriptive.fasta')) == 1, out_dir + '/' + READS_BN + '.zotu_descriptive.fasta' + ' does not exist'

    # return unique sequence FASTQ and total number of merged reads
    return os.path.join(out_dir, READS_BN + '.zotu_descriptive.fasta')


def parse_unoise_output(denoise_stats, zotu_tmp_fasta, reads_bn, out_dir):
    '''UNOISE3 creates a new FASTA file without size annotations
    Parse zotu stats file to add size annotations to zotus so read count is preserved
    '''

    import re

    # make list of Zotu identifiers

    descriptive_names = []

    # read unoise_tabbed_output
    with open(denoise_stats) as fp:
        for line in fp:
            if 'zotu' in line:  # save only lines that correspond to zotus
                descriptive_names.append(line.split("\t")[0])  # save descriptive name, such as Uniq413;size=9;

    # read zotu_tmp_fasta
    # replace non-informative Zotu name with descriptive name
    ct = 0  # initialize counter
    with open(out_dir + '/' + reads_bn + '.zotu_descriptive.fasta', 'w') as w:  # save zotus to new file
        with open(zotu_tmp_fasta) as fp:
            for line in fp:
                if '>' in line:  # save only lines that correspond to zotus
                    w.write(re.sub('>Zotu[0-9]*', '>' + descriptive_names[ct], line))  # write descriptive fasta header
                    ct = ct + 1
                else:
                    w.write(line)  # write sequence lines as-is

    return out_dir + '/' + reads_bn + '.zotu_descriptive.fasta'


def count_reads(reads):
    '''count number of reads in a FASTQ/FASTA file
    use bbmap stats.sh
    return count as integer
    '''

    import subprocess

    # path to stats.sh

    # ensure stats.sh can be run from path
    test_executable(STATS_PATH)

    # make sure FASTQ input file exists
    assert os.path.exists(reads) == 1, 'Read file ' + reads + ' does not exist'

    # stats command
    stats_cmd = [STATS_PATH,
                 'in=' + reads,
                 'format=4']

    # run command
    output = subprocess.check_output(stats_cmd).decode("utf-8")

    # filter output
    # get first entry on second line, which will always be sequence length
    read_ct = (output.split('\n')[1]).split('\t')[0]

    return int(read_ct)


def calculate_threshold(min_freq, reads, read_ct):
    '''determine the -minuniquesize parameter for usearch unique
    count number of reads in FASTQ file and then return threshold
    0.001 would save unique sequences greater than 0.1% of all sequences in dataset'''

    # multiply read threshold by read count
    minuniquesize = int(min_freq * read_ct)

    # return threshold and number of total reads)
    return (minuniquesize)


"""## Remove primers with bbduk

Need to remove primer sequences from left and right ends of individual reads before merging. It is essential to do this in multiple steps (left-end trimming first; right end trimming second) to prevent removal of primers from incorrect locations within the sequences. The default k-mer size for bbduk is 50 so if this is not set to a number as small or smaller than the smallest PCR primer trimming will not work correctly. 
"""


def remove_primers(R1, R2, primers, out_dir):
    '''expect gzip-compressed R1 and R2 FASTQ files
    run bbduk to remove primers
    first remove from left end of sequence
    then remove from right end of sequence
    save results to temporary directory'''

    import os

    # path to bbduk

    # ensure bbduk can be run from path
    test_executable(BBDUK_PATH)

    # make sure primer file exists
    assert os.path.exists(primers) == 1, 'Specified primer file ' + primers + ' does not exist'

    # get basename for R1 and R2
    R1_BN = os.path.splitext(os.path.basename(R1))[0]
    R2_BN = os.path.splitext(os.path.basename(R2))[0]

    # trim from left
    left_primer_cmd = [BBDUK_PATH,
                       'in=' + R1,
                       'in2=' + R2,
                       'ref=' + "'" + primers + "'",
                       'ktrim=l',
                       'k=15',
                       'restrictleft=30',
                       'out=' + out_dir + '/' + R1_BN + '_l.fastq.gz',
                       'out2=' + out_dir + '/' + R2_BN + '_l.fastq.gz']

    # next trim from right
    right_primer_cmd = [BBDUK_PATH,
                        'in=' + out_dir + '/' + R1_BN + '_l.fastq.gz',
                        'in2=' + out_dir + '/' + R2_BN + '_l.fastq.gz',
                        'ref=' + "'" + primers + "'",
                        'ktrim=r',
                        'restrictright=30',
                        'k=15',
                        'out=' + out_dir + '/' + R1_BN + '_lr.fastq.gz',
                        'out2=' + out_dir + '/' + R2_BN + '_lr.fastq.gz']

    # run commands
    run_command(left_primer_cmd)

    # make sure output files from left primer trim exit
    assert os.path.exists(
        out_dir + '/' + R1_BN + '_l.fastq.gz') == 1, out_dir + '/' + R1_BN + '_l.fastq.gz' + ' does not exist'
    assert os.path.exists(
        out_dir + '/' + R2_BN + '_l.fastq.gz') == 1, out_dir + '/' + R2_BN + '_l.fastq.gz' + ' does not exist'

    run_command(right_primer_cmd)
    assert os.path.exists(
        out_dir + '/' + R1_BN + '_lr.fastq.gz') == 1, out_dir + '/' + R1_BN + '_lr.fastq.gz' + ' does not exist'
    assert os.path.exists(
        out_dir + '/' + R2_BN + '_lr.fastq.gz') == 1, out_dir + '/' + R2_BN + '_lr.fastq.gz' + ' does not exist'

    # return output files as tuple
    return (out_dir + '/' + R1_BN + '_lr.fastq.gz', out_dir + '/' + R2_BN + '_lr.fastq.gz')


"""## Merge reads

After trimming, merge overlapping reads and use these merged reads for miSeq genotyping. bbmerge defaults are inconsistent when merging reads from longer amplicons (such as the 280bp MHC class II amplicons) that do not have a long overlap. Performance can be improved by first merging with the default parameters but then repeating merging after applying 3' trims of increasing stringency. Roger and I empirically tested different trim qualities up to 40 and discovered that performance does not improve beyond a trip of 30. 
"""


def merge_reads(R1, R2, out_dir):
    '''expect gzip-compressed R1 and R2 FASTQ files
    run bbmerge'''

    import os

    # path to bbduk

    # ensure bbduk can be run from path
    test_executable(BBMERGE_PATH)

    # make sure R1 and R2 files exists
    assert os.path.exists(R1) == 1, 'R1 file ' + R1 + ' does not exist'
    assert os.path.exists(R2) == 1, 'R1 file ' + R2 + ' does not exist'

    # get basename for R1
    R1_BN = [x for x in map(str.strip, (os.path.basename(R1)).split('.')) if x][0]

    # merge read command
    merge_cmd = [BBMERGE_PATH,
                 'in=' + R1,
                 'in2=' + R2,
                 'qtrim2=r',
                 'trimq=10,15,20,25,30',
                 'pfilter=0.1',
                 'out=' + out_dir + '/' + R1_BN + '.merged.fastq.gz']

    # run command
    run_command(merge_cmd)

    # test that output file exists before exiting function
    assert os.path.exists(
        out_dir + '/' + R1_BN + '.merged.fastq.gz') == 1, out_dir + '/' + R1_BN + '.merged.fastq.gz' + ' does not exist'

    # return output merged FASTQ
    return out_dir + '/' + R1_BN + '.merged.fastq.gz'


"""## Find unique sequences and remove chimeras using USEARCH
A key performance enhancement in this workflow is recognizing that it is easier to map an amplicon sequence that occurs 5,000 times once, after adding a header indicating that there are 5,000 identical copies of the sequence, than mapping each of the 5,000 sequences individually. USEARCH fastx_uniques condenses a set of FASTQ sequences into its unique members. This tool requires decompressed FASTQ files, so the merged FASTQ files are decompressed prior to running USEARCH. To maximize performance, decompression uses the multithreaded pigz tool that should be installed in the user's path.

Another USEARCH tool, UNOISE3, purports to find "authentic" sequences by removing sequencing artifacts and chimeric sequences. Empirically, this tool performs quite well and using it to remove false positives contributes to improved genotyping accuracy and haplotying imputing. It is important to remove incomplete-length amplicon sequences before running UNOISE3 -- in the workflow, unique reads are first mapped to the reference sequence and any that are partial length are removed prior to running UNOISE3.

As described in the introduction, another innovation in this workflow is not reporting poorly supported genotypes below a specified threshold. This is calculated by determining the total number of merged reads in each sample and calculating a minimum read_abundance to report as UNOISE3 output. I initially set this to 0.1%, but Roger thought 0.02% of total reads is a more appropriate threshold. This could be adjusted in the future if there is a consensus this value is too high or too low.

## Map unique reads to reference

bbmap in semiperfect mode is used for read mapping. According to the documentation:

```
semiperfectmode=f       Allow only perfect and semiperfect (perfect except for 
                        N's in the reference) mappings.
```

This function is used for two different purposes. First, it is used to output FASTA files following mapping unique reads and before denoising to remove incomplete length amplicons. Second, it is used to output SAM files after mapping unique, denoised sequences. The SAM file is what is parsed to determine a sample's genotype. In both cases, ambiguous mappings are discarded since there should only be one reference sequence to which each read maps (or else the reference sequences themselves are ambiguous, which shouldn't happen).
"""


def map_semiperfect(reads, ref, out_dir, out_fmt):
    '''expect non-compressed FASTQ file of unique sequences
    run bbmap in semiperfect mode to map reads
    return mapped reads in specified format (out_fmt)
    this is needed because mapped reads are returned as .fasta when removing partial-length mathc
    and .sam when reporting allele calls
    Also need to set ordered=t to keep reads in size-descending order for UNOISE3 to work correctly
    '''

    import os

    # path to bbmap

    # ensure bbmap can be run from path
    test_executable(BBMAP_PATH)

    # make sure FASTQ input file and reference file exists
    assert os.path.exists(reads) == 1, 'Unique FASTA file ' + reads + ' does not exist'
    assert os.path.exists(ref) == 1, 'Reference FASTA file ' + ref + ' does not exist'

    # get basename for reads
    READS_BN = [x for x in map(str.strip, (os.path.basename(reads)).split('.')) if x][0]

    # bbmap command
    # do not make SAM header because this interferes with pandas parsing
    # toss all ambiguously mapped reads since the reference database should be unambiguous
    bbmap_cmd = [BBMAP_PATH,
                 'in=' + reads,
                 'outm=' + out_dir + '/' + READS_BN + '.' + out_fmt,
                 'noheader=t',
                 'nodisk=t',
                 'ambiguous=toss',
                 'ref=' + "'" + ref + "'",
                 'ordered=t',
                 'semiperfectmode=t']

    # run command
    run_command(bbmap_cmd)

    # test that output file exists before exiting function
    assert os.path.exists(
        out_dir + '/' + READS_BN + '.' + out_fmt) == 1, out_dir + '/' + READS_BN + '.' + out_fmt + ' does not exist'

    # return BAM file of mapped reads
    return out_dir + '/' + READS_BN + '.' + out_fmt


"""## Parse SAM file

Use Pandas to parse the SAM output from bbmap. This is necessary to extract the columns necessary for genotyping and to aggregate read counts from multiple unique FASTQ sequences that map uniquely to a single reference target. This can happen when there are different soft clips outside of the reference sequence.

The genotypes are used to determine which haplotypes are present in a sample. Each locus is considered individually against a dictionary of haplotype/genotype definitions. It is designed to be conservative -- if the algorithm cannot confidently assign a haplotype, it will report an error instead of guessing, possibly incorrectly. It is better to have manual intervention than incorrect automatic haplotype assignments.
"""


def parse_sam(mapped_sam, sample_id, total_read_ct, client_id, experiment, run_id, species, comments):
    '''parse SAM output file to determine number of reads per reference
    required because bbmap nor usearch respect read abundance when mapping
    '''

    # read SAM file
    df = pd.read_csv(mapped_sam, sep='\t', header=None)

    # extract read count and allele columns
    genotyping_df = df[[0, 2]]

    # rename columns
    genotyping_df = genotyping_df.rename(columns={0: 'read_ct', 2: 'allele'})

    # extract size value from first column
    genotyping_df['read_ct'] = genotyping_df['read_ct'].str.replace('Uniq[0-9]*;size=', '')
    genotyping_df['read_ct'] = genotyping_df['read_ct'].str.replace(';', '')

    # convert count column to numeric
    genotyping_df[['read_ct']] = genotyping_df[['read_ct']].apply(pd.to_numeric)

    # get total number of mapped reads
    mapped_reads = genotyping_df['read_ct'].sum()

    # group by allele and aggregate read_ct for identical alleles
    # likely due to differently trimmed unique sequences mapping to the same reference
    genotyping_df = genotyping_df.groupby(['allele']).agg({'read_ct': 'sum'})

    # reset index
    genotyping_df = genotyping_df.reset_index()

    # add sample_id as column
    genotyping_df['gs_id'] = sample_id

    # add mapped_read_ct as column
    genotyping_df['mapped_read_count'] = mapped_reads

    # add total_read_ct as column
    genotyping_df['total_read_count'] = total_read_ct

    # add percent_unmapped as column
    genotyping_df['percent_reads_unmapped'] = round(100 - (mapped_reads / total_read_ct * 100), 1)

    # add client_id as column
    genotyping_df['client_id'] = client_id

    # add experiment as column
    genotyping_df['experiment'] = experiment

    # add run_id as column
    genotyping_df['run_id'] = run_id

    # evaluate haplotypes
    # MHC-A
    MHC_A = call_haplotypes(locus=str(species['PREFIX'] + '-A'),
                            locus_haplotype_definitions=species['MHC_A_HAPLOTYPES'], df=genotyping_df)
    genotyping_df['MHC-A Haplotype 1'] = MHC_A[0]
    genotyping_df['MHC-A Haplotype 2'] = MHC_A[1]

    # MHC-B
    MHC_B = call_haplotypes(locus=str(species['PREFIX'] + '-B'),
                            locus_haplotype_definitions=species['MHC_B_HAPLOTYPES'], df=genotyping_df)
    genotyping_df['MHC-B Haplotype 1'] = MHC_B[0]
    genotyping_df['MHC-B Haplotype 2'] = MHC_B[1]

    # MHC-DRB
    MHC_DRB = call_haplotypes(locus=str(species['PREFIX'] + '-DRB'),
                              locus_haplotype_definitions=species['MHC_DRB_HAPLOTYPES'], df=genotyping_df)
    genotyping_df['MHC-DRB Haplotype 1'] = MHC_DRB[0]
    genotyping_df['MHC-DRB Haplotype 2'] = MHC_DRB[1]

    # MHC-DQA
    MHC_DQA = call_haplotypes(locus=str(species['PREFIX'] + '-DQA'),
                              locus_haplotype_definitions=species['MHC_DQA_HAPLOTYPES'], df=genotyping_df)
    genotyping_df['MHC-DQA Haplotype 1'] = MHC_DQA[0]
    genotyping_df['MHC-DQA Haplotype 2'] = MHC_DQA[1]

    # MHC-DQB
    MHC_DQB = call_haplotypes(locus=str(species['PREFIX'] + '-DQB'),
                              locus_haplotype_definitions=species['MHC_DQB_HAPLOTYPES'], df=genotyping_df)
    genotyping_df['MHC-DQB Haplotype 1'] = MHC_DQB[0]
    genotyping_df['MHC-DQB Haplotype 2'] = MHC_DQB[1]

    # MHC-DPA
    MHC_DPA = call_haplotypes(locus=str(species['PREFIX'] + '-DPA'),
                              locus_haplotype_definitions=species['MHC_DPA_HAPLOTYPES'], df=genotyping_df)
    genotyping_df['MHC-DPA Haplotype 1'] = MHC_DPA[0]
    genotyping_df['MHC-DPA Haplotype 2'] = MHC_DPA[1]

    # MHC-DPB
    MHC_DPB = call_haplotypes(locus=str(species['PREFIX'] + '-DPB'),
                              locus_haplotype_definitions=species['MHC_DPB_HAPLOTYPES'], df=genotyping_df)
    genotyping_df['MHC-DPB Haplotype 1'] = MHC_DPB[0]
    genotyping_df['MHC-DPB Haplotype 2'] = MHC_DPB[1]

    # add comments field
    genotyping_df['Comments'] = comments

    return genotyping_df


def call_haplotypes(locus, locus_haplotype_definitions, df):
    '''specify locus (e.g., Mamu-A) to haplotype and provide dictionary of haplotype definitions
    update specified pandas dataframe
    '''
    # convert alleles in dataframe to string
    # makes it possible to search for values more easily
    allele_str = df['allele'].to_string(header=False, index=False)

    # create list to store haplotypes for a sample
    sample_haplotypes = []

    # loop through haplotypes
    for haplotype, alleles in locus_haplotype_definitions.items():

        # if all diagnostic alleles for a haplotype are present
        # save haplotype name to list
        if all((x) in allele_str for x in alleles):
            sample_haplotypes.append(haplotype)

        # evaluate haplotypes
        if len(sample_haplotypes) == 0:  # if there are no haplotypes, that isn't possible
            # special exception for MCM A1*063 which is often undercalled but is very important
            if locus == 'Mafa-A' and '05_M1M2M3_A1_063g' in allele_str:
                h1 = 'A1_063'
                h2 = '-'
            else:
                h1 = 'ERR: NO HAP'
                h2 = 'ERR: NO HAP'
        elif len(sample_haplotypes) == 1:
            # special exception for MCM A1*063 which is often undercalled but is very important
            # if A1_063 is present in genotypes but M1, M2, and M3 are not present in single called haplotype,
            # add A1_063 to haplotype2
            if locus == 'Mafa-A' and '05_M1M2M3_A1_063g' in allele_str and not any(
                    y in sample_haplotypes[0] for y in ('M1A', 'M2A', 'M3A')):
                h1 = sample_haplotypes[0]
                h2 = 'A1_063'
            else:
                h1 = sample_haplotypes[0]
                h2 = '-'
        elif len(sample_haplotypes) == 2:
            h1 = sample_haplotypes[0]
            h2 = sample_haplotypes[1]
        elif len(sample_haplotypes) > 2:
            h1 = 'ERR: TMH (' + ', '.join(sample_haplotypes) + ')'
            h2 = 'ERR: TMH (' + ', '.join(sample_haplotypes) + ')'

    # DPA, DPB, DQA, DQB can only have two genotypes though other loci can have more
    # for these loci, error if more than two genotypes are reported

    # test number of rows that match locus, if > 2 set h1 and h2 to error
    # need to explicitly cast allele as string
    diploid_loci = ['DPA', 'DPB', 'DQA', 'DQB']

    for i in diploid_loci:
        if (i in locus):
            if df.allele.astype(str).str.contains(i).sum() > 2:
                h1 = 'ERR: TMG'
                h2 = 'ERR: TMG'

    return (h1, h2)


"""## Data reporting

A major challenge of the current system is how it takes to generate a finalized report from genotyping data. 
By leveraging integration with LabKey and Pandas, hopefully this can be improved. 
The cell below specifies haplotype defintiions in an editable form. 
Because this entire file is included in the output folder, 
it provides an unambiguous way of showing which definition set is used for which analysis.

The second cell takes a columnal list of genotypes and converts it to a PivotTable int he format that genotypes and 
haplotypes are typically reported to Genetics Services clients. 
"""


def pivot_pandas(df):
    """
    Create pivot table from genotyping Pandas dataframe
    """
    # add locus identifier column to genotypes

    df['locus'] = df['allele'].str.replace('_.*', '')

    # pivot data
    DF_PIVOT = pd.pivot_table(df,
                              index=['locus', 'allele'],
                              columns=[
                                  'gs_id',
                                  'client_id',
                                  'mapped_read_count',
                                  'total_read_count',
                                  'percent_reads_unmapped',
                                  'MHC-A Haplotype 1',
                                  'MHC-A Haplotype 2',
                                  'MHC-B Haplotype 1',
                                  'MHC-B Haplotype 2',
                                  'MHC-DRB Haplotype 1',
                                  'MHC-DRB Haplotype 2',
                                  'MHC-DQA Haplotype 1',
                                  'MHC-DQA Haplotype 2',
                                  'MHC-DQB Haplotype 1',
                                  'MHC-DQB Haplotype 2',
                                  'MHC-DPA Haplotype 1',
                                  'MHC-DPA Haplotype 2',
                                  'MHC-DPB Haplotype 1',
                                  'MHC-DPB Haplotype 2',
                                  'Comments',
                                  'experiment',
                                  'run_id'],
                              values='read_ct',
                              aggfunc=sum)

    # to manipulate the pandas dataframe with multiindexes a lot of steps are involved
    # we want to:
    # break apart the database allele name into the allele group name and ambiguous individual alleles
    # add a column showing how many times each allele occurs in a dataset
    # add a column showing how many reads support each allele call in a dataset
    # reorder the columns to match the desired Excel display

    # create list of sample names
    # this will be needed when reordering columns

    SAMPLE_COLUMNS = []  # list to store data columns
    for count, i in enumerate(DF_PIVOT.columns):  # iterate over data columns
        SAMPLE_COLUMNS.append(DF_PIVOT.columns[count][0])  # add sample name to list

    # extract database allele names to split into allele groups and groups of ambiguous sequences comprising groups
    t = DF_PIVOT.index.values  # get all values in dataframe, yields list where allele is in position [1]
    l = []  # initalize list to store allele
    for i in t:  # iterate over values
        l.append(i[1])  # save allele name to list

    # create column with allele groupings
    DF_PIVOT['allele_group'] = l  # create new column containing full allele sequences from database
    DF_PIVOT['allele_group'] = DF_PIVOT['allele_group'].astype(str).str.replace('\|.*',
                                                                                '')
    # remove ambiguous allele names to yield Mamu_A1_004g

    # create column with ambiguous alleles for each grouping
    DF_PIVOT['ambiguous_alleles'] = l  # create new column containing ambiguous alleles
    DF_PIVOT['ambiguous_alleles'] = DF_PIVOT['ambiguous_alleles'].astype(str).str.replace('.*\|',
                                                                                          '')
    # remove allele grouping information to yield A1_004_01_01,_
    # A1_004_01_02,_A1_004_02_01,_A1_004_02_02,_A1_004_05,_A1_004_06

    # get number of data columns
    TOTAL_COLUMNS = (len(DF_PIVOT.columns))  # count number of data columns

    # select columns containing genotype data
    # operate specifically on these columns in the next steps
    DATA_COLUMNS = DF_PIVOT.iloc[:, 0:TOTAL_COLUMNS]

    # get number of NaN values per genotype
    NULL_COLUMNS = DATA_COLUMNS.isnull().sum(axis=1).tolist()  # count number of NaN columns

    # subtract number of null columns from the number of data columns to get number of columns with genotype
    DF_PIVOT['obs_count'] = NULL_COLUMNS
    DF_PIVOT['obs_count'] = TOTAL_COLUMNS - DF_PIVOT['obs_count']

    # sum number of reads per genotype
    DF_PIVOT['sum_genotype_read_ct'] = DATA_COLUMNS.sum(axis=1).tolist()

    # flatten dataframe
    # this removes the allele index
    DF_PIVOT = DF_PIVOT.reset_index(drop=True)

    # create list with columns as desired in output Excel file
    REORDERED_COLUMNS = ['allele_group',
                         'sum_genotype_read_ct',
                         'obs_count',
                         ]

    for i in SAMPLE_COLUMNS:  # iterate over sample columns
        REORDERED_COLUMNS.append(i)  # add samples to ordered columns

    REORDERED_COLUMNS.append('ambiguous_alleles')

    DF_PIVOT = DF_PIVOT[REORDERED_COLUMNS]  # change column order

    return DF_PIVOT




def generate_excel_report(df, out_dir):
    '''make report format in Excel comparable to what GS currently uses
    though I'd rather see the pandas dataframe loaded into LabKey directly,
    generating an Excel report directly is the best way to troubleshoot this program. After people
    accept the Excel report, we can load the Pandas data into LabKey and then run the Excel reporting
    function to generate the same report from LabKey data.

    Try to decompose pandas dataframe to a series of lists for more flexible control over Excel formatting
    '''



    ## Configure Excel Report ##
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    workbook = xlsxwriter.Workbook(out_dir + '/pivot.xlsx')
    worksheet = workbook.add_worksheet()

    # freeze top 20 rows and left three columns
    worksheet.freeze_panes(22, 3)

    # add decimal formatting to columns with signals
    worksheet.set_column('A:A', 50)

    # get dataframe as list
    df_header = df.columns.tolist()

    ## create Excel formats ##
    # add header_id format that will be applied to top rows of table (rows 1-22)
    header_format = workbook.add_format({
        'bold': True,
        'align': 'left',
        'border': 0})

    # format for header number fields (e.g. mapped read_count)
    # uses thousands separator to make these fields easier to read
    header_number_format = workbook.add_format({
        'bold': True,
        'align': 'left',
        'border': 0,
        'num_format': '#,###'})

    if 'Mafa' in SPECIES.values():
        # conditional formatting for MCM haplotype coloring
        # this is pretty painful to read - sorry for this
        # the key of mcm_haplotype_formats is the haplotype identifier - it is also used as a search term for conditional formatting
        # the value of mcm_haplotype_formats contains a list of background colors [0] and font colors [1] for haplotype highlighting
        # these background values are also used to highlight alleles:
        # M5 has a yellow background color which is almost impossible to read as a font color, so a secondary font color is used in this case

        mcm_haplotype_formats = {'M1': ['#0c0000', '#ffffff'],
                                 'M2': ['#FF0000', '#ffffff'],
                                 'M3': ['#0000FF', '#ffffff'],
                                 'M4': ['#008000', '#ffffff'],
                                 'M5': ['#FFFF00', '#0c0000', '#b2b200'],
                                 'M6': ['#808080', '#ffffff'],
                                 'M7': ['#800080', '#ffffff']}

        # create conditional formats for the 'header' portion of the Excel file that contains per-sample summary statistics and genotype information
        for key, value in mcm_haplotype_formats.items():
            haplotype_highlight_format = workbook.add_format({'bg_color': value[0], 'font_color': value[1]})

            worksheet.conditional_format('D6:ZZ19', {'type': 'text',
                                                     'criteria': 'containing',
                                                     'value': key,
                                                     'format': haplotype_highlight_format})

            # create conditional formats for the 'allele' column of the Excel file (the first column)

        # two types of conditional formats:
        # 1. Color cells containing diagnostic genotypes
        # This requires a separate conditional format for every diagnostic genotype
        # 2. Change font colors of all alleles corresponding to each haplotype
        # This requires changing colors based on matching the 'MX' (e.g., 'M1') keys from mcm_haplotype_formats
        # first extract all values from nested dictionaries containi
        # from https://stackoverflow.com/questions/5164642/python-print-a-generator-expression
        def NestedDictValues(d):
            for v in d.values():
                if isinstance(v, dict):
                    yield from NestedDictValues(v)
                else:
                    yield v

        # save haplotype list values
        haplotype_value_list = (list(NestedDictValues(mcm)))

        # create conditional formats
        # note that this will only apply to first 2000 rows
        # this should be enough but using much larger values (e.g., 9999) caused Excel to break
        for i in haplotype_value_list:
            if 'Mafa' not in i:  # ignore first list value
                for j in i:  # unpack list of diagnostic genotypes for each haplotype
                    for key, value in mcm_haplotype_formats.items():  # create formats for each haplotype
                        if key in j:  # only create highlights for genotypes that correspond to haplotypes
                            diagnostic_highlight_format = workbook.add_format(
                                {'bg_color': value[0], 'font_color': value[1]})

                            worksheet.conditional_format('A23:A2000', {'type': 'text',
                                                                       'criteria': 'containing',
                                                                       'value': j,
                                                                       'format': diagnostic_highlight_format,
                                                                       'stop_if_true': True})

        for key, value in mcm_haplotype_formats.items():
            if key == 'M5':  # use secondary, darker color for M5
                genotype_highlight_format = workbook.add_format({'font_color': value[2]})
            else:
                genotype_highlight_format = workbook.add_format({'font_color': value[0]})

                worksheet.conditional_format('A23:A2000', {'type': 'text',
                                                           'criteria': 'containing',
                                                           'value': key,
                                                           'format': genotype_highlight_format})

                ## write headers  - per sample information above data ##

    # write header identifiers to first column
    header_id = ['gs_id',
                 'client_id',
                 'mapped_read_count',
                 'total_read_count',
                 'percent_reads_unmapped',
                 'MHC-A Haplotype 1',
                 'MHC-A Haplotype 2',
                 'MHC-B Haplotype 1',
                 'MHC-B Haplotype 2',
                 'MHC-DRB Haplotype 1',
                 'MHC-DRB Haplotype 2',
                 'MHC-DQA Haplotype 1',
                 'MHC-DQA Haplotype 2',
                 'MHC-DQB Haplotype 1',
                 'MHC-DQB Haplodtype 2',
                 'MHC-DPA Haplotype 1',
                 'MHC-DPA Haplotype 2',
                 'MHC-DPB Haplotype 1',
                 'MHC-DPB Haplotype 2',
                 'Comments',
                 'experiment',
                 'run_id']

    # write header id column

    for row, i in enumerate(header_id):
        worksheet.write(row, 0, i, header_format)

    # write header columns

    for col, k in enumerate(df_header):  # loop over header columns

        # ignore allele_group column (first columnn in df_header) to put allele_group data underneath other header identifiers
        if col > 0:
            worksheet.write(0, col, k[0], header_format)  # gs_id
            worksheet.write(1, col, k[1], header_format)  # client_id
            worksheet.write(2, col, k[2], header_number_format)  # mapped_read_count
            worksheet.write(3, col, k[3], header_number_format)  # total_read_count
            worksheet.write(4, col, k[4], header_format)  # percent_reads_unmapped
            worksheet.write(5, col, k[5], header_format)  # MHC-A Haplotype 1
            worksheet.write(6, col, k[6], header_format)  # MHC-A Haplotype 2
            worksheet.write(7, col, k[7], header_format)  # MHC-B Haplotype 1
            worksheet.write(8, col, k[8], header_format)  # MHC-B Haplotype 2
            worksheet.write(9, col, k[9], header_format)  # MHC-DRB Haplotype 1
            worksheet.write(10, col, k[10], header_format)  # MHC-DRB Haplotype 2
            worksheet.write(11, col, k[11], header_format)  # MHC-DQA Haplotype 1
            worksheet.write(12, col, k[12], header_format)  # MHC-DQA Haplotype 2
            worksheet.write(13, col, k[13], header_format)  # MHC-DQB Haplotype 1
            worksheet.write(14, col, k[14], header_format)  # MHC-DQB Haplotype 2
            worksheet.write(15, col, k[15], header_format)  # MHC-DPA Haplotype 1
            worksheet.write(16, col, k[16], header_format)  # MHC-DPA Haplotype 2
            worksheet.write(17, col, k[17], header_format)  # MHC-DPB Haplotype 1
            worksheet.write(18, col, k[18], header_format)  # MHC-DPB Haplotype 2
            worksheet.write(19, col, k[19], header_format)  # Comments
            worksheet.write(20, col, k[20], header_format)  # experiment
            worksheet.write(21, col, k[21], header_format)  # run_id

    ## write data ##

    # remove NaN values from dataframe
    remove_nan = df.replace(np.nan, '', regex=True)

    # iterate over each dataframe column
    for current_column, column in enumerate(remove_nan):
        data = (remove_nan[column].tolist())  # convert dataframe column to list

        # iterate over each row and write data
        for idx, i in enumerate(data):
            # start printing data at row 23 in excel file
            starting_row = 22  # 0-based
            current_row = idx + starting_row
            worksheet.write(current_row, current_column, i)

    workbook.close()


"""# Workflow

This is the data processing workflow that specifies that order of the commands. It also copies files from the analysis to the output folder so there is a bundle of all the files necessary to reprecreate the analysis.

## Get meta data from SampleSheet.csv file
1.   List item
2.   List item

"""

# Sample_ID	Sample_Name	Sample_Project	Description	Species

lk_species_dict = {'PIG_TAILED': 'MANE',
                   'CYNOMOLGUS': 'MCM',
                   'RHESUS': 'MAMU',
                   'MAMU': 'MAMU',
                   'MCM': 'MCM',
                   'MANE': 'MANE'}
lk_species_list = list(lk_species_dict.keys())

# need an open raw
# scan for

j = 0
with open(sample_path) as fh:
    for line in fh:
        j += 1
        if line.startswith("[Data]"):
            break

df_samples = pd.read_csv(sample_path, skiprows=j)

# df_samples = pd.read_csv(sample_path)
# optional filter list set project_list to [] to allow all
if len(project_list) > 0:
    df_samples = df_samples[df_samples[project_column].isin(project_list)]

if species_column not in df_samples.columns:
    raise ('The sample sheet must have a "Species" Column')

df_samples[species_column] = df_samples[species_column].str.upper()
# This is where it gets the sample name from the sheet and renames it to merge the species at project to keep different Species, separate
df_samples['PROJECT_NAME'] = ['{0}_{1}'.format(y, lk_species_dict[x]) for x, y in
                              zip(df_samples[species_column], df_samples[project_column])]

# This is the path to the reference this file will use
df_samples['REF_PATH'] = [ref_dict['REF'][lk_species_dict[x]] for x in df_samples[species_column]]
# This is the path to the primers this file will use it is currently hardcoded to use the same primers for all samples and species.
df_samples['PRIMER_PATH'] = ref_dict['PCR_PRIMERS']
#
# # Get the lis of fastq.gz (or fq.gz) files

file_list = os.listdir(input_dir)
fastq_list = [x for x in file_list if (x.endswith('.fastq.gz') or x.endswith('.fq.gz')) and (
            (len(x.split('_')) > 2) & (x.split('_')[-2] in ['R1', 'R2']))]
fastq_list = [x for x in file_list if (x.endswith('.fastq.gz') or x.endswith('.fq.gz'))]
fastq_list.sort()


df_file = pd.DataFrame()
for fastq_i in fastq_list:
    sample_id = fastq_i.split('_S')[0]
    new_fastq_path = os.path.join(input_dir, fastq_i)
    df_file_i = pd.DataFrame({'Sample_ID': [sample_id],
                              'FILEPATH': [new_fastq_path]})
    df_file = pd.concat([df_file, df_file_i], ignore_index=True)
df_samples = df_samples.merge(df_file, on='Sample_ID')
df_samples['Run'] = RUN_ID

"""## Main Pipeline"""

# %%time
# need a way to log the output? And to new files.


# Workflow
print('start Pipeline')

now = datetime.now()
date_time = now.strftime('%Y_%m_%d__%H_%M')
# Setup run pick up all the sample projectes so they are given in separate pivot tables/result folders
# It also separate by species if the project has multiple species.
project_species_list = list(df_samples['PROJECT_NAME'].unique())

for project_i in project_species_list:
    # Clear output so  customers don't see the outputs from other customers
    # clear_output()
    print('start Pipeline for project: {0}'.format(project_i))
    df_samples_i = df_samples[df_samples['PROJECT_NAME'] == project_i]
    REF = list(df_samples_i['REF_PATH'])[0]
    PCR_PRIMERS = list(df_samples_i['PRIMER_PATH'])[0]
    READS = get_read_dict(df_samples_i)
    # create temporary folder
    TMP_DIR = os.path.join(output_dir, 'tmp', project_i)
    # os.makedirs(OUTPUT_FOLDER + '/tmp', exist_ok=True)
    os.makedirs(TMP_DIR, exist_ok=True)

    species_name = list(df_samples_i[species_column])[0]
    SPECIES = haplotype_dict[lk_species_dict[species_name]]
    print(SPECIES)

    ### make dataframe to store all genotypes
    ALL_GENOTYPES_DF = pd.DataFrame(columns=['read_ct',
                                             'allele',
                                             'gs_id',
                                             'mapped_read_count',
                                             'total_read_count',
                                             'percent_reads_unampped',
                                             'client_id',
                                             'MHC-A Haplotype 1',
                                             'MHC-A Haplotype 2',
                                             'MHC-B Haplotype 1',
                                             'MHC-B Haplotype 2',
                                             'MHC-DRB Haplotype 1',
                                             'MHC-DRB Haplotype 2',
                                             'MHC-DQA Haplotype 1',
                                             'MHC-DQA Haplotype 2',
                                             'MHC-DQB Haplotype 1',
                                             'MHC-DQB Haplotype 2',
                                             'MHC-DPA Haplotype 1',
                                             'MHC-DPA Haplotype 2',
                                             'MHC-DPB Haplotype 1',
                                             'MHC-DPB Haplotype 2',
                                             'Comments',
                                             'experiment',
                                             'run_id'
                                             ])

    ### Create folder for output files

    # DRIVE_OUTPUT_FOLDER = os.path.join(genotyper_root_dir, 'output',folder_name, '{0}_{1}'.format(project_i, date_time))
    OUTPUT_FOLDER = create_output_folder(os.path.join(output_dir, '{0}_{1}'.format(project_i, date_time)) + '/')

    os.makedirs(OUTPUT_FOLDER + '/SAM', exist_ok=True)
    os.makedirs(OUTPUT_FOLDER + '/ref', exist_ok=True)
    os.makedirs(OUTPUT_FOLDER + '/dev', exist_ok=True)

    #### Copy reference files to output directory
    shutil.copy2(PCR_PRIMERS, OUTPUT_FOLDER + '/ref')
    shutil.copy2(REF, OUTPUT_FOLDER + '/ref')

    ## Genotype individual samples

    ct = 1  # set counter

    for sample_id, reads in READS.items():

        print('--Processing sample ' + str(ct) + ' of ' + str(len(READS)) + '--')
        print('Sample ID: ' + sample_id)
        print('R1_FASTQ: ' + reads[0])
        print('R2_FASTQ: ' + reads[1])
        comments = reads[4]

        ct = ct + 1
        ### trim primers
        PRIMER_TRIMMED = remove_primers(R1=reads[0], R2=reads[1], primers=PCR_PRIMERS, out_dir=TMP_DIR)

        ### merge reads
        MERGED = merge_reads(R1=PRIMER_TRIMMED[0], R2=PRIMER_TRIMMED[1], out_dir=TMP_DIR)
        if file_size(MERGED) == 0:
            print(sample_id + ' merged FASTQ file is empty')
            continue

        ### count merged reads
        TOTAL_READ_CT = count_reads(MERGED)

        ### cluster identical reads
        UNIQUE = vsearch_unique(reads=MERGED, out_dir=TMP_DIR)
        if file_size(UNIQUE) == 0:
            print(sample_id + ' unique FASTQ file is empty')
            continue

        ## map reads and save potential full-length reads
        ## this avoids partial-matches obscuring full-length matches following UNOISE3
        REMOVE_PARTIAL = map_semiperfect(reads=UNIQUE, ref=REF, out_dir=TMP_DIR, out_fmt='fasta')
        if file_size(REMOVE_PARTIAL) == 0:
            print(sample_id + ' mapped FASTA file is empty')
            continue

        ### denoise artifacts and chimeric sequences
        ZOTU = vsearch_denoise(reads=REMOVE_PARTIAL, read_ct=TOTAL_READ_CT, out_dir=TMP_DIR)
        if file_size(ZOTU) == 0:
            print(sample_id + ' denoised FASTQ file is empty')
            continue

        ### map reads and copy SAM to output folder
        MAPPED_SAM = map_semiperfect(reads=ZOTU, ref=REF, out_dir=TMP_DIR, out_fmt='sam')
        if file_size(MAPPED_SAM) == 0:
            print(sample_id + ' mapped SAM file is empty')
            continue

        shutil.copy2(MAPPED_SAM, OUTPUT_FOLDER + '/SAM')

        ### parse sam
        # reads[2] is client_id extracted from LabKey
        GENOTYPES = parse_sam(MAPPED_SAM, sample_id, TOTAL_READ_CT, reads[2], EXPERIMENT, reads[3], SPECIES, comments)

        ### append genotypes to all_genotypes_df
        ALL_GENOTYPES_DF = pd.concat([ALL_GENOTYPES_DF, GENOTYPES])
        print(ct)

    ## Copy developer files to output folder - eventually make a flag to only display if necessary
    ALL_GENOTYPES_DF.to_csv(OUTPUT_FOLDER + '/dev/genotypes.csv')

    ## Report genotypes

    ### Make pivot table in pandas
    if len(ALL_GENOTYPES_DF) > 0:
        GENOTYPES_PIVOTED = pivot_pandas(df=ALL_GENOTYPES_DF)

        ### Write Excel report
        generate_excel_report(df=GENOTYPES_PIVOTED, out_dir=OUTPUT_FOLDER)
    else:
        print('All samples failed generate enough reads, no pivot created')
        Path(os.path.join(OUTPUT_FOLDER, 'pivot_is_blank.txt')).touch()
    # Bake up files
    # os.makedirs(DRIVE_OUTPUT_FOLDER, exist_ok=True)
    shutil.copy2(os.path.join(genotyper_root_dir, 'main.py'), os.path.join(OUTPUT_FOLDER, 'main.py'))
    shutil.rmtree(TMP_DIR, ignore_errors=True)
    print('--PIPELINE COMPLETE--')
