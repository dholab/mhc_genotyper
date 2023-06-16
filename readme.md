# MHC Genotyper Illumina Short reads in select non-human primates.
## Illumina Fastq file naming nomenclature
- MHC10001b_s80_L001_R1_001.fastq.gz
- MHC10001b_s80_L001_R2_001.fastq.gz
    - \<Sample_ID\>_s\<SampleNumber\>_L\<LaneNumber\>_R\<PairDirection\>_\<PoolNumber\>.fastq.gz
## Sample Sheet format:
Sample_ID,Sample_Name,Description,I7_Index_ID,index,Sample_Project,Species
## Required and Pivot Table Columns:
- Sample_ID : (required string) Case-Sensitive. This is what illumina puts in the file name.  And will be the "gs_id" in the pivot table
- Sample_Name : (required string) Case-Sensitive. This becomes the "Client_id" in the pivot table
- Species: (required string) Not case sensitive. The following are acceptable: Rhesus, Mamu, Pig-tailed, Mane, Cynomolgus, MCM.
- Comments: (optional, string) A comment that is added to the pivot table if it exists.
- Sample_Project: (require string) Case-Sensitive. Creates a separate pivotable for each Species/Sample_Project combination

# Prerequisites:
## Precompiled Pre Requisites:
- Debian, Ubuntu, or MacOSx required (not windows compatible)
  - Tested on MacOsx + Python 3.6.5 and  Google Colab).
  - Python 3.7 and later (SHOULD) also work.
- download BBMAP_39.01.tar.gz and untar it (i.e. tar -xzvf BBMAP_39.01.tar.gz)
  - https://sourceforge.net/projects/bbmap/files/BBMap_39.01.tar.gz/download
- download vsearch with your appropriate OS, untar and create a symbolic link as needed
  - https://github.com/torognes/vsearch/releases
    - (i.e. tar -xzvf vsearch-2.21.1-macos-x86_64.tar.gz)
```
ln -s ~/vsearch-2.21.1-macos-x86_64/bin/vsearch /usr/local/bin/vsearch
ln -s ~/bbmap/bbduk.sh /usr/local/bin/bbduk.sh
ln -s ~/bbmap/bbmerge.sh /usr/local/bin/bbmerge.sh
ln -s ~/bbmap/bbmap.sh /usr/local/bin/bbmap.sh
ln -s ~/bbmap/stats.sh /usr/local/bin/stats.sh
```
## Install more python libraries and apt-get prereqs.
```bash
pip install xlsxwriter
pip install pandas
pip install numpy
apt install -y pigz
```
# Run The pipeline
```bash
python3 ~/github/mhc_genotyper/main.py \
--run_id 1234 \
--input_dir ~/BaseSpace/fastqReads \
--genotyper_root_dir ~/github/mhc_genotyper
```