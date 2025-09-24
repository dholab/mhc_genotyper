# MHC Genotyper

MHC genotyping pipeline for Illumina short-read sequencing data from non-human
primates.

## Overview

This pipeline performs MHC genotyping on Illumina paired-end sequencing data
from three non-human primate species:

- Rhesus macaques (Mamu)
- Cynomolgus macaques (MCM/Mafa)
- Pig-tailed macaques (Mane)

The pipeline identifies MHC alleles and calls haplotypes based on diagnostic
allele combinations across seven MHC loci (A, B, DRB, DQA, DQB, DPA, DPB).

## Installation

### Prerequisites

- Linux or macOS (not Windows compatible)
- Python 3.9 or later

### Setup with Pixi

This project uses Pixi to manage both conda and Python dependencies:

```bash
# Install pixi if not already installed
curl -fsSL https://pixi.sh/install.sh | bash

# Clone the repository
git clone https://github.com/nicholasrminor/mhc_genotyper.git
cd mhc_genotyper

# Install all dependencies
pixi install
```

This will automatically install:

- Bioinformatics tools: BBMap, VSEARCH, pigz
- Python packages: Jupyter, pandas, numpy, xlsxwriter, pysam

## Usage

This pipeline is designed as an interactive Jupyter notebook following literate
programming principles. Users work through the notebook cell by cell, with
documentation and code interspersed.

### Starting the Pipeline

```bash
# Launch the notebook
pixi run notebook
```

This will open the `mhc_genotyper.ipynb` notebook in your browser.

### Working Through the Notebook

1. **Configuration**: The first sections of the notebook contain configuration
   cells where you specify:
   - Experiment number
   - Input folder path containing FASTQ files
   - Sample sheet location
   - Project list (optional)

2. **Execution**: Follow the notebook cells sequentially. Each section is
   documented with its purpose and expected outputs.

3. **Reference Files**: The pipeline uses reference sequences located in the
   `ref/` directory:
   - `26128_ipd-mhc-mamu-2021-07-09.miseq.RWv4.fasta` - Rhesus macaque MHC
     alleles
   - `MCM_MHC-all_mRNA-MiSeq_singles-RENAME_20Jun16.fasta` - Cynomolgus macaque
     MHC alleles
   - `Mafa_MiSeq-IPD_17.06.01_2.2.0.0_plus_SP.fasta` - Alternative Cynomolgus
     reference
   - `SBT195_MHCII_primers_2Sep13.fasta` - PCR primer sequences

### Important Note for Local Execution

The notebook is currently configured for Google Colab with Google Drive paths.
For local execution, you need to modify the path configuration in the notebook:

1. Locate the cell containing
   `genotyper_root_dir = os.path.join(Shareddrives_path,'dholab/gs/genotyper')`
2. Replace it with your local path, for example:
   ```python
   genotyper_root_dir = '/path/to/your/mhc_genotyper'
   ```
3. Update the reference dictionary paths to match the actual file locations in
   the `ref/` directory:
   ```python
   ref_dict = {
       'REF': {
           'MCM': os.path.join(genotyper_root_dir, 'ref/MCM_MHC-all_mRNA-MiSeq_singles-RENAME_20Jun16.fasta'),
           'MANE': os.path.join(genotyper_root_dir, 'ref/Mafa_MiSeq-IPD_17.06.01_2.2.0.0_plus_SP.fasta'),
           'MAMU': os.path.join(genotyper_root_dir, 'ref/26128_ipd-mhc-mamu-2021-07-09.miseq.RWv4.fasta')
       },
       'PCR_PRIMERS': os.path.join(genotyper_root_dir, 'ref/SBT195_MHCII_primers_2Sep13.fasta')
   }
   ```

## Input Requirements

### FASTQ Files

Paired-end Illumina sequencing files following the naming convention:

- `<Sample_ID>_S<Number>_L<Lane>_R1_001.fastq.gz`
- `<Sample_ID>_S<Number>_L<Lane>_R2_001.fastq.gz`

### Sample Sheet

CSV file with the following required columns:

- `Sample_ID`: Illumina sample identifier (matches FASTQ filename)
- `Sample_Name`: Internal sample name
- `Sample_Project`: Project identifier
- `Description`: Sample description (can be blank)
- `Species`: One of: Rhesus, Cynomolgus, Pig-tailed (or Mamu, MCM, Mane)

Example:

```csv
Sample_ID,Sample_Name,Sample_Project,Description,Species
74564,DW427,Kenyon19,H20C31,Cynomolgus
```

## Output

The pipeline generates:

- Excel pivot tables with genotyping results per project/species combination
- Read counts for each identified allele
- Haplotype calls for each MHC locus
- Quality metrics (mapped/unmapped read percentages)
- Intermediate SAM alignment files

Output files are organized by project and timestamped for reproducibility.

## Pipeline Steps

1. Primer trimming using BBDuk
2. Read merging using BBMerge
3. Unique sequence identification using VSEARCH
4. Mapping to reference alleles using BBMap
5. Chimera and artifact removal using VSEARCH UNOISE
6. Haplotype calling based on diagnostic allele patterns
7. Excel report generation with pivot tables

## Checking Dependencies

Verify that all required tools are installed:

```bash
pixi run check-deps
```

## License

See LICENSE file for details.

## Citation

If you use this pipeline in your research, please cite: [Citation information to
be added]
