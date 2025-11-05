# Meta-analysis of Negative Crossover Interference

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

## Overview

This repository contains data, code, and supplementary materials for the doctoral dissertation:

**"Meta-analysis of negative crossover interference: variation between species, environments, and sexes"**

**Author:** Shaul Sapielkin  
**Supervisors:** Prof. Abraham B. Korol, Dr. Eyal Privman  
**Institution:** University of Haifa, Department of Evolutionary and Environmental Biology  
**Year:** 2025

---

## Abstract

Crossover interference (COI) is a fundamental mechanism regulating meiotic recombination distribution along chromosomes. While positive COI has been extensively studied, negative crossover interference (NCOI)—where crossovers cluster closer than expected by chance—remains poorly understood. This dissertation presents the first comprehensive meta-analysis of NCOI across eukaryotic species, revealing its widespread occurrence, sex-specific patterns, and molecular determinants.

Through analysis of 27 species (2,515 chromosomal intervals) spanning plants, arthropods, and fish, we demonstrate that NCOI occurs in 47.7% of intervals with 96.7% statistical significance. Sex-specific analysis reveals dramatic dimorphism: 65.3% NCOI prevalence in males versus 32.9% in females. Phylogenetic analysis shows extraordinary variation from 14.3% NCOI in maize to 77.8% in wheat.

Methodologically, this work introduces an extended gamma-sprinkled model for simultaneous estimation of positive, negative, and neutral interference, and develops a sex-aware maximum likelihood estimation framework that corrects systematic bias in F2 populations (bias magnitude: 0.21-0.46 CoC units). Molecular analysis of wheat epigenome reveals strong correlations between NCOI and chromatin accessibility (r=0.995, p<0.001), histone modification H3K27ac (r=0.883, p<0.05), and DNA methylation (r=-0.68, p<0.01).

This research establishes NCOI as a biologically significant phenomenon with implications for understanding meiotic regulation, genome evolution, and breeding programs.

---

## Repository Structure

```
NCOI-metaanalysis/
├── Chapter3/ 
├── Chapter5/ 
├── data/                          # All datasets used in the dissertation
├── figures/                       # Figures organized by dissertation chapters  
├── scripts/                       # Python analysis code
├── tables/                        # Supplementary tables and results
└── README.md                      # This file
```

### 📂 `data/`

Contains raw and processed datasets:
- **Segregation data** from 27 species (32 independent studies)
- **Species metadata** (taxonomy, sample sizes, population types)
- **Simulation results** for sex-aware MLE validation
- **Wheat epigenomic features** (chromatin accessibility, histone modifications, DNA methylation, transposons, lncRNAs)

**Key statistics:**
- 27 species analyzed (59% plants, 25% arthropods, 19% fish)
- 2,515 chromosomal intervals
- 1,199 intervals with NCOI (47.7%)
- 1,264 intervals with positive interference (50.3%)

### 📊 `figures/`

High-resolution figures from the dissertation, organized by chapter:
- Publication-quality versions (300+ dpi)
- Vector formats (PDF/SVG) where available

### 💻 `scripts/`

Python analysis code organized by dissertation chapter:

**Data processing and quality control:**
- Quality control using MultiPoint software criteria
- Coefficient of Coincidence (CoC) calculations
- Interference type classification

**Extended gamma-sprinkled model (Chapter 2):**
- Three-component model implementation (ν>1, ν=1, ν<1)
- Maximum likelihood estimation
- Likelihood ratio testing (H₀, H₁, H₂)

**Sex-aware MLE framework (Chapter 3):**
- Sex-aware maximum likelihood estimation for F2 populations
- Simulation study (N=1,000 to 100,000 individuals)
- Bias quantification and empirical validation

**Phylogenetic analysis (Chapter 4):**
- NCOI patterns across phylogenetic groups
- CoC summary statistics by species
- Visualization and statistical comparisons

**Genomic correlations (Chapter 5):**
- Wheat epigenomic feature analysis
- Correlation analysis with NCOI intensity
- Sliding window calculations

### 📋 `tables/`

Supplementary tables including:
- Complete species and dataset information
- Statistical test results
- Correlation matrices
- Simulation validation data

---

## Key Findings

### 1. Widespread NCOI across Eukaryotes

- **NCOI detected in 47.7%** of chromosomal intervals (1,199/2,515)
- **96.7% statistical significance** (p<0.05) among NCOI intervals
- Present across all major phylogenetic groups: plants, arthropods, fish

### 2. Dramatic Sexual Dimorphism

- **65.3% NCOI prevalence in males** (764/1,169 male-informative intervals)
- **32.9% NCOI prevalence in females** (225/683 female-informative intervals)
- **2.0-fold sex bias** in NCOI occurrence

### 3. Phylogenetic Variation

- **Triticum** (wheat): 77.8% NCOI, mean CoC = 3.79
- **Zea mays** (maize): 14.3% NCOI, mean CoC = 1.89
- **Pericentromeric regions** show 2.3-fold elevated NCOI

### 4. Molecular Determinants

Strong correlations with wheat epigenomic features:
- **Chromatin accessibility:** r=0.995, p<0.001
- **H3K27ac (active chromatin):** r=0.883, p<0.05  
- **DNA methylation:** r=-0.68, p<0.01 (negative correlation)
- **Transposable elements:** Positive association with NCOI intensity

### 5. Methodological Advances

- **Sex-aware MLE framework** eliminates systematic bias in F2 populations
- **Bias magnitude:** +0.21 to +0.36 CoC (concordant sex differences), -0.22 to -0.46 CoC (discordant)
- **Extended gamma-sprinkled model** enables simultaneous estimation of all interference types

---

## Installation

### Requirements

- Python 3.8 or higher
- See `requirements.txt` for complete list of dependencies

### Setup

```bash
# Clone the repository
git clone https://github.com/[username]/NCOI-metaanalysis.git
cd NCOI-metaanalysis

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### Using Conda

```bash
# Create conda environment
conda env create -f environment.yml

# Activate environment
conda activate ncoi-metaanalysis
```

---

## Usage

### Reproducing Main Analyses

**Data processing and QC:**
```python
# Example: Quality control and CoC calculation
python scripts/01_data_processing/quality_control.py \
    --input data/raw_segregation_data.csv \
    --output data/processed_intervals.csv \
    --min_individuals 100
```

**Extended gamma-sprinkled model:**
```python
# Estimate interference parameters
python scripts/02_gamma_model/extended_gamma_sprinkled.py \
    --data data/processed_intervals.csv \
    --output results/model_estimates.csv
```

**Sex-aware MLE framework (Chapter 3):**
```python
# Run sex-aware analysis
python scripts/03_sex_aware_MLE/sex_aware_framework.py \
    --input data/F2_datasets.csv \
    --output results/sex_aware_estimates.csv
```

**Phylogenetic analysis (Chapter 4):**
```python
# Analyze phylogenetic patterns
python scripts/04_phylogenetic_analysis/chapter4_main_analysis.py \
    --data data/processed_intervals.csv \
    --output results/phylogenetic_summary.csv
```

**Genomic correlations (Chapter 5):**
```python
# Wheat genomic feature correlations
python scripts/05_genomic_correlations/wheat_feature_analysis.py \
    --features data/wheat_genomic_features/ \
    --ncoi data/wheat_ncoi_data.csv \
    --output results/correlation_analysis.csv
```

---

## Publications

**Chapter 2:** Extended gamma-sprinkled model for crossover interference

### Published / In Press
**Sapielkin, S., Korol, A.B., Privman, E.** (2025). Accounting for sex differences in recombination rates reduces bias in crossover interference estimation in F2 populations. *BMC Genomics* (in press).

### Chapters as Preprints

- **Chapter 3:** Sex-aware maximum likelihood estimation framework
- **Chapter 4:** Phylogenetic patterns of NCOI across eukaryotes
- **Chapter 5:** Genomic determinants of NCOI in wheat

---

## Citation

If you use this data or code, please cite:

### Dissertation

```bibtex
@phdthesis{sapielkin2025ncoi,
  author = {Sapielkin, Shaul},
  title = {Meta-analysis of negative crossover interference: variation between species, environments, and sexes},
  school = {University of Haifa},
  year = {2025},
  doi = {10.5281/zenodo.XXXXXX}
}
```

### Published Article

```bibtex
@article{sapielkin2025sexaware,
  author = {Sapielkin, Shaul and Korol, Abraham B. and Privman, Eyal},
  title = {Accounting for sex differences in recombination rates reduces bias in crossover interference estimation in F2 populations},
  journal = {BMC Genomics},
  year = {2025},
  note = {in press}
}
```

---

## Data Availability

All data used in this dissertation are publicly available:

- **Raw segregation data:** Available in `data/` directory and permanently archived on Zenodo (DOI: 10.5281/zenodo.XXXXXX)
- **Wheat genomic features:** Sourced from IWGSC RefSeq v2.1 and URGI wheat genome browser
- **Simulation datasets:** Complete simulation data available in `data/simulation_data/`

### Data Sources

Original segregation data were compiled from 32 published studies spanning 27 species. Complete citations and data sources are provided in `tables/data_sources.csv` and dissertation Chapter 2 (Materials and Methods).

---

## Software Dependencies

### Core Python Packages

- **numpy** (≥1.21.0) - Numerical computations
- **pandas** (≥1.3.0) - Data manipulation
- **scipy** (≥1.7.0) - Statistical analysis
- **statsmodels** (≥0.13.0) - Statistical modeling
- **scikit-learn** (≥1.0.0) - Machine learning utilities

### Visualization

- **matplotlib** (≥3.4.0) - Plotting
- **seaborn** (≥0.11.0) - Statistical visualization

### Bioinformatics

- **biopython** (≥1.79) - Sequence analysis

See `requirements.txt` for complete list and exact versions.

---

## Repository Organization

### Data Files

- **`data/raw_data/`** - Original segregation data from literature sources
- **`data/processed_data/`** - Quality-controlled and standardized datasets
- **`data/simulation_data/`** - Sex-aware MLE simulation results
- **`data/wheat_genomic_features/`** - Epigenomic data for correlation analysis

### Analysis Scripts

- **`scripts/01_data_processing/`** - Data cleaning, QC, and CoC calculation
- **`scripts/02_gamma_model/`** - Extended gamma-sprinkled model implementation
- **`scripts/03_sex_aware_MLE/`** - Sex-aware MLE framework (Chapter 3)
- **`scripts/04_phylogenetic_analysis/`** - Phylogenetic patterns (Chapter 4)
- **`scripts/05_genomic_correlations/`** - Genomic feature analysis (Chapter 5)
- **`scripts/utils/`** - Utility functions and helpers

### Documentation

Each subdirectory contains a `README.md` with detailed information about:
- File contents and formats
- Usage instructions
- Expected inputs and outputs
- Relevant dissertation chapters

---

## License

**Code:** MIT License - see [LICENSE](LICENSE) file for details

**Data:** Creative Commons Attribution 4.0 International (CC-BY-4.0)

You are free to use, modify, and distribute this code and data with appropriate attribution.

---

## Contact

**Shaul Sapielkin**  
PhD Candidate  
Department of Evolutionary and Environmental Biology  
University of Haifa, Israel

**Email:** shaul@evo.haifa.ac.il  
**ORCID:** https://orcid.org/0000-0001-9375-3850

**Supervisors:**
- Prof. Abraham B. Korol ([email@example.com])
- Dr. Eyal Privman ([email@example.com])

---

## Acknowledgments

I express profound gratitude to my supervisor, Professor Abraham B. Korol, for his unwavering guidance and mentorship. I thank Dr. Eyal Privman for his invaluable insights and collaboration. I am grateful to all researchers who made their segregation data publicly available, enabling this comprehensive meta-analysis. This work was supported by [funding sources].

Special thanks to the University of Haifa and the Department of Evolutionary and Environmental Biology for providing excellent research facilities and environment.

---

## Related Resources

- **University of Haifa:** https://evolution.haifa.ac.il
- **Lab Website:** [link to lab website if available]
- **Extended Methods:** See `tables/supplementary_methods.pdf`
- **Interactive Visualizations:** [link if available]

---

## Version History

- **v1.0.0** (2025-10-28) - Initial release with complete dissertation data and code
- **v1.1.0** (TBD) - Post-defense updates and corrections

---

## Issues and Contributions

Found a bug or have a suggestion? Please open an issue on GitHub or contact the author directly.

Contributions are welcome! Please feel free to submit pull requests with improvements or corrections.

---

**Last updated:** October 28, 2025
