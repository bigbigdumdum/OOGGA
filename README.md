# Overhang Optimizer for golden gate assembly (OOGGA)

Python package for fragmenting DNA sequences optimally for Golden Gate assembly cloning.

How the OOGGA works and OOGGA citation:

> Overhang Optimizer for Golden Gate Assembly (OOGGA) S. Mukundan, M. S. Madhusudhan bioRxiv 2025.06.16.659877; doi: https://doi.org/10.1101/2025.06.16.659877

---

## Installation

```bash
pip install oogga
```

---

## Getting the overhang scores

OOGGA requires a CSV scoring table from Potapov et al. 2018. Please use the data from the journal website and cite them. 

> Potapov V, Ong JL, Kucera RB, Langhorst BW, Bilotti K, Pryor JM, Cantor EJ, Canton B, Knight TF, Evans TC Jr, Lohman GJS. Comprehensive Profiling of Four Base Overhang Ligation Fidelity by T4 DNA Ligase and Application to DNA Assembly. ACS Synth Biol. 2018 Nov 16;7(11):2665-2674. doi: 10.1021/acssynbio.8b00333. Epub 2018 Oct 29. PMID: 30335370.

1. Download the zip from:  
   https://pubs.acs.org/doi/suppl/10.1021/acssynbio.8b00333/suppl_file/sb8b00333_si_002.zip

2. Extract it and open **`FileS04_T4_18h_37C.xlsx`** in Microsoft Excel or LibreOffice Calc.

3. Save as **`FileS04_T4_18h_37C.csv`** (in CSV format).

4. Pass the path to this CSV file with the `-data_file` argument or copy to ./lib/.

---

## Usage

### See all command line options

```bash
oogga -h 
```
A normal example of splitting a plasmid into fragments of length range 500-1000
```bash
oogga inputplasmid.fasta 500 1500 output
```
---

### Evaluate overhangs (`eval-frags`)
This script is used in the manuscript
```bash
eval-frags AACC TTGG ACGT TGCA
```


