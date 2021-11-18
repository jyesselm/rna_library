# RNA Library 
[![Documentation Status](https://readthedocs.org/projects/rna-library/badge/?version=latest)](https://rna-library.readthedocs.io/en/latest/?badge=latest)

## Installation
1. Dependencies
	`$ pip install -r requirements.txt`
1. package installation
	`$ python setup.py install`

## Basic Usage

### Parsing Structures
The folllwing code:
```
import rna_library as rl

for m in rl.SecStruct('(((...(((...(((...)))...)))...)))', 
			'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'):
	print(m)

```
Should produce this:
```
Helix,NNN&NNN,(((&)))
Junction,NNNNN&NNNNN,(...(&)...)
Helix,NNN&NNN,(((&)))
Junction,NNNNN&NNNNN,(...(&)...)
Helix,NNN&NNN,(((&)))
Hairpin,NNNNN,(...)
```
### Analyzing Mutational Histograms
```
import rna_library as rl



# this step creates the population averages from the mutational histogram 
# file and saves the .csv files into './outdir/'
rl.process_histos('/path/to/mutational-histo.p', 'outdir')
reactivity_df = rl.build_react_df(
	out_dir='outdir',
	start_seq='GGGCUUCGGCCC',
	end_seq="AAAGAAACAACAACAACAAC",
	fasta_file='/path/to/fasta.fasta',
	histos_file='/path/to/mutational-histo.p'
	)

# normalization via hairpin
reactivity_df['reactivity'] = rl.processing.normalize_hairpin(reactivity_df, 'CGAGUAG','(.....)')
# creating the motif dataframe
motif_df = rl.build_motif_df( reactivity_df )
```
