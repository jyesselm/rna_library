# RNA Library 
[![Documentation Status](https://readthedocs.org/projects/rna-library/badge/?version=latest)](https://rna-library.readthedocs.io/en/latest/?badge=latest)

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
```
