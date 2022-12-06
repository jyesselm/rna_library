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

