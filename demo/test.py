#!/usr/bin/env python3

import rna_library


#for k in rna_library.__dict__.keys():
#    print( k ) 

ss = rna_library.SecStruct(
    '((((((...)))(((...)))(((...))))))',
    'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')

for idx,m in ss.itermotifs():
    print( idx,m) 



#print( rna_library.build_barcodes( '(((((....)))))') ) 
