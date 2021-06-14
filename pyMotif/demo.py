import motif

root = motif.parse_to_motifs('(((...)))..(((.(((...))).)))...',
                        'GGGAAACCCAAGGGAGGGAAACCCACCCAAA'
                        )


stats = {"singlestrand": 0, "hairpin": 0, "junction":0, "helix": 0}

def traverse( m : motif.Motif, stats ):
    """Function that returns number of motifs in the structure"""
    print(m.type(),  m.strands()) 
    
    for c in m.children():
        traverse( c, stats )

traverse( root, stats )


