"""Contains the JunctionEntry and JunctionData objects"""
from typing import List, Tuple

class JunctionEntry:
    """Represents a single junction entry from an RNA construct""" 
    def __init__(self, **kwargs ):
		#sequence, structure, reactivity, construct, sn, reads, score 
        # TODO add some checks here... I think it will make 
        # it better if the checks occur before the actual variable
        # assignment? not sure about that one though 
        assert len( structure ) == len( sequence ) and len( sequence ) == len( reactivity )
        assert len( sequence ) 
        self.__sequence = kwargs.get(sequence, None)
        self.__structure = kwargs.get(structure, None)
        self.__reactivity = kwargs.get(reactivity, None)
        self.__construct = kwargs.get(construct, None)
        self.__sn = kwargs.get(sn, None)
        self.__reads = kwargs.get(reads, None)
        self.__score = kwargs.get(score, None)
        self.__symmetrical = is_symmetrical( self.sequence )

    def key( self ) -> Tuple[str,str]:
        return ( self.__sequence, self.__structure )

    def is_symmetrical( self ) -> bool:
        return sel.symmetrical

    def __getitem__( self, idx : int):
        return self.reactivity[ idx ]


class JunctionData:

    def __init__( self, sequence, structure, entries ):
        self.sequence = sequence
        self.structure = structure
        self.entries = entries
        self.symmetrical = is_symmetrical( sequence )
        self.data = [ [] for _ in self.sequence ]
        
        for entry in self.entries:
            for idx, nt in enumerate( self.sequence ):
                if nt == 'A' or nt == 'C':
                    self.data[ idx ].append( entry[idx] )

        self.data = [ np.array(vals) for vals in self.data ]

    def rebuild_data( self ):
        self.data = [ [] for _ in self.sequence ]
        
        for entry in self.entries:
            for idx, nt in enumerate( self.sequence ):
                if nt == 'A' or nt == 'C':
                    self.data[ idx ].append( entry[idx] )

        self.data = [ np.array(vals) for vals in self.data ]


    def is_symmetrical( self ):
        return self.symmetrical


    def legal_pairs( self ):
        LEGAL = []
        for idx, nt in enumerate( self.sequence ):
            # TODO should I change this so that we only consider
            if nt == 'A' or nt == 'C':
                LEGAL.append( idx )
        pairs = set()
        for i1 in LEGAL:
            for i2 in LEGAL:
                if i1 != i2:
                    pairs.add( (min(i1, i2),max(i1,i2) ) )
        return list( pairs ) 



    def significant_pairs( self, cutoff=0.05 ):
        
        result = []
        for (i1, i2) in self.legal_pairs():
            
            _, p_val = wilcoxon( self.data[i1], self.data[i2], alternative='greater' )
            if p_val <= cutoff:
                result.append((i1, i2))
                continue
            
            _, p_val = wilcoxon( self.data[i2], self.data[i1], alternative='greater' )
            if p_val <= cutoff:
                result.append((i2, i1))
                continue
       
        return result


    def plot( self, plot_dir ):
        fname = f"{plot_dir}/{self.structure}_{self.sequence}.png"
        if os.path.isfile( fname ):
            # if it's already been made, then skip it
            return

        fig, axes = plt.subplots(1, 1, figsize=(15,10))
        data = { 'pos': [], 'react': [] }
        outliers = 0
        for idx, vals in enumerate( self.data ):
            if not len(vals):
                data['pos'].append( idx )
                data['react'].append( -5 )

            outliers += count_outliers( vals )
            for v in vals:
                data['pos'].append( idx )
                data['react'].append( v )
        cmap = dict()
        for ii, nt in enumerate( self.sequence ):
            if nt == 'A':
                cmap[ii] = 'red'
            elif nt == 'T' or nt == 'U':
                cmap[ii] = 'green'
            elif nt == 'G':
                cmap[ii] = 'orange'
            elif nt == 'C':
                cmap[ii] = 'blue'
            else:
                cmap[ii] = 'black'

        df = pd.DataFrame( data )
        x = list(range(len(self.sequence))) 
        ymax = df['react'].max()
        sns.boxplot(x='pos', y='react', data=data, palette=cmap )
        sns.stripplot(x='pos', y='react', data=data, color='black')
        axes.set_ylabel('reactivity (a.u.)')
        axes.set_xticks( x )
        axes.set_xticklabels([ f"{s}\n{nt}" for s,nt in zip(self.sequence, self.structure ) ] )
        name = f"{self.structure}_{self.sequence} (N = {len(self.entries)}), outliers = {outliers}"
        #axes.set_title( name ) 
        fig.suptitle( name )
        axes.set_ylim((0, ymax*1.25)) 
        plt.savefig( fname )
        plt.clf()
        plt.close( fig )
    
    def show( self ):
        #fname = f"{plot_dir}/{self.structure}_{self.sequence}_{int(temperature)}.png"
        #if os.path.isfile( fname ):
        #    # if it's already been made, then skip it
        #    return

        fig, axes = plt.subplots(1, 1, figsize=(9,6))
        data = { 'pos': [], 'react': [] }
        for idx, vals in enumerate( self.data ):
            if not len(vals):
                data['pos'].append( idx )
                data['react'].append( -5 )

            for v in vals:
                data['pos'].append( idx )
                data['react'].append( v )
        cmap = dict()
        for ii, nt in enumerate( self.sequence ):
            if nt == 'A':
                cmap[ii] = 'red'
            elif nt == 'T' or nt == 'U':
                cmap[ii] = 'green'
            elif nt == 'G':
                cmap[ii] = 'orange'
            elif nt == 'C':
                cmap[ii] = 'blue'
            else:
                cmap[ii] = 'black'

        df = pd.DataFrame( data )
        x = list(range(len(self.sequence))) 
        ymax = df['react'].max()
        sns.boxplot(x='pos', y='react', data=data, palette=cmap )
        sns.stripplot(x='pos', y='react', data=data, color='black')
        axes.set_ylabel('reactivity (a.u.)')
        axes.set_xticks( x )
        axes.set_xticklabels([ f"{s}\n{nt}" for s,nt in zip(self.sequence, self.structure ) ] )
        #name = f"{self.structure}_{self.sequence} (N = {len(self.entries)}), T={temperature}C"
        #axes.set_title( name ) 
        fig.suptitle( f"N={len(self.entries)}" )
        axes.set_ylim((0, ymax*1.25)) 
        plt.show()
        #plt.savefig( fname )
        #plt.clf()
        #plt.close( fig )
    
    def bind( self, ax ):
        #fname = f"{plot_dir}/{self.structure}_{self.sequence}_{int(temperature)}.png"
        #if os.path.isfile( fname ):
        #    # if it's already been made, then skip it
        #    return

        #fig, axes = plt.subplots(1, 1, figsize=(9,6))
        data = { 'pos': [], 'react': [] }
        for idx, vals in enumerate( self.data ):
            if not len(vals):
                data['pos'].append( idx )
                data['react'].append( -5 )

            for v in vals:
                data['pos'].append( idx )
                data['react'].append( v )
        cmap = dict()
        for ii, nt in enumerate( self.sequence ):
            if nt == 'A':
                cmap[ii] = 'red'
            elif nt == 'T' or nt == 'U':
                cmap[ii] = 'green'
            elif nt == 'G':
                cmap[ii] = 'orange'
            elif nt == 'C':
                cmap[ii] = 'blue'
            else:
                cmap[ii] = 'black'

        df = pd.DataFrame( data )
        x = list(range(len(self.sequence))) 
        ymax = df['react'].max()
        sns.boxplot(x='pos', y='react', data=data, palette=cmap, ax=ax )
        sns.stripplot(x='pos', y='react', data=data, color='black', ax=ax)
        ax.set_ylabel('reactivity (a.u.)')
        ax.set_xticks( x )
        ax.set_xticklabels([ f"{s}\n{nt}" for s,nt in zip(self.sequence, self.structure ) ] )
        #name = f"{self.structure}_{self.sequence} (N = {len(self.entries)}), T={temperature}C"
        #axes.set_title( name ) 
        #fig.suptitle( f"N={len(self.entries)}" )
        ax.set_ylim((0, ymax*1.25)) 
        #plt.show()

