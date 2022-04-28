"""Contains the JunctionEntry() and JunctionData() objects used for storing parsed
DMS data. In addition to storing data, these objects provide information on where
their data comes from. A single JunctionEntry() object holds the reactivity data
for a single junction, as well as the reads, construct, and signal to noise for 
that construct. The JunctionData() object is composed of all JunctionEntry()'s 
with the same sequence/structure and provides composite analysis and functions
to visualize the underlying data.


Author: Chris Jurich <cjurich2@huskers.unl.edu>
Date: 2022-04-27
"""
from __future__ import annotations
import os
import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
from typing import Dict, List
from .stats import EXCESS_MAPPER, comparative_variance
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict
from rna_library.core import is_symmetrical, InvalidArgument, safe_mkdir


from collections import namedtuple

QualitySummary = namedtuple('QualitySummary', 'sn reads')
"""Named tuple for holding the signal-to-noise (sn) and reads for each JunctionEntry()'s parent object."""

class JunctionEntry:
    """Represents a single junction entry from an RNA construct. Holds the reactivity data
		and describes the structure of the junction and its originating construct.

		Attributes:
			__sequence : A str() of the structure's sequence including '&'. Ex: "GAAC&GAUAC".
			__structure : A str() dot-bracket of the structure. Ex: "(..(&)...)".
			__construct : The str() of the original design its from.
			__reactivity : An array of the reactivity data. Same len() as __sequence and __structure with a -1 for the '&' index.
			__sn : The signal to noise of the original construct.
			__reads : Number of reads in the original construct as an int().
			__score : DSCI score for the original structure.
			__indices : A list() of indices of where the junction is in the original structure().
			__external : boolean flag indicating if the JunctionEntry() is the least nested motif in its structure.
	"""

    def __init__(self, **kwargs):
        """Initializes the JunctionEntry() object with a number of keyword data members. Performs basic validation."""
        # now do the variable assigment
        self.__sequence = kwargs.get("sequence", None)
        self.__structure = kwargs.get("structure", None)
        self.__reactivity = kwargs.get("reactivity", None)
        self.__construct = kwargs.get("construct", None)
        self.__sn = kwargs.get("sn", None)
        self.__reads = kwargs.get("reads", None)
        self.__score = kwargs.get("score", None)
        self.__symmetrical = is_symmetrical(self.__sequence)
        self.__indices = kwargs.get('indices', None)
        self.__external = kwargs.get('external', None)
        # argument validation
        self.validate_arguments_()

    
    def reads(self) -> int:
        """
		Getter for the number of reads in the original construct.
		:rtype: int
		"""
        return self.__reads

    def sn(self) -> int:
        """
		Getter for the signal to noise ratio in the original construct.
		:rtype: float
		"""
        return self.__sn

    def is_external(self) -> bool:
        """
		Getter for whether the junction is the least nested motif in its structure.

		:rtype: bool
		"""
        return self.__external

    def indices(self) -> List[int]:
        """Gets the indices for where the junction is found in the original structure.
        :rtype: List[int]
		"""
        return self.__indices

    def flip(self) -> None:
        """Method that flips the reactivity values in the event that the sequence is flipped.
		Also flips the structure and sequence.
		
		>>> je = JunctionEntry(sequence='GG&CAAC', structure='((&)..)', reactivity=[0, 0, -1, 0.25, 1.1, 1.2, 0.3])
		>>> je.key()
		('GG&CAAC', '((&)..)')
		>>> je.reactivity()
		[0, 0, -1, 0.25, 1.1, 1.2, 0.3]
		>>> je.flip()
		>>> je.key()
		('GG&CAAC', '(..(&))')
		>>> je.reactivity()
		[0.25, 1.1, 1.2, 0.3, -1, 0, 0]

		:rtype: None
		"""
        it = self.__sequence.find('&')
        left, right = self.__reactivity[0:it], self.__reactivity[it+1:]
        new_reactivity = right + [0.0] + left
        assert len(new_reactivity) == len(self.__reactivity)
        assert abs(np.sum(new_reactivity) - np.sum(self.__reactivity)) <= 0.05, f"{abs(np.sum(new_reactivity) - np.sum(self.__reactivity))}"
        self.__reactivity = new_reactivity
        seq_left, seq_right = self.__sequence.split('&')
        self.__sequence = seq_right + '&' + seq_left
        ss_left, ss_right = self.__structure.split('&')
        self.__structure = ss_right + '&' + ss_left

    def validate_arguments_(self) -> None:
        """
        Helper method that validates arguments in the constructor. Creates a detailed error message for the user if there are 
		missing data members.
        """
        # some basic size checks
        assert len(self.__structure) == len(self.__sequence)
        assert len(self.__sequence) == len(self.__reactivity)
        assert len(self.__sequence)
        if self.__indices:
            assert len(self.__sequence) == len(self.__indices)
        invalid_args, error_msg = 0, ""
        if not self.__sequence:
            invalid_args += 1
            error_msg += f"no value supplied for sequence\n"
        if not self.__structure:
            invalid_args += 1
            error_msg += f"no value supplied for structure\n"
        if not self.__reactivity:
            invalid_args += 1
            error_msg += f"no value supplied for reactivity\n"
        if not self.__construct:
            invalid_args += 1
            error_msg += f"no value supplied for construct\n"
        if not self.__sn:
            invalid_args += 1
            error_msg += f"no value supplied for sn\n"
        if not self.__reads:
            invalid_args += 1
            error_msg += f"no value supplied for reads\n"
        if not self.__score:
            invalid_args += 1
            error_msg += f"no value supplied for score\n"

        if invalid_args:
            raise InvalidArgument(error_msg)


    def key(self) -> Tuple[str, str]:
        """
        Getter that accesses the (sequence, structure) key for the JunctionEntry

        :rtype: Tuple[str,str]
        """
        return (self.__sequence, self.__structure)

    def is_symmetrical(self) -> bool:
        """
        Getter that checks if the JunctionEntry is for a symmetrical unction.
        """
        return sel.__symmetrical

    def reactivity( self ) -> List[float]: 
        """
        Getter for the reactivity values.
        """
        return self.__reactivity
    
    def sequence( self ) -> str: 
        """
        Getter for the sequence value.
        """
        return self.__sequence
    
    def construct( self ) -> str: 
        """
        Getter for the construct value.
        """
        return self.__construct


    def get_dms_active( self ) -> np.array:
        """Collapses and returns only the DMS active data. Loses significant information about the data and indices.
		>>> je = JunctionEntry(sequence='GAC&GAC', structure='(.(&).)',reactivity=[0, 0.5, 0.2, -1, 0, 0.6, 0.15])
		>>> je.get_dms_active()
		[0.5, 0.2 0.6, 0.15]

		:rtype: np.array[float]
		"""
        result = []
        assert len(self.__reactivity) == len(self.__sequence)
        for react, nt in zip( self.__reactivity, self.__sequence ):
            if nt == 'A' or nt == 'C':
                result.append( react )

        return np.array( result )

    def __getitem__(self, idx: int) -> float:
        """Directly indexes into the JunctionEntry.__reactivity array for retrieval ONLY.
		:rtype: float
		""" 
        return self.__reactivity[idx]

class JunctionData:
    """
    Composite class that represents a collection of JunctionEntry() objects in an experiment that have
	the same sequence. Provides functionality to visualize and describe the composite of that junction in 
	an experiment.

	Attributes:
		sequence : The dot-bracket str() of the junction including '&'. Ex: 'CAC&GAG'.
		structure : The secondary structure str() of the junction including '&'. Ex: '(.(&).)'.
		entries : A list() of JunctionEntry() objects.
		symmetrical : Whether the junction is symmetrical in terms of junction sequence size.
		data : Composite of all reactivities in the JunctionEntry() objects TODO(FINISH)..
    """

    def __init__(self, **kwargs):
        """Initializes the JunctionData() object with a sequence, structure and list() of JunctionEntry() objects."""
        # sequence, structure, entries
        self.sequence = kwargs.get("sequence", None)
        self.structure = kwargs.get("structure", None)
        self.entries = kwargs.get("entries", None)
        self.symmetrical = is_symmetrical(self.sequence)
        self.data = [[] for _ in self.sequence]
        self.error = [[] for _ in self.sequence]

        invalid_args = 0  # TODO create a validation function for this
        error_msg = ""
        if not self.sequence:
            invalid_args += 1
            error_msg += f"no value supplied for sequence\n"
        if not self.structure:
            invalid_args += 1
            error_msg += f"no value supplied for structure\n"
        if not self.entries:
            invalid_args += 1
            error_msg += f"no value supplied for entries\n"

        if invalid_args:
            raise InvalidArgument(error_msg)

        self.rebuild_data()

    def quality_summary(self) -> Dict[str,QualitySummary]:
        """Provides a quality summary of the JunctionData()'s JunctionEntry()'s. Returns a dict()
		with (key,value) pairs of (construct name, QualitySummary()) where QualitySummary() is
		a named tuple holding both sn and read values for the original constructs.

		:rtype: Dict[str, QualitySummary]
		"""
        result = dict()
        ee : JunctionEntry
        for ee in self.entries:
            result[ee.construct()] = QualitySummary(sn=ee.sn(), reads=ee.reads())
        return result
            

    def flip( self ) -> None:
        """Flips the two strands in the data and updates the sequence, structure, data and child
		JunctionEntry() objects.
		
		see JunctionEntry.flip() for additional information.

		:rtype: None
		"""
        for idx in range(len(self.entries)):
            self.entries[idx].flip()
        
        seq_left, seq_right = self.sequence.split('&')	 
        ss_left, ss_right = self.structure.split('&')
        self.sequence = seq_right + '&' + seq_left
        self.structure = '(' + (len(ss_right)-2)*'.' + '(&)' + (len(ss_left)-2)*'.' + ')'
        self.rebuild_data()


    def merge(self, other : JunctionData ) -> None:
        """Merges the entries from another JunctionData() object into this one and re-builds the data.
		Checks if sequences and structures are same and throws if not.
		
		>>> #initialized elsewhere
		>>> jd1 : JunctionData
		>>> jd2 : JunctionData
		>>> jd1.num_entries()
		10
		>>> jd2.num_entries()
		23
		>>> jd1.merge( jd2 )
		>>> jd1.num_entries()
		33
		>>> jd2.num_entries()
		23


		:rtype: None
		"""
	
        if self.sequence != other.sequence:
            raise Exception(f"Cannot merge JunctionData() objects with different sequences: '{self.sequence}' vs '{other.sequence}'")
        
        if self.structure != other.structure:
            raise Exception(f"Cannot merge JunctionData() objects with different structures: '{self.structure}' vs '{other.structure}'")
        
        self.entries.extend( other.entries )
        self.rebuild_data()

    def get_active_data(self) -> List[List[float]]: 
        """
		Aggregates on the active data (where nt's are 'A' or 'C') and returns it as a list() of list()'s with 
		dimensions MxN where M is the JunctionData.num_entries() and M is the len(JunctionEntry.get_dms_active()).

		see JunctionEntry.get_dms_active() for more information.

		:rtype: List[List[float]]
		"""
        active = []
        for row, nt in zip(self.data, self.sequence):
            if nt == "A" or nt == "C":
                active.append(row)
        return active

    def rebuild_data(self) -> None:
        """
        Method that rebuilds the internal data representation from the JunctionEntry objects.

        :rtype: NoneType
        """
        self.data = [[] for _ in self.sequence]

        for entry in self.entries:
            for idx, nt in enumerate(self.sequence):
                if nt == "A" or nt == "C":
                    self.data[idx].append(entry[idx])

        self.data = [np.array(vals) for vals in self.data]

    def is_symmetrical(self) -> bool:
        """
        Getter that tells if the current JunctionData object models a symmetrical junction.
        :rtype: bool
        """
        return self.symmetrical

    def measure_variance( self ) -> Dict[str,float]:
        """
	    Measures the pairwise variance between the different value series for each JunctionEntry(). Not used
		for any current protocols.

		:rtype: Dict[str,float]
        """
        result = dict()
        constructs = list(map(lambda e: e.construct(), self.entries))
        reacts = list(map(lambda e: e.get_dms_active(), self.entries))
        for idx1, con1 in enumerate( constructs ):
            for idx2, con2 in enumerate( constructs ):
                if idx1 == idx2:
                    continue
                variance = comparative_variance( reacts[idx1], reacts[idx2] )
                result[con1] = variance
                result[con2] = variance
        return result


    def plot(self, plot_dir: str, overwrite=False) -> None:
        """
        Method that saves a plot of the JunctionData's data points to the supplied directory. 

        :param: str plot_dir: The directory where the plot will be saved. Does not have to exist.
        :rtype: NoneType
        """
        safe_mkdir(plot_dir)
        fname = f"{plot_dir}/{self.structure}_{self.sequence}.png"
        if not overwrite and os.path.isfile(fname):
            return
        name = f"{self.structure}_{self.sequence}, N={len(self.entries)}"
        fig, axes = plt.subplots(1, 1, figsize=(15, 10))
        fig.suptitle(name)
        self.bind(axes)
        plt.savefig(fname)
        plt.clf()
        plt.close(fig)


    def error_plots(self, dirname=None) -> None:
        """
        Method that brings up a plot of the JunctionData's data points
        
        :rtype: NoneType
        """
        fig, ax = plt.subplots(1, 1, figsize=(9, 6))
        fig.suptitle(f"N={len(self.entries)}")
        data = {"pos": [], "react": []}
        for idx, vals in enumerate(self.data):
            if not len(vals):
                data["pos"].append(idx)
                data["react"].append(-5)

            for v in vals:
                data["pos"].append(idx)
                data["react"].append(v)
        cmap = dict()
        for ii, nt in enumerate(self.sequence):
            if nt == "A":
                cmap[ii] = "red"
            elif nt == "T" or nt == "U":
                cmap[ii] = "green"
            elif nt == "G":
                cmap[ii] = "orange"
            elif nt == "C":
                cmap[ii] = "blue"
            else:
                cmap[ii] = "black"

        df = pd.DataFrame(data)
        x = list(range(len(self.sequence)))
        ymax = df["react"].max()
        #sns.boxplot(x="pos", y="react", data=data, palette=cmap, ax=ax)
        sns.stripplot(x="pos", y="react", data=data, color="black", ax=ax)
        avgs = np.nan_to_num( np.array(list(map(np.mean, self.data))),-5)
        ax.errorbar(list(range(len(self.data))),avgs, yerr=self.error)
        ax.set_ylabel("reactivity (a.u.)")
        ax.set_xticks(x)
        ax.set_xticklabels(
            [f"{s}\n{nt}" for s, nt in zip(self.sequence, self.structure)]
        )
        ax.set_ylim((0, ymax * 1.25))
        fname = f"{dirname}/{self.structure}_{self.sequence}.png"
        if not dirname:
            plt.show()
        else:
            plt.savefig(fname)
            plt.clf()
            plt.close(fig)


    def show(self) -> None:
        """
        Method that brings up a plot of the JunctionData's data points
        
        :rtype: NoneType
        """

        fig, axes = plt.subplots(1, 1, figsize=(9, 6))
        fig.suptitle(f"N={len(self.entries)}")
        self.bind(axes)
        plt.show()

    def bind(self, ax: matplotlib.axes.Axes) -> None:
        """
        Method that binds the JunctionData points to a supplied matplotlib Axes object.

        :param: matplotlib.axes.Axes ax: the Axes object which the plot will be bound to
        :rtype: NoneType
        """
        data = {"pos": [], "react": []}
        for idx, vals in enumerate(self.data):
            if not len(vals):
                data["pos"].append(idx)
                data["react"].append(-5)

            for v in vals:
                data["pos"].append(idx)
                data["react"].append(v)
        cmap = dict()
        for ii, nt in enumerate(self.sequence):
            if nt == "A":
                cmap[ii] = "red"
            elif nt == "T" or nt == "U":
                cmap[ii] = "green"
            elif nt == "G":
                cmap[ii] = "orange"
            elif nt == "C":
                cmap[ii] = "blue"
            else:
                cmap[ii] = "black"

        df = pd.DataFrame(data)
        x = list(range(len(self.sequence)))
        ymax = df["react"].max()
        sns.boxplot(x="pos", y="react", data=data, palette=cmap, ax=ax)
        sns.stripplot(x="pos", y="react", data=data, color="black", ax=ax)
        ax.set_ylabel("reactivity (a.u.)")
        ax.set_xticks(x)
        ax.set_xticklabels(
            [f"{s}\n{nt}" for s, nt in zip(self.sequence, self.structure)]
        )
        ax.set_ylim((0, ymax * 1.25))

    def num_entries(self) -> int:
        """
		Getter that shows how many entries are in the JunctionData object.
		
		:rtype: int
		"""
        return len(self.entries)


    def get_indices(self) -> Dict[str, List[int]]:
        """
		Aggregates original construct indices into a dict() where the (key, value) pairs are
		(construct name, junction indices) with types of (str, List[int]) for each JunctionEntry() 
		entry in the JunctionData() object.

		:rtype: Dict[str,List[int]]
		"""
        result : Dict[str, List[int]] = dict()
        ee : JunctionEntry
        for ee in self.entries:
            result[ee.construct()] = ee.indices()
        return result

    def get_reactivities(self) -> Dict[str, List[float]]:
        """
		Aggregates original construct reactivities into a dict() where the (key, value) pairs are
		(construct name, reactivities) with types of (str, List[float]) for each JunctionEntry() 
		entry in the JunctionData() object.

		:rtype: Dict[str,List[float]]
		"""
        result : Dict[str, List[float]] = dict()
        ee : JunctionEntry
        for ee in self.entries:
            result[ee.construct()] = ee.reactivity()
        return result


    def get_external(self) -> Dict[str, bool]:
        """
		Aggregates whether each junction is an external junction for each of the JunctionEntry() objects
        insdie of the JunctionData() object.

		:rtype: Dict[str,bool]
		"""
        result : Dict[str, List[float]] = dict()
        ee : JunctionEntry
        for ee in self.entries:
            result[ee.construct()] = ee.is_external()
        return result

