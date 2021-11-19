"""Contains the JunctionEntry and JunctionData objects"""
from __future__ import annotations
import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
from .stats import EXCESS_MAPPER
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict
from rna_library.core import is_symmetrical, InvalidArgument, safe_mkdir
#TODO finish documentaiton here

class JunctionEntry:
    """Represents a single junction entry from an RNA construct"""

    def __init__(self, **kwargs):
        # now do the variable assigment
        self.__sequence = kwargs.get("sequence", None)
        self.__structure = kwargs.get("structure", None)
        self.__reactivity = kwargs.get("reactivity", None)
        self.__construct = kwargs.get("construct", None)
        self.__sn = kwargs.get("sn", None)
        self.__reads = kwargs.get("reads", None)
        self.__score = kwargs.get("score", None)
        self.__symmetrical = is_symmetrical(self.__sequence)
        # argument validation
        self.validate_arguments_()

    def validate_arguments_(self):
        """
        Helper method that validates arguments in the constructor.
        """
        # some basic size checks
        assert len(self.__structure) == len(self.__sequence)
        assert len(self.__sequence) == len(self.__reactivity)
        assert len(self.__sequence)
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

    def __getitem__(self, idx: int) -> float:
        return self.__reactivity[idx]


class JunctionData:
    """
    Composite class that represents a collection of JunctionEntry objects in an experiment.
    """

    def __init__(self, **kwargs):
        # sequence, structure, entries
        self.sequence = kwargs.get("sequence", None)
        self.structure = kwargs.get("structure", None)
        self.entries = kwargs.get("entries", None)
        self.symmetrical = is_symmetrical(self.sequence)
        self.data = [[] for _ in self.sequence]

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


    def get_active_data( self ):
        active = []
        for row, nt in zip( self.data, self.sequence ):
            if nt == 'A' or nt == 'C':
                active.append( row )	
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

    def is_symmetrical(self):
        """
        Getter that tells if the current JunctionData object models a symmetrical junction.
        :rtype: bool
        """
        return self.symmetrical

    def plot(self, plot_dir: str) -> None:
        """
        Method that saves a plot of the JunctionData's data points to the supplied directory. 

        :param: str plot_dir: The directory where the plot will be saved. Does not have to exist.
        :rtype: NoneType
        """
        safe_mkdir(plot_dir)
        fname = f"{plot_dir}/{self.structure}_{self.sequence}.png"
        if os.path.isfile(fname):
            return

        fig, axes = plt.subplots(1, 1, figsize=(15, 10))
        fig.suptitle(name)
        self.bind(axes)
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


    def measure_variance( self ) -> Dict[str,float]:
        
        result = dict()		
        active_data = self.get_active_data()
        for method, func in EXCESS_MAPPER.items():
            temp = 0
            for row in active_data:
                temp += func( row )
            result[ method ] = temp
            
        return result 



