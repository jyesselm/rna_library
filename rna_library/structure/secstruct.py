from rna_library.structure.motif import *
from rna_library.core import *
from .parser import parse_to_motifs

#TODO documentation
class SecStruct:
    """
    Represents a 
    """

    # what do we want this to do?
    # 1. serve as a place to hold the motifs
    # 2. allow for insertion and changing of motifs

    # dont use named variables when there was only 2 variables
    # switch sequence and secstruct
    def __init__(self, secstruct: str, sequence: str):
        self.root_ = None
        self.structure_ = ""
        self.sequence_ = ""
        self.id_mapping_ = {}
        self.total_structures_ = 1

        assert len(secstruct) == len(sequence) and len(secstruct) > 0
        assert valid_db( secstruct )
        assert len(re.sub("[ACGUTN]", "", sequence)) == 0
        assert secstruct.count("(") == secstruct.count(")")

        self.structure_ = secstruct
        self.sequence_ = sequence
        self.root_ = parse_to_motifs(self.structure_, self.sequence_)
        # helpe for recusion
        self.counter_ = 0
        self.set_ids_(self.root_)

        self.it_ = self.root_
        self.__end_id = highest_id(self.root_)

    def set_ids_(self, m: Motif):
        m.id(self.counter_)
        self.id_mapping_[self.counter_] = m

        self.counter_ += 1

        for child in m.children():
            self.set_ids_(child)

    def display(self):
        print(self.root_.str())

    @property
    def sequence(self):
        return self.sequence_

    @property
    def structure(self):
        return self.structure_

    def helix_replace_(self, id, secstruct, sequence):
        helix = self.id_mapping_[id]
        helix.structure(secstruct)
        helix.sequence(sequence)

    def motif_replace_(self, id, new_secstruct, new_sequence):
        motif = self.id_mapping_[id]
        new_motif = parse_to_motifs(new_secstruct, new_sequence)

        if motif.has_parent():
            # reset the parent of the new motif
            new_motif.parent(motif.parent())
            motif.parent(None)

            cleaned_children = []
            for child in motif.parent().children():
                if child.id() == id:
                    cleaned_children.append(new_motif)
                else:
                    cleaned_children.append(child)

            new_motif.parent().set_children(cleaned_children)

        else:
            # this is the case where you are actually just replacing the root? seems like a dumb thing to do
            # but who knows what people will try these days
            pass

    def change_motif(self, id, new_secstruct, new_sequence):
        # performing some basic checks
        assert len(new_sequence) == len(new_secstruct)
        assert len(re.sub("[AUCGT&]", "", new_sequence)) == 0
        assert len(re.sub("[.)(&]", "", new_secstruct)) == 0
        assert (
            new_secstruct.count("&") == new_sequence.count("&")
            and new_sequence.count("&") <= 1
        )
        assert id in self.id_mapping_

        helix_replace = new_secstruct.find("&") != -1

        if helix_replace:
            self.helix_replace_(id, new_secstruct, new_sequence)
        else:
            self.motif_replace_(id, new_secstruct, new_sequence)

        full_secstruct, full_sequence = (
            self.root_.recursive_structure(),
            self.root_.recursive_sequence(),
        )

        self.structure_ = full_secstruct
        self.sequence_ = full_sequence

        self.root_ = parse_to_motifs(full_secstruct, full_sequence)
        # helpe for recusion
        self.counter_ = 0
        self.id_mapping_ = dict()
        self.set_ids_(self.root_)
        self.total_structures_ += 1

    def get_sequence_structure(self):
        return self.sequence_, self.structure_

    def _get_ids_internal(self, m, ids, mtype):
        if m.type() == mtype:
            ids.append((m.id(), m.token(), m.sequence()))

        for c in m.children():
            self._get_ids_internal(c, ids, mtype)

    def get_ids(self, motif_type):
        ids = []
        for mtype, string in TYPE_MAPPER.items():
            if string == motif_type:
                self._get_ids_internal(self.root_, ids, mtype)
                return ids
        else:
            print(
                f"{motif_type} is an invalid motif type. Allowed values are:\n{', '.join(list(TYPE_MAPPER.values()))}"
            )

    def get_motif(self, id) -> Motif:
        motif = self.id_mapping_[id]
        return motif

    def get_substructure(self, id1, id2=None):
        # get substructure between id1 and id2
        # or get all children of id1
        # make
        pass

    def get(self, id):
        return self.id_mapping_[id]

    def __iter__(self):
        for m in self.id_mapping_.values():
            yield m

    def hairpins(self, **kwargs):
        for m in self.id_mapping_.values():
            if m.is_hairpin():
                yield m

    def helix(self):
        for m in self.id_mapping_.values():
            if m.is_helix():
                yield m

    def junctions(self):
        for m in self.id_mapping_.values():
            if m.is_junction():
                yield m

    def singlestrands(self):
        for m in self.id_mapping_.values():
            if m.is_singlestrand():
                yield m

    def set_barcode(self, m_id, bc_seq):
        m: Motif
        m = self.id_mapping_[m_id]
        if not m.is_singlestrand() and not m.is_helix():
            raise TypeError(f"Barcode motif must be either a singlestrand or helix")

        if not m.same_pattern(bc_seq):
            raise TypeError(
                f'The supplied barcode "{bc_seq}" is not compatible with the selected motif'
            )

        m.sequence(bc_seq)
        self.sequence_ = self.root_.recursive_sequence()

    def __add__(self, other):
        seq_total = self.sequence + other.sequence
        ss_total = self.structure + other.structure
        # TODO bring over barcode info
        return SecStruct(ss_total, seq_total)
