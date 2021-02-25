# Copyright 2016 Autodesk Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# 2021.02.25.: moved atom, molecule to seperate file @Elija, Feigl

"""
This module defines the classes used to create an atomic representation of a
DNA structure.

A Molecule object is created for each AtomicStructureStrand object. It contains
a list of atoms for a strand.

"""
from collections import OrderedDict
import numpy as np
import os

# temp code to handle objects as they are being transitioned into  main package

dna_pdb_templates_dir = os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '../../res/'))


class Atom(object):
    """ This class stores the data for an atom.

        Attributes:
            chainID (string): The ID of the chain the atom belongs to.
            coords (NumPy 3x1 array of floats): The atom coordinates.
            element (string): The atom element name (e.g. 'P').
            id (int): The atom ID.
            name (string): The atom name (e.g. 'OP1').
            res_name (string): The name of the residue the atom belongs to
                (e.g. 'DT').
            res_seq_num (int): The residue sequence number.
    """

    def __init__(self, id, name, res_name, chainID, res_seq_num, x, y, z,
                 element):
        self.id = int(id)
        self.res_name = res_name
        self.res_seq_num = int(res_seq_num)
        self.chainID = chainID
        self.name = name
        self.element = element
        self.coords = np.array([float(x), float(y), float(z)], dtype=float)

    def dupe(self):
        """ Create a new Atom object with the same data. """
        return Atom(self.id, self.name, self.res_name, self.chainID,
                    self.res_seq_num, self.coords[0], self.coords[1],
                    self.coords[2], self.element)


class MoleculeType:
    """ This class defines atomic structure types. """
    UNKNOWN = 'unknown'
    ION = 'ion'
    NUCLEIC_ACID = 'nucleic_acid'
    PROTEIN = 'protein'
    WATER = 'water'


class Molecule(object):
    """ This class stores data for a molecular structure.

        Attributes:
            atoms (List[Atom]): The atoms in this molecule.
            chains (Set[string]): The chains in this molecule.
            model_id (int): The ID of the model this molecule belongs to.
    """

    def __init__(self, model_id, type=MoleculeType.UNKNOWN):
        self.id = model_id
        self.model_id = model_id
        self.type = type
        self.atoms = []
        self.chains = set([])
        self.residues = OrderedDict()

    def add_atom(self, atom):
        """ Add an atom to the list of atoms.
        """
        self.atoms.append(atom)
        self.chains.add(atom.chainID)
        if atom.res_seq_num not in self.residues:
            # self.residues[atom.res_seq_num] = []
            self.residues[atom.res_seq_num] = {}
        self.residues[atom.res_seq_num][atom.name.strip()] = atom
