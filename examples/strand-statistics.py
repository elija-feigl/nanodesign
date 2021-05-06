#!/usr/bin/env python

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

"""This module demonstrates using the Nanodesign interface to calculate some statistics
   for the strands in a design.

    Using the length of each strand the number of instances of each strand length
    are counted. The counts for each length are then sorted and printed.

    The number of helices each strand visits is then calculated using its list of bases.
    The counts for each helix are then sorted and printed.
"""
from collections import OrderedDict
from nanodesign.converters.converter import Converter


def read_file(file_name, seq_name):
    """ Read in a cadnano file. """
    converter = Converter()
    seq_file = None
    converter.read_cadnano_file(file_name, seq_file, seq_name)
    return converter


def main():
    # Set caDNAno file name.
    file_name = "../tests/samples/Nature09_squarenut_no_joins.json"

    # Set sequence to assign to scaffold.
    seq_name = "M13mp18"

    # Read cadnano file and create dna structure.
    converter = read_file(file_name, seq_name)
    dna_structure = converter.dna_structure

    # Count the number of instances of each strand length.
    strand_lengths = OrderedDict()
    for strand in dna_structure.strands:
        num_bases = len(strand.tour)
        if num_bases not in strand_lengths:
            strand_lengths[num_bases] = 0
        strand_lengths[num_bases] += 1

    print("\nStrand length counts:")
    for length, count in strand_lengths.items():
        print(f'Length {length}  Count {count}')

    # Calculate the number of helices each strand visits.
    strand_helix = OrderedDict()
    for strand in dna_structure.strands:
        helix_ids = set()
        for base in strand.tour:
            helix_ids.add(base.h)
        num_helices = len(helix_ids)
        if num_helices not in strand_helix:
            strand_helix[num_helices] = 0
        strand_helix[num_helices] += 1

    print("\nStrand helix counts:")
    for num_helices, count in strand_helix.items():
        print(f'Number of helices {num_helices}  Count {count}')


if __name__ == '__main__':
    main()
