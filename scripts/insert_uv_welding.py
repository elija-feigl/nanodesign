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

#
# This file was originally contributed by Jean-Philippe Sobczak of Tilibit Nanosystems
# modified by Elija Feigl Dietzlab TUM Munich 2022

import sys
from pathlib import Path
from typing import List
import logging
from nanodesign.converters.converter import Converter
from nanodesign.converters.dna_sequence_data import dna_sequence_data

logger = logging.getLogger(__name__)


def process_input(argv):
    if len(argv) < 3:
        logger.error("%s CADNANO_FILE SEQUENCE [fill_character]", argv[0])

    # Create path to caDNAno file to load
    file_path: Path = Path(sys.argv[1]).absolute()

    # Set sequence to assign to scaffold.
    scaffold_seq_source: str = sys.argv[2]
    if Path(scaffold_seq_source).exists():
        seq_file_name = str(Path(scaffold_seq_source).absolute())
        seq_name = None
    elif scaffold_seq_source in dna_sequence_data:
        seq_file_name = None
        seq_name = scaffold_seq_source
    else:
        raise FileNotFoundError

    # Set fill character
    if len(argv) == 4:
        fill_char = argv[3]
    else:
        fill_char = "T"

    return file_path, (seq_file_name, seq_name), fill_char


def welding(dna_structure, fill_char: str) -> List[str]:
    """determine scaffold strand id, alert if multiple scaffold strands
    prints starting location if not circular strand
    strand id of last found scaffold strand
    """
    scaffold_id = -1
    for strand in dna_structure.strands:
        if strand.is_scaffold:
            if scaffold_id != -1:
                logger.error("multiple scaffold strands!")
            scaffold_id = strand.id
            if not strand.is_circular:
                logger.info(
                    "scaffold strand starts in helix %s base number %s",
                    strand.tour[0].h,
                    strand.tour[0].p,
                )
            else:
                logger.info(
                    "scaffold strand is circular, passes through helix %s base number %s",
                    strand.tour[0].h,
                    strand.tour[0].p,
                )

    staple_end_sequence_insert = fill_char
    crossover_sequence_insert = 2 * fill_char

    output_lines = list()
    output_lines.append("id, h-p, sequence\n")
    # print double-stranded domain lengths of all staples and starting locations of staples
    for strand in dna_structure.strands:
        if not strand.is_scaffold:

            # new staple sequence with inserted T's
            # insert T in the beginning if no single-stranded overhang
            if strand.domain_list[0].connected_domain != -1:
                modified_sequence = staple_end_sequence_insert
            else:
                modified_sequence = ""

            # add crossover T's between all domains,
            #   if domain and next domain are not single-stranded
            # also check if it is a staple or a scaffold crossover,
            #   do not add T's in scaffold crossovers
            for index, domain in enumerate(strand.domain_list[:-1]):
                if (
                    (domain.connected_domain != -1)
                    and (strand.domain_list[index + 1].connected_domain != -1)
                    and (domain.helix != strand.domain_list[index + 1].helix)
                ):

                    modified_sequence = (
                        modified_sequence + domain.sequence + crossover_sequence_insert
                    )
                else:
                    modified_sequence = modified_sequence + domain.sequence

            # insert T in the end if no single-stranded overhang
            if strand.domain_list[-1].connected_domain != -1:
                modified_sequence = (
                    modified_sequence
                    + strand.domain_list[-1].sequence
                    + staple_end_sequence_insert
                )
            else:
                modified_sequence = modified_sequence + strand.domain_list[-1].sequence

            line = f"{strand.id}, {strand.tour[0].h}-{strand.tour[0].p}, {modified_sequence}\n"
            output_lines.append(line)
            logger.info(line)
    return output_lines


def main():
    # parse command line arguments
    logger.info("Starting UV Welding modifications")
    file_path, (seq_file_name, seq_name), fill_char = process_input(sys.argv)

    # Read cadnano file and create dna structure.
    converter = Converter()
    converter.modify = True
    if not file_path.exists():
        raise FileNotFoundError
    converter.read_cadnano_file(
        file_name=str(file_path),
        seq_file_name=seq_file_name,
        seq_name=seq_name,
    )
    dna_structure = converter.dna_structure
    dna_structure.get_domains()

    output_lines = welding(dna_structure, fill_char=fill_char)

    with open("UVw_sequence.csv", mode="w") as seq_file:
        seq_file.writelines(output_lines)


if __name__ == "__main__":
    main()
