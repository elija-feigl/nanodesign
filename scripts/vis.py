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

"""This module is used to visualize a 3D DNA structure defined by a caDNAno
   design file.

   Virtual helix and strand information is read from a caDNAno design file and
   used to generate 3D coordinates of virtual helices and their base positions
   based on the design lattice type. These coordinates are then used to
   visualize the geometry of virtual helices as cylinders, scaffold and staple
   strands as they wind through the virtual helices, base coordinates and
   coordinate frames, and domains calculated from base pairing. The backbone P
   atoms and atom bonds of the atomic model generated for the DNA structure can
   also be visualized.

   The 3D geometry of the DNA structure can be visualized by virtual helix or
   strand basis by selecting an ID from a popup menu.  Virtual helices are
   selected from a menu using their caDNAno numbering; strands are selected
   using a name derived from their type (scaffold or staple), the virtual helix
   number they start from and the base position within the virtual helix they
   start from. Strand names are sorted by virtual helix number.

   Selecting an ID from a popup menu displays its geometry in the grapics
   window. Items selected from a menu display a '+' after it to indicate it has
   been selected and is visible. Selecting an item with a '+' will then deselect
   it and remove it from the graphics window. All menus have a 'All' and 'None'
   entries that are used to select all or none of the items in the menu,
   respectively.

   Geometry displayed in the graphics window can be interactively rotated,
   translated and scaled using the mouse. Geometry can also be interactively
   selected using the mouse. The selection is a 3D operation that highlights and
   prints any meaningful component of the geometry. For example, selecting a
   bond on a atom bonds geometry of the atomic model highlights all the bonds
   belonging to the base and prints its strand name, sequence number and base
   name.

   Each menu operation generates a visualization command that is written to a
   file named 'vis.cmd'. These commands can be saved to a file and read in to
   perform operations when the visualizer starts. Commands can also be executed
   from the command line as a semicolon-separated string.
"""
import sys
import logging
import argparse

from nanodesign.visualizer.model import VisModel
from nanodesign.converters.converter import Converter
from nanodesign.converters.pdbcif.atomic_structure import AtomicStructure


def parse_args():
    """ Parse command-line arguments. """
    parser = argparse.ArgumentParser()
    parser.add_argument("-a",  "--atomic_model",
                        help="generate atomic model")
    parser.add_argument("-c",  "--commands",       help="commands")
    parser.add_argument("-cf", "--cmdfile",        help="command file")
    parser.add_argument("-i",  "--infile",         help="input file")
    parser.add_argument("-is", "--inseqfile",      help="input sequence file")
    parser.add_argument("-isn", "--inseqname",      help="input sequence name")
    parser.add_argument("-s",   "--staples",       help="staple operations")
    parser.add_argument("-x",   "--transform",
                        help="apply a transformation to a set of helices")
    return parser.parse_args()


def read_file(args, logger):
    """ Read in a cadnano file. """
    converter = Converter()
    logger = logging.getLogger('nanodesign.visualizer')
    converter.logger = logger

    if args.infile is None:
        logger.error("No input file name given.")
        sys.exit(1)
    else:
        logger.info("Input file name: %s" % args.infile)
        # infile = args.infile

    converter.read_cadnano_file(args.infile, args.inseqfile, args.inseqname)

    return converter


def main():
    logger = logging.getLogger("vis")

    # Get command-line arguments.
    args = parse_args()

    # Read cadnano file and create dna structure.
    converter = read_file(args, logger)
    dna_structure = converter.dna_structure

    # appy a 3D transformation to the geometry of selected helices.
    if args.transform:
        converter.transform_structure(args.transform)

    # perform staple operations (e.g., delete, generate maximal set, etc.)
    if args.staples:
        converter.perform_staple_operations(args.staples)

    # Generate atomic models of the dna structure.
    if not args.atomic_model or (args.atomic_model == "true"):
        atomic_structure = AtomicStructure(dna_structure)
    else:
        atomic_structure = None

    # Initialize visualization.
    vis_model = VisModel(args.infile, args.cmdfile,
                         args.commands, dna_structure, atomic_structure)

    # Start interactive visualization.
    vis_model.start_interactive()


if __name__ == '__main__':
    main()
