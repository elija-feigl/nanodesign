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

"""
This module is used to write DNA Design viewer JSON files.
"""
import collections
import json
import logging
import numpy as np
from math import pi
from ..cadnano.common import CadnanoLatticeType


class ViewerWriter(object):
    """ The ViewerWriter class writes out a DNA Design viewer JSON file.
    """

    def __init__(self, dna_structure, dna_parameters):
        self.dna_structure = dna_structure
        self.dna_parameters = dna_parameters
        self._logger = logging.getLogger(__name__)

    def write(self, file_name):
        """Write a viewer JSON file.

        Arguments:
            file_name (string): The name of a viewer JSON file to write.

        """

        self._logger.info(
            "Writing DNA Design Viewer JSON file: %s " % file_name)
        dna_structure = self.dna_structure
        dna_structure.compute_aux_data()
        lattice_type = CadnanoLatticeType.names[dna_structure.lattice_type]
        strands_info = self._get_strand_info(dna_structure)
        helices_info = self._get_helices_info(dna_structure)
        domains_info = self._get_domain_info(dna_structure)

        vis_model = {'model_name': dna_structure.name,
                     'lattice_type': lattice_type,
                     'strands': strands_info,
                     'virtual_helices': helices_info,
                     'domains': domains_info
                     }

        with open(file_name, 'w') as outfile:
            json.dump(vis_model, outfile, indent=4, separators=(',', ': '))
            # json.dump(vis_model, outfile)

    def _get_domain_info(self, dna_structure):
        """ Get JSON serialized data for all the domains. """
        domains_info = []
        for domain in dna_structure.domain_list:
            point1, point2 = domain.get_end_points()
            base_info = [base.id for base in domain.base_list]
            if (domain.strand):
                # strand_id = domain.strand.id
                strand_base_ids = [base.id for base in domain.strand.tour]
                start_base_index = strand_base_ids.index(base_info[0])
                end_base_index = strand_base_ids.index(base_info[-1])
            else:
                # strand_id = -1
                start_base_index = -1
                end_base_index = -1
            frame = domain.helix.end_frames[:, :, 0]
            info = {'id': domain.id,
                    'color': domain.color,
                    'strand_id': domain.strand.id,
                    'start_position': [point1[0], point1[1], point1[2]],
                    'end_position': [point2[0], point2[1], point2[2]],
                    'orientation': [frame[0, 2], frame[1, 2], frame[2, 2]],
                    'number_of_bases': len(domain.base_list),
                    'bases': base_info,
                    'start_base_index': start_base_index,
                    'end_base_index': end_base_index,
                    'connected_strand': domain.connected_strand,
                    'connected_domain': domain.connected_domain,
                    # "{:.2f}".format(domain.melting_temperature())
                    'melting_temperature': domain.melting_temperature()
                    }
            domains_info.append(info)

        return domains_info

    def _get_helices_info(self, dna_structure):
        """ Get JSON serialized data for helix objects. """
        self._logger.debug(
            "==================== get helix information ===================")
        helices_info = []
        # Need helices sorted by ID for indexing into helix arryay.
        helix_list = sorted(
            list(dna_structure.structure_helices_map.values()), key=lambda x: x.id)

        for helix in helix_list:
            self._logger.debug(
                "---------- helix  id %d  num %d ----------" % (helix.id, helix.lattice_num))
            point1 = helix.end_coordinates[0]
            point2 = helix.end_coordinates[1]
            length = np.linalg.norm(point1-point2)
            frame = helix.end_frames[:, :, 0]
            domain_ids = sorted(helix.get_domain_ids())
            connectivity_info = self._get_helix_conn_info(helix)
            possible_staple_crossovers = helix.possible_staple_crossovers
            possible_scaffold_crossovers = helix.possible_scaffold_crossovers

            info = {'id': helix.id,
                    'length': length,
                    'helix_distance': self.dna_parameters.helix_distance,
                    'base_pair_rise': self.dna_parameters.base_pair_rise,
                    'start_position': list(point1),
                    'orientation': [frame[0, 2], frame[1, 2], frame[2, 2]],
                    'end_position': list(point2),
                    'scaffold_polarity': helix.scaffold_polarity,
                    'cadnano_info': {'row': helix.lattice_row, 'col': helix.lattice_col,
                                     'num': helix.lattice_num},
                    'domains': domain_ids,
                    'helix_connectivity': connectivity_info,
                    'num_possible_staple_crossovers': len(possible_staple_crossovers),
                    'num_possible_scaffold_crossovers': len(possible_scaffold_crossovers)
                    }

            helices_info.append(info)

        return helices_info

    def _get_strand_info(self, dna_structure):
        """ Get JSON serialized data for strands objetcs.
        """
        self._logger.debug(
            "==================== get strand information p ===================")
        strand_info_list = []
        for strand in dna_structure.strands:
            self._logger.debug("---------- strand %d ----------" % strand.id)
            self._logger.debug("Is scaffold %s" % str(strand.is_scaffold))
            self._logger.debug("Is circular %s" % str(strand.is_circular))
            strand_start_base = strand.tour[0]
            strand_end_base = strand.tour[-1]
            self._logger.debug("Start helix %d  pos %d " %
                               (strand_start_base.h, strand_start_base.p))
            self._logger.debug("End helix %d  pos %d " %
                               (strand_end_base.h, strand_end_base.p))
            base = strand.tour[1]
            self._logger.debug("Next helix %d  pos %d " % (base.h, base.p))

            self._logger.debug("Domains:")
            for domain in strand.domain_list:
                # point1, point2 = domain.get_end_points()
                start_base = domain.base_list[0]
                end_base = domain.base_list[-1]
                self._logger.debug("id %d  vhelix %d  start %d  end %d" % (domain.id, start_base.h, start_base.p,
                                                                           end_base.p))

            domain_ids = [domain.id for domain in strand.domain_list]

            # We are now modifying circular strands to eliminate crossovers
            # at the first strand base; we should not need to modify
            # domain lists anymore.
            # TODO (Davep) Remove this when we are confidant that modifying
            # circular strands solves all of the domain ordering issues.
            # if strand.is_circular:
            #    domain_ids = self._modify_domain_ids(strand, domain_ids)

            base_coords = strand.get_base_coords()
            base_info = []
            for i in range(0, len(strand.tour)):
                base = strand.tour[i]
                coord = base_coords[i]
                base_info.append({'id': base.id, 'coordinates': list(coord), 'sequence': base.seq,
                                  'h': base.h,  'p': base.p
                                  })

            info = {'id': strand.id,
                    'is_scaffold': strand.is_scaffold,
                    'is_circular': strand.is_circular,
                    'number_of_bases': len(strand.tour),
                    'virtual_helices': list(strand.helix_list.keys()),
                    'domains': domain_ids,
                    'bases': base_info,
                    'color': strand.color
                    }

            strand_info_list.append(info)
        return strand_info_list

    def _modify_domain_ids(self, strand, domain_ids):
        """ Modify the list of strand domains IDs so that the order of domain bases
            follows that of its parent strand.
            # TODO (Davep) Remove this when we are confidant that modifying
            # circular strands solves all of the domain ordering issues.
        """
        strand_start_base = self.dna_structure.base_connectivity[strand.tour[0]-1]
        circular_domains = {}
        # Create a list of domains at the start of the strand that are in the
        # strand's start virtual helix.
        start_domains = {}
        start_end_ids = set()
        for domain in strand.domain_list:
            start_base = domain.base_list[0]
            if start_base.h == strand_start_base.h:
                start_domains[start_base.p] = domain
                circular_domains[start_base.p] = domain
                start_end_ids.add(domain.id)
            else:
                break

        # Create a list of domains at the end of the strand that are in the
        # strand's start virtual helix.
        end_domains = {}
        for domain in reversed(strand.domain_list):
            start_base = domain.base_list[0]
            if start_base.h == strand_start_base.h:
                end_domains[start_base.p] = domain
                start_end_ids.add(domain.id)
                circular_domains[start_base.p] = domain
            else:
                break

        if (not start_domains) or (not end_domains):
            return domain_ids

        # Sort the start and ends domains by virtual helix position.
        self._logger.debug("**** Modify domain IDs list ****")
        self._logger.debug("Add start domains to end.")
        domain = start_domains.values()[0]
        if domain.base_list[0].p > domain.base_list[-1].p:
            sorted_domains = sorted(circular_domains.keys(), reverse=True)
        else:
            sorted_domains = sorted(circular_domains.keys())
        self._logger.debug("Sorted domains:")
        # Create domain ID list.
        domain_ids = []
        for domain in strand.domain_list:
            if domain.id not in start_end_ids:
                domain_ids.append(domain.id)

        for p in sorted_domains:
            domain = circular_domains[p]
            start_base = domain.base_list[0]
            end_base = domain.base_list[-1]
            self._logger.debug("id %d  vhelix %d  start %d  end %d" % (domain.id, start_base.h, start_base.p,
                                                                       end_base.p))
            domain_ids.append(domain.id)

        self._logger.debug("Domain IDs %s " % str(domain_ids))
        return domain_ids

    def _get_helix_conn_info(self, helix):
        """ Get the information to write for helix connectivity.
        """
        self._logger.debug(
            "==================== get conn info for helix num %d ===================" % helix.lattice_num)
        dna_structure = self.dna_structure
        lattice = dna_structure.lattice
        num_neigh = lattice.number_of_neighbors
        self._logger.debug("Number of lattice directions: %d " % num_neigh)
        helix_connectivity = helix.helix_connectivity
        possible_staple_crossovers = helix.possible_staple_crossovers
        possible_scaffold_crossovers = helix.possible_scaffold_crossovers
        self._logger.debug("Number of possible staple crossovers: %d " % len(
            possible_staple_crossovers))
        self._logger.debug("Number of possible scaffold crossovers: %d " % len(
            possible_scaffold_crossovers))
        self._logger.debug("Number of connected helices: %d " %
                           len(helix_connectivity))

        # Create an array with size the number of possible helix neighors (4 for square lattice, 3 for honeycomb).
        helix_conn_array = [None]*num_neigh
        ang_inc = pi / num_neigh

        # Iterate over the helix connections filling in the appropriate element in helix_conn_array[].
        for connection in helix_connectivity:
            to_helix = connection.to_helix
            nindex = lattice.get_neighbor_index(
                helix.lattice_row, helix.lattice_col, to_helix.lattice_row, to_helix.lattice_col)
            angle = ang_inc * nindex
            dir = connection.direction
            num_staple = 0
            num_scaffold = 0
            crossovers = connection.crossovers
            for crossover in crossovers:
                if crossover.crossover_base.is_scaf:
                    num_scaffold += 1
                else:
                    num_staple += 1
            self._logger.debug(
                f"Connected helix: num: {to_helix.lattice_num}  row: {to_helix.lattice_row} col:{to_helix.lattice_col}")
            self._logger.debug("    Nindex: %d " % nindex)
            self._logger.debug("    Direction: (%g %g %g) " %
                               (dir[0], dir[1], dir[2]))

            crossover_info = self._get_crossover_info(connection)
            self._logger.debug("    crossover_info %s  " % str(crossover_info))

            conn_info = {'helix_id': to_helix.id,
                         'helix_num': to_helix.lattice_num,
                         'angle': angle,
                         'direction': list(dir),
                         'crossovers': crossover_info
                         }

            helix_conn_array[nindex] = conn_info

        # self._logger.debug(" Helix conn array: %s" % str(helix_conn_array))
        return helix_conn_array

    def _get_crossover_info(self, connection):
        """ Get the helix crossover information to be written to a file.
        """
        from_helix = connection.from_helix
        # to_helix = connection.to_helix
        crossovers = connection.crossovers
        start_pos = from_helix.get_start_pos()
        # Strands that are not base-paired don't have a start position.
        if start_pos is None:
            return []
        self._logger.debug("    Start position: %d " % start_pos)

        # Create a list of staple and scaffold crossovers.
        staples = {}
        scaffolds = {}
        for crossover in crossovers:
            base = crossover.crossover_base
            if base.is_scaf:
                scaffolds[base.p] = crossover
            else:
                staples[base.p] = crossover
        self._logger.debug(
            f"    Number of design crossovers: {len(crossovers)}  staple: {len(staples)}  scaffold: {len(scaffolds)}")

        # Add staple information. Double crossovers will occur in pairs separated by
        # a single position.
        crossover_info = []
        sorted_staples = collections.OrderedDict(sorted(staples.items()))
        self._logger.debug("    Staple crossovers: ")
        pos = sorted_staples.keys()
        self._get_crossover_strand_info(
            start_pos, pos, staples, crossover_info)

        # Add scaffold information.
        sorted_scaffolds = collections.OrderedDict(sorted(scaffolds.items()))
        self._logger.debug("    Scaffold crossovers: ")
        pos = sorted_scaffolds.keys()
        self._get_crossover_strand_info(
            start_pos, pos, scaffolds, crossover_info)

        return crossover_info

    def _get_crossover_strand_info(self, start_pos, pos, crossovers, crossover_info):
        num_pos = len(pos)
        n = 0
        while n < num_pos:
            pos1 = pos[n]
            crossover1 = crossovers[pos1]
            base1 = crossover1.crossover_base
            strand1_id = crossover1.strand.id
            strand1_index = crossover1.strand.get_base_index(base1)
            add_pair = False
            if n != num_pos-1:
                pos2 = pos[n+1]
                if pos2-pos1 == 1:
                    crossover2 = crossovers[pos2]
                    base2 = crossover2.crossover_base
                    strand2_id = crossover2.strand.id
                    strand2_index = crossover2.strand.get_base_index(base2)
                    self._logger.debug(
                        "        n: %d  pos1: %d  pos2: %d " % (n, pos1, pos2))
                    self._logger.debug(
                        "               base1 id %d   base2 id %d " % (base1.id, base2.id))
                    self._logger.debug("               strand1: %d  index: %d " % (
                        strand1_id, strand1_index))
                    self._logger.debug("               strand2: %d  index: %d " % (
                        strand2_id, strand2_index))
                    add_pair = True
                    n += 2

            if not add_pair:
                self._logger.debug(
                    "        n: %d  pos1: %d  str: %d" % (n, pos1, strand1_id))
                self._logger.debug("               strand1: %d  index: %d " % (
                    strand1_id, strand1_index))
                pos2 = None
                strand2_id = None
                strand2_index = None
                n += 1

            info = {"vhelix_base_index": pos1 - start_pos,
                    "first_strand_ID": strand1_id,
                    "first_strand_base_index": strand1_index,
                    "second_strand_ID": strand2_id,
                    "second_strand_base_index": strand2_index
                    }

            crossover_info.append(info)

        return crossover_info
