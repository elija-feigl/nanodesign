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


"""This module is used to interactively visualize a DNA model strand.

   A strand of a DNA nanodesign is an abstraction used to represent ssDNA. It is
   represented as a series of bases. Strands wind their way through the virtual
   helices of a DNA nanodesign.

   A strand ID has the following format:

      <Scaffold|staple>.<vhelixNum>.<startPos>

      where vhelixNUm = the virtual helix number from cadnano
             startPos = the position in the virtual helix of the first base in
                        the strand.
"""
from collections import OrderedDict
import logging
import numpy as np
from .geometry import (VisGeometryAxes, VisGeometryCylinder,
                       VisGeometryPath, VisGeometrySphere, VisGeometryLines, vector_norm)


class VisStrandRepType:
    """ This class defines the strand visualization representation types. """
    UNKNOWN = 'unknown'
    CONNECTORS = 'connectors'
    DOMAINS = 'domains'
    FRAMES = 'frames'
    PATH = 'path'
    TEMPERATURE = 'temperature'


class VisStrand(object):
    """ This class is used to visualize a strand from a DNA design.

        Attributes:
            color (List[float]): The default color of the representation geometry.
            dna_strand (DnaStrand): The structure strand from a DNA structure. This object contains all the data needed
                to visualize a strand.
            dna_structure (DnaStructure): The dna structure derived from a DNA design.
            graphics (VisGraphics): The VisGraphics object.
            id (int): The strand id from 0 to the number of strands in the design - 1.
            name (String): The string representation of the strand name for visualization.
            representations (Dict[List[VisGeometry]): The dictionary storing the list of geometry for a representation.
            tour (List[int]): The list of base IDs defining the stand's path through the DNA design.
            start_helix (int): The number where the strand starts.
            start_pos (int): The base position in the virtual helix where the strand starts.
    """

    def __init__(self, model, graphics, dna_structure, dna_strand):
        """ Initialize a VisStrand object.

            Arguments:
                graphics (VisGraphics): The VisGraphics object.
                dna_strand (DnaStrand): The structure strand from a DNA structure.
                    This object contains all the data needed to visualize a strand.
                dna_structure (DnaStructure): The dna structure derived from a DNA design.
                id (int): The strand id from 0 to the number of strands in the design - 1.
        """
        self.id = dna_strand.id
        self.model = model
        self.graphics = graphics
        self.dna_structure = dna_structure
        self.dna_strand = dna_strand
        self.tour = dna_strand.tour
        self.color = dna_strand.color
        self.color.append(1.0)
        self.representations = {}
        # Set the strand starting helix ann position within that helix.
        tour = dna_strand.tour
        start_base = tour[0]
        self.start_helix = start_base.h
        self.start_pos = start_base.p
        # Set the strand name.
        self.name = VisStrand.get_strand_name(self.dna_strand)
        self._logger = logging.getLogger(__name__ + ":" + str(self.name))
        # Set the methods to create geometry for the different representations.
        self.create_rep_methods = {
            VisStrandRepType.CONNECTORS: self.create_connectors_rep,
            VisStrandRepType.DOMAINS: self.create_domains_rep,
            VisStrandRepType.FRAMES: self.create_frames_rep,
            VisStrandRepType.PATH: self.create_path_rep,
            VisStrandRepType.TEMPERATURE: self.create_temperature_rep
        }

    @staticmethod
    def compare(a, b):
        """ The compare function used to sort strands by helix number and then position. """
        return (a.start_helix != b.start_helix) or (a.start_pos != b.start_pos)

    @staticmethod
    def get_strand_name(dna_strand):
        start_base = dna_strand.tour[0]
        if dna_strand.is_scaffold:
            name = "Scaffold"
        else:
            name = "staple"
        name += "_%d_%d" % (start_base.h, start_base.p)
        return name

    def show(self, rep, show, display=True):
        """ Show the geometry for the given representation.

            Arguments:
                rep (String): The representation name from VisStrandRepType.
                show (bool): If true then show the geometry for the representation, else hide it.
                dispay (bool): If true then redisplay all the graphics geometry.

            If the geometry for the representation has not been created then create and store it.
        """
        if rep not in self.representations:
            self.create_rep_methods[rep]()
        for geom in self.representations[rep]:
            geom.visible = show
        if display:
            self.graphics.display()

    def create_rep(self, rep):
        """ Create and store the geometry for the representation if it has not been created. """
        if rep not in self.representations:
            self.create_rep_methods[rep]()

    def set_color(self, rep, color, display=True):
        """ Set the color for the representation. """
        if len(color) == 3:
            color.append(0.5)
        self.create_rep(rep)
        for geom in self.representations[rep]:
            geom.color[:] = color[:]
        if display:
            self.graphics.display()

    def set_line_width(self, rep, line_width, display=True):
        """ Set the line width for the representation. """
        self.create_rep(rep)
        for geom in self.representations[rep]:
            geom.line_width = line_width
        if display:
            self.graphics.display()

    def print_info(self):
        """ Print strand information. """
        start_base = self.tour[0]
        end_base = self.tour[-1]
        self._logger.info("Scaffold %s  Circular %s " % (
            str(self.dna_strand.is_scaffold), str(self.dna_strand.is_circular)))
        self._logger.info("Start vhelix %d position %d  End vhelix %d position %d " % (start_base.h, start_base.p,
                                                                                       end_base.h, end_base.p))
        self._logger.info("Number of bases %d  Number of domains %d " % (
            len(self.tour), len(self.dna_strand.domain_list)))

    def create_domains_rep(self):
        """ Create the geometry for the strand domain representation. """
        self.representations[VisStrandRepType.DOMAINS] = []
        radius = 0.2
        for i, domain in enumerate(self.dna_strand.domain_list):
            point1, point2 = domain.get_end_points()
            point1 = point1.copy()
            point2 = point2.copy()
            if (domain.strand.is_scaffold):
                point1[2] += radius
                point2[2] += radius
            else:
                point1[2] -= radius
                point2[2] -= radius

            # Create a sphere to visually mark the first domain.
            if i == 0:
                v = vector_norm([point2[i] - point1[i] for i in range(0, 3)])
                s = 0.1
                point3 = point1
                point1 = [point3[i] + s*v[i] for i in range(0, 3)]
                ends_highlight_color = [0.0, 0.0, 0.0, 1.0]
                name = "StrandDomain:%s.%d.%d.start" % (
                    self.name, i, domain.id)
                geom = VisGeometryCylinder(name, radius, point3, point1)
                geom.color = domain.color
                geom.ends_highlight_color = ends_highlight_color
                self.graphics.add_render_geometry(geom)
                self.representations[VisStrandRepType.DOMAINS].append(geom)

            name = "StrandDomain:%s.%d.%d" % (self.name, i, domain.id)
            geom = VisGeometryCylinder(name, radius, point1, point2)
            geom.color = domain.color
            geom.selected_callback = self.select_domains
            self.representations[VisStrandRepType.DOMAINS].append(geom)
            self.graphics.add_render_geometry(geom)

    def select_domains(self, geom, index):
        """ Process strand domain selection.

            Arguments:
                geom (VisGeometry): The geometry selected.
                index (int): The index into the geometry selected.
        """
        # base_conn = self.dna_structure.base_connectivity
        tokens = geom.name.split(".")
        domain_num = int(tokens[1])
        domain_id = int(tokens[2])
        domain = self.dna_structure.domain_list[domain_id]
        num_bases = len(domain.base_list)
        start_base = domain.base_list[0]
        end_base = domain.base_list[-1]
        self._logger.info("Selected Strand %s domain. Location in strand %d" % (
            self.name, domain_num+1))
        self._logger.info("Domain ID %d  Number of bases %d  Start pos %d  End pos %d" % (domain_id, num_bases,
                                                                                          start_base.p, end_base.p))
        self.print_info()

    def create_temperature_rep(self):
        """ Create the geometry for the strand melting temperature representation.

            The melting temperature from domains is determined for base positions along a strand.
            The strand geometry is displayed as a sequence of lines with endpoints colored by mapping
            temperatures to a color spectrum.
        """
        # Create the strand path vertices.
        base_coords = self.dna_strand.get_base_coords()
        self._logger.debug("Number of point %d" % len(base_coords))
        verts = []
        for coords in base_coords:
            verts.append([coords[0], coords[1], coords[2]])
        if self.dna_strand.is_circular:
            coords = base_coords[0]
            verts.append([coords[0], coords[1], coords[2]])

        # Create temperature colors for each strand vertex (base position).
        base_domain_map = OrderedDict()
        for base in self.dna_strand.tour:
            base_domain_map[(base.h, base.p)] = None

        for domain in self.dna_strand.domain_list:
            for base in domain.base_list:
                base_domain_map[(base.h, base.p)] = domain

        spectrum_colors = self.graphics.get_spectrum_colors()
        tmin, tmax = self.model.get_domains_temperature_range()
        colors = []
        temp_data = []
        for domain in base_domain_map.values():
            temp = domain.melting_temperature()
            temp_data.append(temp)
            if temp == -500.0:
                color = [0.5, 0.5, 0.5]
            else:
                color = self.graphics.map_value_to_color(
                    spectrum_colors, tmin, tmax, temp)
            color.append(1.0)
            colors.append(color)

        if self.dna_strand.is_circular:
            colors.append(colors[0])

        # Create the strand path geometry.
        name = "StrandTemperature:%s" % self.name
        show_verts = False
        show_arrows = True
        geom = VisGeometryPath(name, verts, show_verts,
                               show_arrows, colors=colors)
        geom.start_marker = True
        geom.select_vertex = True
        geom.line_width = 3.0
        geom.entity_indexes = range(0, len(base_coords))
        geom.data = temp_data
        geom.selected_callback = self.select_temperature
        self.representations[VisStrandRepType.TEMPERATURE] = [geom]
        self.graphics.add_render_geometry(geom)

    def select_temperature(self, geom, index):
        """ Process strand domain selection.

            Arguments:
                geom (VisGeometry): The geometry selected.
                index (int): The index into the geometry selected.
        """
        base = self.tour[index]
        self._logger.info("Selected Strand %s. Domain melting temperature %g " % (
            self.name, geom.data[index]))
        self._logger.info("Location in path %d  Vhelix %d  Position %d " % (
            index+1, base.h, base.p))
        self.print_info()

    def create_frames_rep(self):
        """ Create the geometry for the strand coordinate frames representation. """
        self._logger.debug("Create strand frames rep.")
        # base_conn = self.dna_structure.base_connectivity
        base_coords = self.dna_strand.get_base_coords()
        # helix_map = self.dna_structure.structure_helices_map
        self._logger.debug("Number of base coordinates %d" % len(base_coords))
        self._logger.debug("Tour size %d" % len(self.tour))
        self._logger.debug("Bases: ")
        tour_size = len(self.tour)
        origins = np.empty(shape=(tour_size, 3), dtype=float)
        directions = np.empty(shape=(3, 3, tour_size), dtype=float)
        for i, base in enumerate(self.tour):
            frame = base.ref_frame
            origins[i, :] = base_coords[i]
            directions[:, :, i] = frame

        scale = 0.2
        name = "StrandFrames:%s" % self.name
        geom = VisGeometryAxes(name, origins, directions, scale)
        geom.color = self.color
        geom.entity_indexes = range(0, len(base_coords))
        geom.selected_callback = self.select_frames
        self.representations[VisStrandRepType.FRAMES] = [geom]
        self.graphics.add_render_geometry(geom)

    def select_frames(self, geom, index):
        """ Process strand frames selection.

            Arguments:
                geom (VisGeometry): The geometry selected.
                index (int): The index into the geometry selected.
        """
        base = self.tour[index]
        self._logger.info(
            f"Selected Strand {self.name}. Location in strand {index+1} Base id {base.id} h {base.h} p {base.p}")
        self.print_info()

    def create_path_rep(self):
        """ Create the geometry for the strand geometry representation. """
        base_coords = self.dna_strand.get_base_coords()
        self._logger.debug("Number of point %d" % len(base_coords))
        name = "Strand:%s" % self.name
        # show_verts = True
        show_verts = False
        show_arrows = True
        # If the strand is a scaffold then offset it a bit in Z.
        if self.dna_strand.is_scaffold:
            offset = 0.3
        else:
            offset = 0.0
        verts = []
        for coords in base_coords:
            verts.append([coords[0], coords[1], coords[2]+offset])
        if self.dna_strand.is_circular:
            coords = base_coords[0]
            verts.append([coords[0], coords[1], coords[2]+offset])
        geom = VisGeometryPath(name, verts, show_verts, show_arrows)
        geom.color = self.color
        geom.start_marker = True
        geom.select_vertex = True
        geom.line_width = 2.0
        geom.entity_indexes = range(0, len(base_coords))
        geom.selected_callback = self.select_path
        self.representations[VisStrandRepType.PATH] = [geom]
        self.graphics.add_render_geometry(geom)

    def select_path(self, geom, index):
        """ Process strand path selection.

            Arguments:
                geom (VisGeometry): The geometry selected.
                index (int): The index into the geometry selected.
        """
        base = self.tour[index]
        self._logger.info(
            f"Selected Strand {self.name} path. Location in path {index+1} Vhelix {base.h}  Position {base.p}")
        self.print_info()

    def create_connectors_rep(self):
        """ Create the geometry for the strand connectors representation. """
        # dist_bp = self.dna_structure.dna_parameters.base_pair_rise
        max_dist = 2.0*self.dna_structure.dna_parameters.helix_distance
        points1 = []
        points1_map = set()
        points2 = []
        points2_map = set()
        tour = self.dna_strand.tour
        for i in range(0, len(tour)-1):
            base1 = tour[i]
            pt1 = base1.coordinates
            base2 = tour[i+1]
            pt2 = base2.coordinates
            dist = np.linalg.norm(pt1 - pt2)
            if dist > max_dist:
                self._logger.debug(" ")
                self._logger.debug(
                    f"Base 1  h {base1.h}  p {base1.p}  ({base1.coordinates})")
                self._logger.debug(
                    f"Base 2  h {base2.h}  p {base2.p}  ({base2.coordinates})")
                if (base1.h, base1.p) not in points1_map:
                    points1_map.add((base1.h, base1.p))
                    points1.append(pt1)

                if (base2.h, base2.p) not in points2_map:
                    points2_map.add((base2.h, base2.p))
                    points2.append(pt2)

        verts = []
        for i in range(0, len(points1)):
            verts.append(points1[i])
            verts.append(points2[i])

        arrows = False
        arrows = True
        name = "StrandConnectors:" + str(self.id)
        geom = VisGeometryLines(name, verts, arrows)
        geom.line_width = 2.0
        geom.color = [0.8, 0.0, 0.8, 1]
        # geom.selected_callback = self.select_path
        self.representations[VisStrandRepType.CONNECTORS] = [geom]
        self.graphics.add_render_geometry(geom)

        if False:
            # if True:
            radius = 0.1
            for pt in points1:
                name = "StrandConnectorsPt1:%s" % self.name
                sphere = VisGeometrySphere(name, pt, radius)
                sphere.color = [0.0, 0.8, 0.0, 1]
                self.representations[VisStrandRepType.CONNECTORS].append(
                    sphere)
                self.graphics.add_render_geometry(sphere)
            for pt in points2:
                name = "StrandConnectorsPt2:%s" % self.name
                sphere = VisGeometrySphere(name, pt, radius)
                sphere.color = [0.8, 0.0, 0.0, 1]
                self.representations[VisStrandRepType.CONNECTORS].append(
                    sphere)
                self.graphics.add_render_geometry(sphere)

        if len(self.dna_structure.connector_points) != 0:
            radius = 0.1
            points1 = self.dna_structure.connector_points[0]
            points2 = self.dna_structure.connector_points[1]
            if False:
                verts = []
                for i in range(0, len(points1)):
                    if i == 1:
                        verts.append(points1[i])
                        verts.append(points2[i])

                arrows = True
                arrows = False
                name = "StrandConnectors:" + str(self.id)
                geom = VisGeometryLines(name, verts, arrows)
                geom.line_width = 3.0
                geom.color = [0.8, 0.0, 0.8, 1]
                # geom.selected_callback = self.select_path
                self.representations[VisStrandRepType.CONNECTORS] = [geom]
                self.graphics.add_render_geometry(geom)

            for pt in points1:
                name = "StrandConnectorsPt1:%s" % self.name
                sphere = VisGeometrySphere(name, pt, radius)
                sphere.color = [0.0, 0.8, 0.0, 1]
                self.representations[VisStrandRepType.CONNECTORS].append(
                    sphere)
                self.graphics.add_render_geometry(sphere)
            for pt in points2:
                name = "StrandConnectorsPt2:%s" % self.name
                sphere = VisGeometrySphere(name, pt, radius)
                sphere.color = [0.8, 0.0, 0.0, 1]
                self.representations[VisStrandRepType.CONNECTORS].append(
                    sphere)
                self.graphics.add_render_geometry(sphere)
