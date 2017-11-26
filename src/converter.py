#!/usr/bin/env python
from __future__ import print_function

from sys import stderr
from bisect import bisect_left
from sys import stderr

class Converter:
    """
    Superclass for converting formats
    """
    haplo_individual_buffer = "%-8d %8d %8d %8d %2d %2d  "
    lod_line_buffer = "%8s %12s %8s %8s  %-s"

    def __init__(self, pedin, mapin):
        self.out_lod = "linkage.allegro_lod"
        self.out_haplo = "linkage.allegro_haplo"
        self.out_descent = "linkage.allegro_descent"

        self.haplo_map = {}  # fam_id -> indiv_id -> (all1,all2)
        self.descent_map = {}
        self.pos_marker = {}  # genpos -> marker
        self.marker_order = []
        self.lod_array = []
        self.pedigree = {}

        self.__populatePedigree(pedin)
        self.__populateMarkerMap(mapin)

    def __populatePedigree(self, input_ped):
        with open(input_ped, "r") as pio:
            for line in pio:
                tokens = [int(x) for x in line.split()]
                f_id, p_id, father, mother, gender, affect = tokens[:6]

                if f_id not in self.pedigree:
                    self.pedigree[f_id] = {}
                if p_id not in self.pedigree:
                    self.pedigree[f_id][p_id] = (
                        father, mother, gender, affect
                    )
                else:
                    print("Duplicate individual:", f_id, p_id, file=stderr)
                    exit(-1)

    def __populateMarkerMap(self, mapin):
        with open(mapin, "r") as mio:
            mio.readline()  # chomp header

            for line in mio:
                chrom, gpos, marker, ppos, nr = Converter.tokenizer(line)
                marker = marker.strip()
                self.pos_marker[float(gpos)] = marker
                self.marker_order.append(marker)

    def __annotateClosestMarker(self, pos_lod):

        # Produce sorted list of gpos
        pos_rsid_array = sorted(self.pos_marker.keys())
        pos_lod_array = sorted(pos_lod.keys())

        keys_to_annotate = {}

        # bisect and find closest marker
        for marker_pos in pos_rsid_array:
            marker = self.pos_marker[marker_pos]
            lm_ind = bisect_left(pos_lod_array, marker_pos)
            lm_pos = pos_lod_array[lm_ind - 1]

            diff1 = marker_pos - lm_pos
            diff1 = diff1 if diff1 > 0 else -diff1

            if diff1 > 0.1:
                continue

            if lm_pos not in keys_to_annotate:
                keys_to_annotate[lm_pos] = marker
            else:
                diff2 = marker_pos - lm_pos
                diff2 = diff2 if diff2 > 0 else -diff2

                if diff1 < diff2:
                    # update with closer marker
                    keys_to_annotate[lm_pos] = marker

        return keys_to_annotate, pos_lod_array

    def makeLODArray(self, pos_lod):
        (
            postns_with_closest_markers, sorted_poslod
        ) = self.__annotateClosestMarker(pos_lod)

        # update map with marker info
        for pos in sorted_poslod:
            marker = ""
            if pos in postns_with_closest_markers:
                marker = postns_with_closest_markers[pos]

            lod, alpha, hlod = pos_lod[pos]
            self.lod_array.append((pos, lod, alpha, hlod, marker))

    def __generateHeaders(self, npad_left=10):
        marker_order = self.marker_order
        max_len = -1
        markerpadd = []
        for marker in marker_order:
            nmark = len(marker)
            if nmark > max_len:
                max_len = nmark

        # paddleft
        for marker in marker_order:
            markerpadd.append(("%%-%ds" % max_len) % marker)

        # transpose
        buffer_left = ("%%%ds" % npad_left) % " "
        return '\n'.join(
            [
                buffer_left + "  ".join(x)
                for x in zip(*markerpadd)
            ][::-1]
        )

    @staticmethod
    def tokenizer(line):
        return line.splitlines()[0].split()

    # Override
    def extractLOD(self, file1):
        return ""

    # Override
    def extractHaplo(self, file1, extrafile=None):
        return ""

    # Override
    def extractDescent(self, file1):
        return ""

    def writeLOD(self):

        out_lod = open(self.out_lod, "w")

        header = Converter.lod_line_buffer % (
            "location", "LOD", "alpha", "HLOD", "marker"
        )
        print(header, file=out_lod)

        for pos, lod, alpha, hlod, marker in self.lod_array:
            out_line = Converter.lod_line_buffer % (
                "%.4f" % pos,
                "%.4f" % lod,
                "%.4f" % alpha,
                "%.4f" % hlod,
                marker
            )
            print(out_line, file=out_lod)
        out_lod.close()
        print("Wrote: ", out_lod.name, file=stderr)

    def writeDescent(self):
        self.writeHaplo(True)

    def writeHaplo(self, descent=False):
        out_file = open(self.out_haplo
                        if not descent
                        else self.out_descent, "w")

        map_map = self.haplo_map if not descent else self.descent_map

        dummy_l = Converter.haplo_individual_buffer % (
            1, 1, 1, 1, 1, 1
        )  # *range(6) (unpack in py3 only...)
        headers = self.__generateHeaders(len(dummy_l)) + '\n'
        print(headers, file=out_file)

        for fam_id in map_map:
            for indiv_id in map_map[fam_id]:
                alleles = map_map[fam_id][indiv_id]
                ped_data = self.pedigree[fam_id][indiv_id]
                father, mother, gender, affect = ped_data

                indiv_data = Converter.haplo_individual_buffer % (
                    fam_id, indiv_id, father, mother, gender, affect
                )

                print(indiv_data + " ".join(map(
                    lambda x: "%-2d" % x, alleles[0]
                )), file=out_file)

                print(indiv_data + " ".join(map(
                    lambda x: "%-2d" % x, alleles[1]
                )), file=out_file)
        out_file.close()

        print("Wrote: ", out_file.name, file=stderr)
