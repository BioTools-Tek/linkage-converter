#!/usr/bin/env python

from __future__ import print_function
from linkage2allegro.converter import Converter

class Simwalk(Converter):

    def extractHaplo(self, hef):
        self.extractHaploAndDescent(hef)

    def extractDescent(self, hef):
        self.extractHaploAndDescent(hef)

    def extractHaploAndDescent(self, hef):

        # clear marker map
        self.pos_marker = {}
        self.marker_order = []

        hio = open(hef, 'r')

        # extract markerloci
        line = ""
        while not line.startswith(" Name    Female   Male   Alleles"):
            line = hio.readline()
        line = hio.readline()

        while not line.startswith("Pedigree Name,       ID   Father"):
            line = hio.readline()
            if not line.startswith("  "):
                tokens = Converter.tokenizer(line)

                if len(tokens) != 4:
                    break

                gpos = float(tokens[2])
                marker = tokens[0].strip()
                self.pos_marker[gpos] = marker
                self.marker_order.append(marker)

        tmp = {
            "_fam": None,
            "_perc": None,
            "_allpat": [],  # alleles
            "_allmat": [],  #
            "_decpat": [],  # descent
            "_decmat": []
        }

        def insertDat():
            if len(tmp["_allpat"]) > 0:
                fid = tmp["_fam"]
                pid = tmp["_perc"]

                all1 = tmp["_allpat"]
                all2 = tmp["_allmat"]
                dec1 = tmp["_decpat"]
                dec2 = tmp["_decmat"]

                if fid not in self.haplo_map:
                    self.haplo_map[fid] = {}
                    self.descent_map[fid] = {}

                self.haplo_map[fid][pid] = (all1, all2)
                self.descent_map[fid][pid] = (dec1, dec2)

                tmp["_perc"] = None
                tmp["_allmat"] = []
                tmp["_allpat"] = []
                tmp["_decmat"] = []
                tmp["_decpat"] = []

        dashedlines_found = False

        for line in hio:
            if line.startswith("____"):
                if tmp["_perc"] is not None:
                    insertDat()

                dashedlines_found = True
                continue

            if dashedlines_found and not line.startswith(" "):
                fam = line.split("(")[0].strip()
                tmp["_fam"] = int(fam)
                dashedlines_found = False
                continue

            tokens = line.splitlines()[0].split()

            if tmp["_fam"] is not None and len(tokens) == 5:
                insertDat()
                tmp["_perc"] = int(tokens[0])

            # Allele data
            if (
                    tmp["_fam"] is not None and
                    tmp["_perc"] is not None and
                    len(tokens) == 6
            ):
                tokens = map(int, tokens)

                tmp["_allpat"].append(tokens[0])
                tmp["_allmat"].append(tokens[1])
                tmp["_decpat"].append(tokens[2])
                tmp["_decmat"].append(tokens[3])

        hio.close()
        insertDat()

    def extractLOD(self, scorefile):
        # This code is terrible, but simwalk is a horrible format

        pos_marker_lod = {}
        sio = open(scorefile, 'r')

        line = sio.readline()
        while not line.startswith("  NAME  , Haldane cM ,   alpha=1.00   ,"):
            line = sio.readline()

        sio.readline()  # skip empty

        recomb_fracts = False
        map_headers = False

        current_marker = None

        # Now add offsets and calc absolute cM
        pos_lod = {}

        for line in sio:
            if not recomb_fracts:
                if line[:2] == "__":
                    recomb_fracts = True
                    continue

                tokens = [x for x in line.split() if x != ","]
                if len(tokens) == 1 and not line.startswith(" "):
                    current_marker = tokens[0].splitlines()[0].strip()
                    if current_marker in pos_marker_lod:
                        print("Duplicate", current_marker, file=stderr)
                        exit(-1)

                    pos_marker_lod[current_marker] = {}
                    continue

                if len(tokens) != 4:
                    continue

                gpos, lod, hlod, alpha = map(float, tokens)
                if gpos < 0 or gpos > 5:
                    continue

                pos_marker_lod[current_marker][gpos] = (lod, 1, hlod)

            else:
                if line.startswith(
                        "Haldane cM   NAME     FRACTION    OBSERVED & EXPECTED"
                ):
                    map_headers = True
                    continue

                if not map_headers:
                    continue

                if line[:2] == "__":
                    break

                tokens = Converter.tokenizer(line)
                if len(tokens) == 2:
                    gpos, marker = tokens
                    gpos = float(gpos)
                    self.pos_marker[gpos] = marker

                    offset_map = pos_marker_lod[marker]

                    for relative_gpos in offset_map:
                        abs_gpos = gpos + relative_gpos
                        if abs_gpos in pos_lod:
                            print("Duplicate abs cM", abs_gpos, file=stderr)
                            exit(-1)

                        pos_lod[abs_gpos] = offset_map[relative_gpos]

        sio.close()
        self.makeLODArray(pos_lod)
