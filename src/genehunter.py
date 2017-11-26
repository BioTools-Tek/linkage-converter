#!/usr/bin/env python

from __future__ import print_function
from linkage2allegro.converter import Converter

class Genehunter(Converter):

    def extractHaplo(self, hapfile):
        with open(hapfile, 'r') as hio:
            line = "   "
            fam_id = -1

            try:
                while True:
                    line = hio.readline()
                    fam_id = int(line.split()[1])

                    if fam_id not in self.haplo_map:
                        self.haplo_map[fam_id] = {}

                    line = ""
                    max_line = 100000
                    while not line.startswith("*****"):
                        max_line -= 1
                        if max_line < 0:
                            return False

                        line = hio.readline()
                        indiv_all1 = [int(x) for x in line.split()]
                        line = hio.readline()
                        indiv_all2 = [int(x) for x in line.split()]

                        if indiv_all1 == [] and indiv_all2 == []:
                            return True

                        indiv_id = int(indiv_all1[0])

                        if len(indiv_all1) != len(indiv_all2) + 4:
                            print("Allele mismatch for indiv", fam_id, indiv_id, file=stderr)
                            exit(-1)

                        self.haplo_map[fam_id][indiv_id] = (
                            indiv_all1[4:],
                            indiv_all2
                        )

            except IOError:
                hio.close()
                return True

    def extractLOD(self, file):
        pos_lod = {}  # pos -> lod

        with open(file, 'r') as fio:

            line_found = {
                ("position", "NPL_score", "p-value", "information") : False,                                  # NPL
                ("position", "LOD_score", "(alpha,", "HLOD)", "information") : False,                         # LOD
                ("position", "LOD_score", "(alpha,", "HLOD)", "NPL_score", "p-value", "information") : False  # both
            }
            posline_found = False

            for line in fio:
                tokens = tuple(Converter.tokenizer(line))
                
                # Jump to LOD or NPL scoring
                if tokens in line_found:
                    line_found[tokens] = True
                    posline_found = True
                    continue

                # Reached pos line
                if posline_found:

                    lod, alpha, hlod = "-INFINITY", "(0,", "0)"
                    parsed = False

                    if len(tokens) == 4:
                        gpos, lod, pval, info = tokens
                    elif len(tokens) == 5:
                        gpos, lod, alpha, hlod, info = tokens
                    elif len(tokens) == 7:
                        gpos, lod, alpha, hlod, npl, pval, info = tokens
                    else:
                        continue
                        
                    if lod == "-INFINITY":
                        lod = -10000
                        
                    gpos = float(gpos)
                    lod = float(lod)
                    alpha = float(alpha.split("(")[-1].split(",")[0])
                    hlod = float(hlod.split(")")[0])

                    # insert
                    if gpos not in pos_lod:
                        pos_lod[gpos] = [lod, alpha, hlod]
                    else:
                        # update if new lod is larger
                        old_lod = pos_lod[gpos][0]
                        if lod > old_lod:
                            pos_lod[gpos][0] = lod

                #else:if not all_found:
                #    raise Exception("Could not find any linkage data")

            fio.close()

        self.makeLODArray(pos_lod)
