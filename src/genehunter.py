from __future__ import print_function
from converter import Converter

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

                    while not line.startswith("*****"):
                        line = hio.readline()
                        indiv_all1 = map(int, line.split())
                        line = hio.readline()
                        indiv_all2 = map(int, line.split())

                        if indiv_all1 == [] and indiv_all2 == []:
                            return True

                        indiv_id = int(indiv_all1[0])

                        if len(indiv_all1) != len(indiv_all2) + 4:
                            print >> stderr, (
                                "Allele mismatch for indiv",
                                fam_id, indiv_id)
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
            # Jump to marker positions
            line = ""
            while not line.startswith("Current map ("):
                line = fio.readline()
                # iterate

            # # Extract gts
            # while not line.startswith("npl"):
            #     line = fio.readline()
            #     gt_line += line.splitlines()[0] + " "

            # Get num pedigrees
            while not line.startswith("Totalling pedigrees"):
                line = fio.readline()
                # iterate

            # Jump to LOD scoring
            while not line.startswith("position LOD_score  (alpha, HLOD)"):
                line = fio.readline()

            # Extract lod
            line = fio.readline()

            while len(line) != 1:
                # print len(line), line
                (gpos, lod, alpha, hlod,
                 npl, pval, info) = Converter.tokenizer(line)
                gpos = float(gpos)

                if lod == "-INFINITY":
                    lod = -10000
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

                line = fio.readline()
            fio.close()

        self.makeLODArray(pos_lod)
