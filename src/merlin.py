from linkage2allegro.converter import Converter

class Merlin(Converter):    

    def extractDescent(self, file1):
        self.extractHaplo(file1, True)

    def extractHaplo(self, file1, use_flow=False):
        ped_map = {}

        tmp = [
            None,  # fam id
            [],    # [individuals]
            []     # [assosciated allele pairs]
        ]

        def flushTmpData(tmp):
            # finish populating alleles
            if len(tmp[1]) != len(tmp[2]):
                print("length mismatch", file=stderr)
                exit(-1)

            for tpa in range(len(tmp[1])):

                perc_alleles = tmp[2][tpa]
                perc_id = tmp[1][tpa]
                fam_id = tmp[0]

                if fam_id not in ped_map:
                    ped_map[fam_id] = {}

                if perc_id not in ped_map[fam_id]:
                    ped_map[fam_id][perc_id] = (
                        perc_alleles[0], perc_alleles[1]
                    )

            # clear
            tmp[1] = []
            tmp[2] = []

        with open(file1, 'r') as fio:

            for line in fio:

                if line.startswith("FAMILY"):
                    flushTmpData(tmp)

                    tmp[0] = int(line.split()[1])
                    continue

                if len(tmp[1]) == 0:  # hunt names after a flush
                    if line.find("(") != -1 and line.find(")") != -1:

                        people = [
                            x.strip() for x in line.splitlines()[0].split("  ")
                            if x.strip() != ""
                        ]

                        for p in range(len(people)):
                            perc = people[p].split(" ")
                            perc_id = int(perc[0])

                            # Add new perc to tmp array with blank alleles
                            tmp[1].append(perc_id)
                            tmp[2].append([[], []])
                    #
                    continue

                # Allele pairs
                trimmed = line.strip()
                if len(trimmed) == 0:
                    flushTmpData(tmp)
                    continue

                multiple_alleles = [
                    x.strip() for x in trimmed.split("   ")
                    if x.strip() != ""
                ]

                if len(multiple_alleles) != len(tmp[1]):
                    print("Num alleles and num percs mismatch", file=stderr)
                    exit(-1)

                for a in range(len(multiple_alleles)):
                    alleles = multiple_alleles[a]
                    left_b_right = alleles.split()

                    if not use_flow:
                        # pick first phasing
                        left_b_right[0] = int(
                            left_b_right[0]
                            .split(",")[0]
                            .replace("A", "").replace("?", "0")
                        )
                        left_b_right[2] = int(
                            left_b_right[2]
                            .split(",")[0]
                            .replace("A", "").replace("?", "0")
                        )
                    else:
                        # convert letters to numbers (A - Z)
                        left_b_right[0] = ord(left_b_right[0]) - 64
                        left_b_right[2] = ord(left_b_right[2]) - 64

                    tmp[2][a][0].append(left_b_right[0])
                    tmp[2][a][1].append(left_b_right[2])

            flushTmpData(tmp)

        if use_flow:
            self.descent_map = ped_map
        else:
            self.haplo_map = ped_map

    def extractLOD(self, file1):
        pos_lod = {}

        with open(file1, 'r') as fio:
            line = ""
            while not line.startswith(
                    "       POSITION        LOD      ALPHA       HLOD"
            ):
                line = fio.readline()

            tokens = Converter.tokenizer(fio.readline())

            while len(tokens) == 4:
                gpos = float(tokens[0])

                lod = tokens[1]
                if lod == "-INFINITY":
                    lod = -10000
                lod = float(lod)

                alph = float(tokens[2])
                hlod = float(tokens[3])

                pos_lod[gpos] = (lod, alph, hlod)
                tokens = Converter.tokenizer(fio.readline())

        self.makeLODArray(pos_lod)
