from __future__ import print_function
from converter import Converter

class Swiftlink(Converter):

    def extractLOD(self, lodfile):
        pos_lod = {}
        with open(lodfile, 'r') as lio:
            lio.readline()

            for line in lio:
                if line[0] == "-":
                    pos, lod = map(float, Converter.tokenizer(line)[1:])
                    pos_lod[pos] = (lod, 1, 0)

        self.makeLODArray(pos_lod)
