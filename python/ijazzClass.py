
class ijazzResOP:
    def __init__( self, inputfile, color, legend, isRef = False, isMC = False, lumi = -1 ):
        self.inputfile = inputfile
        self.color     = color
        self.legend    = legend
        self.isRef     = isRef
        self.isMC      = isMC
        self.lumi      = lumi
