class Overlap:
    def __init__(self, line):
        v = line.split()
        self.a = v[0]
        self.aRead = (int(v[1]), int(v[2]))
        self.b = v[3]
        self.bRead = (int(v[4]), int(v[5]))
        if v[8] == "-1" or v[9] == "None":
            self.aOvp=(0,0)
            self.bOvp=(0,0)
        else:
            self.aOvp = (int(v[6]), int(v[7]))
            self.bOvp = (int(v[8]), int(v[9]))

    def HasOverlap(self):
        if self.aOvp[0] == self.aOvp[1] or self.bOvp[0] == self.bOvp[1]:
            return False
        else:
            return True

    def Overlap(self):
        if self.HasOverlap():
            # determine overlap length from the src
            return self.aOvp[1] - self.aOvp[0]
        else:
            return 0

    def DistanceToStart(self):
        if self.HasOverlap():
            return self.aOvp[0] - self.aRead[0]
        else:
            return -1

    def OverlapLength(self):
        if self.HasOverlap():
            return self.aOvp[1] - self.aOvp[0]
        else:
            return 0

    def Contained(self, wiggle=500):
        if self.aOvp[0] - self.aRead[0] < wiggle and self.aRead[1] - self.aOvp[1] < wiggle:
            return True
        if self.bOvp[0] - self.bRead[0] < wiggle and self.bRead[1] - self.bOvp[1] < wiggle:
            return True
        return False
        

def ReadOverlapFile(overlapFileName):
    ovpFile = open(overlapFileName)
    return [ Overlap(line) for line in ovpFile]

def ReadOverlapFile(overlapFileName):
    ovpFile = open(overlapFileName)
    return [ Overlap(line) for line in ovpFile]


def MakeOverlapQuery(overlaps):
    query = { (ovp.a, ovp.b) : ovp for ovp in overlaps }
    return query

