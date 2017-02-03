class Overlap:
    def __init__(self, line):
        v = line.split()
        self.a = v[0]
        self.aRead = [int(v[1]), int(v[2])]
        self.b = v[3]
        self.bRead = [int(v[4]), int(v[5])]
        if v[8] == "-1" or v[9] == "None":
            self.aOvp=[0,0]
            self.bOvp=[0,0]
        else:
            self.aOvp = [int(v[6]), int(v[7])]
            self.bOvp = [int(v[8]), int(v[9])]

        self.indel = int(v[10])
        self.aMidOvp = [int(v[11]), int(v[11])]
        self.bMidOvp = [int(v[12]), int(v[13])]

    def BLength(self):
        return self.bRead[1]-self.bRead[0]

    def ALength(self):
        return self.aRead[1]-self.aRead[0]

    def Print(self, file):
        file.write("\t".join([str(i) for i in [self.a] + self.aRead + [self.b] + self.bRead + self.aOvp + self.bOvp]) + "\n")

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

    def Extends(self, wiggle=500):
        return self.aRead[1] - self.aOvp[1] < wiggle and self.bOvp[0] - self.bRead[0] < wiggle

    def Contained(self, wiggle=500):
        if self.aOvp[0] - self.aRead[0] < wiggle and self.aRead[1] - self.aOvp[1] < wiggle:
            return True
        if self.bOvp[0] - self.bRead[0] < wiggle and self.bRead[1] - self.bOvp[1] < wiggle:
            return True
        return False


def ReadOverlapFile(overlapFileName):

        
    overlapFile = open(overlapFileName)
    overlaps = []
    i=1
    for line in overlapFile:
        # hack to get around blank lines
        if len(line) > 1:
            overlaps.append(Overlap(line))
        i+=1
    
    return overlaps



def MakeOverlapQuery(overlaps):
    query = { (ovp.a, ovp.b) : ovp for ovp in overlaps }
    return query


def GetOverlapPoints(first, second):
    segmentStart = first.bOvp[1]
    segmentEnd = second.aOvp[1]
    return (segmentStart, segmentEnd)
