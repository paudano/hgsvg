def ReadPaths(pathsFileName):
    pathsFile = open(pathsFileName)
    paths = []
    for line in pathsFile:
        paths.append(line.rstrip().split())
    return paths
        
