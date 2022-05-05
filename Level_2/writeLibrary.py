import sys
import os

# Argument complier g++/icpc
ComplierType = sys.argv[1]

sources = []
CurDir =  os.getcwd()

# Path to out directory generated from PyJac 
SourceDir = os.path.join(CurDir,sys.argv[2])

# Recurrsively find all sources 
for root, dirs, files in os.walk(SourceDir):
    for file in files:
        if file.endswith(".c"):
            path  = os.path.join(root, file)
            path  = os.path.join(CurDir, path)
            sources.append(path)

# Compile sources  
for i in sources:
    string = '{} -03 -I{} -I{}/jacobs/ -c {}'.format(ComplierType,SourceDir,SourceDir,i)
    os.system(string)

# make library 
os.system("ar rc libc_pyjac.a *.o")

