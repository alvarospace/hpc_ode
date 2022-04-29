import sys
import os

ComplierType = sys.argv[1]

sources = []
CurDir =  os.getcwd()
SourceDir = os.path.join(CurDir,sys.argv[2])

for root, dirs, files in os.walk(SourceDir):
    for file in files:
        if file.endswith(".c"):
            path  = os.path.join(root, file)
            path  = os.path.join(CurDir, path)
            sources.append(path)

for i in sources:
    string = '{} -03 -I{} -I{}/jacobs/ -c {}'.format(ComplierType,SourceDir,SourceDir,i)
    os.system(string)

os.system("ar rc libc_pyjac.a *.o")

