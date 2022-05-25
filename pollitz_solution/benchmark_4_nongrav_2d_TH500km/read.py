import numpy as np
import sys

in_name = sys.argv[1]
n=int(sys.argv[2])   # number of rows you want to read from the file
out = np.zeros(n*3).reshape(n,3)
with open(in_name,'rt') as f:
    for i in range(n):
        f.read(10)
        f.read(10)
        out[i,0] = f.read(13)
        out[i,1] = f.read(13)
        out[i,2] = f.read(13)
        f.readline()
#print(out) see if the output is correct
out_name = sys.argv[3]
np.savetxt(out_name,out)
