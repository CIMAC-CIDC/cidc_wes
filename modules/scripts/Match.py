#!/usr/bin/env python

import sys
Input_file = sys.argv[1]
f = open(Input_file, 'r')
lines = f.readlines()

for line in lines:
    #print(line)
    if "non-matching" in line:
        print(line)
        a = int(line.split()[1])
    if "main file" in line:
        print(line)
        b = int(line.split()[1])
    if "second file" in line:
        print(line)
        c = int(line.split()[1])
if sum((a,b,c)) == 0:
    print ("match")
else:
    print("mismatch")
    sys.exit()



