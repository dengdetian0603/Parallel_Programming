#!/usr/bin/python 

import sys

triDict = dict()

for line in sys.stdin:
    triangle = line.strip().split()
    # key is the tuple of triangle, value is the count
    triDict[(triangle[0], triangle[1], triangle[2])] = triDict.get((triangle[0], triangle[1], triangle[2]), 0) + 1

for key in triDict.keys():
    if triDict[key] > 1 :
    	a, b, c = key
        print a, b, c
        print b, a, c
        print c, a, b

