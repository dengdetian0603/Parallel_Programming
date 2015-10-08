#!/usr/bin/python

import sys

for line in sys.stdin:
    # take each individual as a vector
    person = line.strip().split()
    keyperson = person[0]
    # find all possible triangles involving key person
    n = person.__len__()
    if n > 2:
        for i in xrange(1,n-1):
            for j in xrange(i+1,n):
                triangle = [int(keyperson), int(person[i]), int(person[j])]
                triangle.sort()
                print triangle[0], triangle[1], triangle[2], 1



