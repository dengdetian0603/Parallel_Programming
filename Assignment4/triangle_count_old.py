"""
Author: DETIAN DENG
JHED: ddeng3
Name: triangle_count.py

Get the list of cycle triangles in the graph
"""

from pyspark import SparkContext
from time import time
import sys, os

#TODO: Possibly define functions up here
from operator import add


def GetTriad(people):
    person = people.split()
    keyperson = person[0]
    print person
    allTriads = []
    n = len(person)
    if n > 2:
        for i in xrange(1,n-1):
            for j in xrange(i+1,n):
                triangle = [int(keyperson), int(person[i]), int(person[j])]
                triangle.sort()
                key = str(triangle[2])+" "+str(triangle[1])+" "+str(triangle[0])
                allTriads.append((key,1))
    return allTriads  



# NOTE: Do not change the name/signature of this function
def count_triangles(data, master="local[2]"):
    """
    @brief: Count triangles using Spark
    @param data: The data location for the input files
    @param master: The master URL as defined at
    https://spark.apache.org/docs/1.1.0/submitting-applications.html#master-urls
    """

    #################  NO EDITS HERE ###################
    assert not os.path.exists("triangles.out"), "File: triangles.out \
    already exists"
    sc = SparkContext(master, "Triangle Count")
    start = time()
    ###############  END NO EDITS HERE  ################
    # TODO: Your code goes here!
    people = sc.textFile(data)
    triad = people.flatMap(GetTriad).reduceByKey(add).filter(lambda x: x[1]>1)
    #triadCount = triad.map(lambda x: (x,1))
    #triadSum = triadCount.reduceByKey(add)
    #triangles = triadSum.filter(lambda x: x[1]>1)
    #output = triangles.collect()
    output = triad.collect()
    #triangles.saveAsTextFile("test1")
    #################  NO EDITS HERE  ###################
    print "\n\n*****************************************"
    print "\nTotal algorithm time: %.4f sec \n" % (time()-start)
    print "*****************************************\n\n""" 
    ###############  END NO EDITS HERE ################
    with open("triangles.out", "wb") as f:
        for friends in output:
            f.write(friends[0]+"\n") # TODO: Loop with f to write your result to file serially
        pass


#################  NO EDITS HERE  ###################
if __name__ == "__main__":
    if len(sys.argv) == 2:
        print "Counting triangles with master as 'local[2]'"
        count_triangles(sys.argv[1])
    elif len(sys.argv) == 3: 
        print "Counting triangles with master as '%s'" % sys.argv[2]
        count_triangles(sys.argv[1], sys.argv[2])
    else:
        sys.stderr.write("\nusage: SPARK_ROOT/bin/spark-submit \
            example/python/tri_count.py data_dir [master-url]")
        exit(1)
############### NO EDITS BELOW EITHER ################