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
def makepair(text):
    words = text.split()
    k = words[0]
    v = set(words[1:])
    return (k,v)


def AintersectB(k1,k2,TheList):
    A = TheList.get(k1)
    B = TheList.get(k2)
    if A != None and B != None:
        return A.intersection(B)
    else:
        return {}



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
    AdjList = people.map(makepair)
    DriverAdj = dict(AdjList.collect())
    WorkerAdj = sc.broadcast(DriverAdj)
    Edges = AdjList.flatMapValues(lambda x: x)
    TriSet = Edges.map(lambda (k,v): ((k,v), 
             AintersectB(k,v,WorkerAdj.value)))
    Triangle = TriSet.flatMapValues(lambda x: x).map(lambda (k,v): 
             tuple(sorted([int(v),int(k[0]),int(k[1])],reverse=True)))
    output = set(Triangle.collect())
    #################  NO EDITS HERE  ###################
    print "\n\n*****************************************"
    print "\nTotal algorithm time: %.4f sec \n" % (time()-start)
    print "*****************************************\n\n""" 
    ###############  END NO EDITS HERE ################
    with open("triangles.out", "wb") as f:
        for friends in output:
            f.write(str(friends[0])+" "+str(friends[1])+" "+str(friends[2])+"\n") # TODO: Loop with f to write your result to file serially
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