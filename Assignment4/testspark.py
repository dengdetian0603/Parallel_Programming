# Word Count example

text = sc.textFile("friends.simple/")
#print text.collect()


from operator import add

def tokenize(text):
    return text.split()


words = text.flatMap(tokenize)
print words

wc = words.map(lambda x: (x,1))
print wc.toDebugString() # to inspect the lineage : 

# wc   |  PythonRDD[3] at RDD at PythonRDD.scala:42 []
# words|  PythonRDD[2] at RDD at PythonRDD.scala:42 []
# text |  friends.simple/100 MapPartitionsRDD[1] at textFile at NativeMethodAccessorImpl.java:-2 []
# ???  |  friends.simple/100 HadoopRDD[0] at textFile at NativeMethodAccessorImpl.java:-2 []

counts = wc.reduceByKey(add)
counts.saveAsTextFile("simple.output/wc.txt")

# ------------------------------------------------------------------------------------------------ # 
# Map/Reduce version of Triangle Count 
import os
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
                allTriads.append(str(triangle[0])+" "+str(triangle[1])+" "+str(triangle[2]))
    return allTriads  



people = sc.textFile("friends.simple")
triad = people.flatMap(GetTriad)
triadCount = triad.map(lambda x: (x,1))
triadSum = triadCount.reduceByKey(add)

triangles = triadSum.filter(lambda x: x[1]>1).map(lambda x: x[0])
triangles.saveAsTextFile("test1")

# ------------------------------------------------------------------------------------------------- #
# Spark version without shuffle
def makepair(text):
    words = text.split()
    k = words[0]
    v = set(words[1:])
    return (k,v)


people = sc.textFile("friends.simple")
AdjList = people.map(makepair)
Edges = AdjList.flatMapValues(lambda x: x)
EdgePool = sc.broadcast(Edges.collect())
Nodes = AdjList.keys()
Triads = Edges.cartesian(Nodes).filter(lambda (k,v): (v!=k[0])&(v!=k[1]))
TriEdges = Triads.filter(lambda (k,v): ((v,k[0]) in EdgePool.value) & ((v,k[1]) in EdgePool.value))
Triangle = TriEdges.map(lambda (k,v): tuple(sorted([int(v),int(k[0]),int(k[1])],reverse=True)))
print set(Triangle.collect())



# 
people = sc.textFile("friends1000")
AdjList = people.map(makepair)
Edges = AdjList.flatMapValues(lambda x: x)
AdjbyS = Edges.join(AdjList)
AdjbyT = AdjbyS.map(lambda (k,v): (v[0],[k,v[1]]))
AdjbyTS = AdjbyT.join(AdjList).map(lambda (k,v): ((k,v[0][0]), v[1].intersection(v[0][1])))

Triangle = AdjbyTS.flatMapValues(lambda x: x).map(lambda (k,v): tuple(sorted([int(v),int(k[0]),int(k[1])],reverse=True)))
print set(Triangle.collect())

# ------------------------------------------------------------------------------------------------ #
def AintersectB(k1,k2,TheList):
    A = TheList.get(k1)
    B = TheList.get(k2)
    if A != None and B != None:
        return A.intersection(B)
    else:
        return {}

from time import time
start = time()
people = sc.textFile("friends1000")
AdjList = people.map(makepair)
DriverAdj = dict(AdjList.collect())
WorkerAdj = sc.broadcast(DriverAdj)
Edges = AdjList.flatMapValues(lambda x: x)
TriSet = Edges.map(lambda (k,v): ((k,v), 
          AintersectB(k,v,WorkerAdj.value)))
Triangle = TriSet.flatMapValues(lambda x: x).map(lambda (k,v): 
          tuple(sorted([int(v),int(k[0]),int(k[1])],reverse=True)))
print set(Triangle.collect())
print "\nTotal algorithm time: %.4f sec \n" % (time()-start)










