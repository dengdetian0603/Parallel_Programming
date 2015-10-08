# 0. set environment variable
export AWS_SECRET_ACCESS_KEY=G/PSTczGcK0/89MWkj2/dbg1RC9jnD8M0uUnCjVZ

export AWS_ACCESS_KEY_ID=AKIAIYXVAOKVM27GMWJA

# 1. create spark EC2 cluster
/Applications/spark-1.3.0/ec2/spark-ec2 -k PP_2015 -i /Users/dengdetian0603/PP_2015.pem -s 2 launch PP_hw4_deng

# command help
/Applications/spark-1.3.0/ec2/spark-ec2 -h 

# 2. ssh into the cluster
/Applications/spark-1.3.0/ec2/spark-ec2 -k PP_2015 -i /Users/dengdetian0603/PP_2015.pem login PP_hw4_deng

# 3.1 put data into cluster
## on local machine
scp -i /Users/dengdetian0603/PP_2015.pem ~/Desktop/Assignment4/data* root@ec2-52-4-234-145.compute-1.amazonaws.com:~
scp -i /Users/dengdetian0603/PP_2015.pem ~/Desktop/Assignment4/triangle_count.py root@ec2-52-4-234-145.compute-1.amazonaws.com:~

# on spark-ec2 master
tar xvf data.tar.gz
ephemeral-hdfs/bin/hadoop fs -copyFromLocal ~/friends1000/ friends1000


# 3.2 submit spark job
spark/bin/spark-submit triangle_count.py friends1000 spark://ec2-52-4-234-145.compute-1.amazonaws.com:7077

# 3.3 collect result
scp -i /Users/dengdetian0603/PP_2015.pem root@ec2-52-4-234-145.compute-1.amazonaws.com:~/triangles.out Desktop/Assignment4/

# 4. terminate cluster
/Applications/spark-1.3.0/ec2/spark-ec2 -k PP_2015 -i /Users/dengdetian0603/PP_2015.pem destroy PP_hw4_deng






# test code locally
/Applications/spark-1.3.0/bin/spark-submit ~/Desktop/Assignment4/triangle_count2.py friends1000/

