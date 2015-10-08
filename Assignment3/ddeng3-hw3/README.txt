To avoid redundant, this readme file only describes how /friends1000/* data were analyzed. 
And since the EMR step does not allow upper case letter in the bucket name, instead of using suggested bucket name  <yourJHED>-PP2015-<RandomString>, I used detian.hw3 as the bucket name.

# Compile java source code
1. use command: cd ./source/java
2. edit the Makefile so that it has the correct directory to hadoop
3. use command: make compile

# Move codes to S3 bucket detian.hw3
1. s3cmd put -r ./source/java/* s3://detian.hw3/java_code/
2. s3cmd put -r ./source/python/* s3://detian.hw3/python_code/
3. to clean up, use command: make clean

# Create AWS EMR cluster:
The EMR cluster was created with the following specifications:
AMI version: 3.0.4 Hadoop Amazon 2.2.0
Master: 1 m1.medium
Core: 4 c3.2xlarge
Task: 0

# Run Python Implementation by adding a streaming program step:
Mapper file: s3://detian.hw3/python_code/fof.mapper.py
Reducer file: s3://detian.hw3/python_code/fof.reducer.py
Input S3 location: s3://detian.hw3/friends1000/
Output S3 location: s3://detian.hw3/py.output/
Arguments: None
Action on failure: Continue

# Run Java Implementation by adding a custom JAR step:
JAR s3 location s3://detian.hw3/java_code/FindTriangle.jar
Arguments: FindTriangle s3://detian.hw3/friends1000/ s3://detian.hw3/jar.output/
Action on failure: Continue