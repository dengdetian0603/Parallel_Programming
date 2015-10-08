import java.io.IOException;
import java.util.*;

import org.apache.hadoop.fs.Path;
import org.apache.hadoop.conf.*;
import org.apache.hadoop.io.*;
import org.apache.hadoop.mapred.*;
import org.apache.hadoop.util.*;

public class FindTriangle {
public static class Map extends MapReduceBase implements Mapper<LongWritable, Text, Text, IntWritable> 
{
    
    private final static IntWritable one = new IntWritable(1); // output value

    public void map(LongWritable key, Text value, OutputCollector<Text, IntWritable> output, Reporter reporter) throws IOException 
    {
        Text outputKey = new Text(); // define key type

        String[] person = value.toString().split("\\s"); // keyperson friend1 friend2 ... friendn
        int[] friends = new int[3];
        int n = person.length;
//TODO: consider n > 2 only? 
        for (int j = 1; j < n; j++) 
        {
            for (int k = j + 1; k < n; k++) 
            {
                friends[0] = Integer.parseInt(person[0]);
                friends[1] = Integer.parseInt(person[j]);
                friends[2] = Integer.parseInt(person[k]);
                Arrays.sort(friends); // create a unique triangle candidates: keyperson friendj friendk

                outputKey.set(friends[0] + " " + friends[1] + " " + friends[2]);
                output.collect(outputKey, one);
            }
        }
    }
}

public static class Reduce extends MapReduceBase implements Reducer<Text, IntWritable, Text, Text> 
{
    public void reduce(Text key, Iterator<IntWritable> values, OutputCollector<Text, Text> output, Reporter reporter) throws IOException 
    {
        int sum = 0;
        while (values.hasNext())
        {
            sum += values.next().get(); // real triangles should appear more than once
        }

        if (sum >1 )
        {
            Text outputValue = new Text();
            Text key2 = new Text();
            Text key3 = new Text();

            String[] triangle = key.toString().split("\\s");
            key2.set(triangle[1] + " " + triangle[0] + " " + triangle[2]);
            key3.set(triangle[2] + " " + triangle[0] + " " + triangle[1]);
            output.collect(key, outputValue);
            output.collect(key2, outputValue);
            output.collect(key3, outputValue);
        }
    }
}

public static void main(String[] args) throws Exception 
{
//  code below follows the WordCount example without using Combiner
    JobConf conf = new JobConf(FindTriangle.class);
    conf.setJobName("Find Triangles");

    conf.setOutputKeyClass(Text.class);
    conf.setOutputValueClass(IntWritable.class);

    conf.setMapperClass(Map.class);
//  conf.setCombinerClass(Reduce.class);
    conf.setReducerClass(Reduce.class);

    conf.setInputFormat(TextInputFormat.class);
    conf.setOutputFormat(TextOutputFormat.class);

    FileInputFormat.setInputPaths(conf, new Path(args[0]));
    FileOutputFormat.setOutputPath(conf, new Path(args[1]));

    JobClient.runJob(conf);
}

}

