////////////////////////////////////////////////////////////////////////////////
// Detian Deng
////////////////////////////////////////////////////////////////////////////////
import java.util.Random;

class CoinStartUp implements Runnable
{
  int thread_id;    // Variable containing specific id of this thread.
  long eachthread;
  long lastthread;
  int num_thread;
  Random toss = new Random();
  
  // Create some variables for testing.
  static long heads = 0L;

 
  // Run: overides Runnabale.Run, thread entry point
  public void run ()
  {
    long heads_inthread = 0;
    if (this.thread_id < this.num_thread)
    {
      for ( int i=0; i<this.eachthread; i++ )
      {
        heads_inthread += this.toss.nextInt(2);
      }
      synchronized(CoinStartUp.class){this.heads += heads_inthread;}
    }
    else
    {
      for ( int i=0; i<lastthread; i++ )
      {
        heads_inthread += this.toss.nextInt(2);
      }
      synchronized(CoinStartUp.class){this.heads += heads_inthread;}
    }
  }

  // Constructor: set thread id
  CoinStartUp ( int id, long n_sample, int n_thread) 
  {
    this.thread_id = id;
    this.eachthread = n_sample/n_thread;
    this.lastthread = n_sample - (n_thread - 1)*this.eachthread;
    this.num_thread = n_thread;
  }

  public static void main ( String[] args )
  {
    if ( 2 != args.length ) 
    {
      System.out.println ("Usage: CoinStartUp #threads #iterations");
      return;
    } 

    
    // Get the number of threads we are going to run from the command line
    int numthreads = Integer.parseInt ( args[0] );
    long numIter = Long.parseLong ( args[1] );

    long starttime = System.nanoTime();
    // Array to hold references to thread objects
    Thread[] threads = new Thread[numthreads];

    // create and start specified thread objects of class SynchronizedWorks
    for ( int i=0; i<numthreads; i++ )
    {
      threads[i] = new Thread ( new CoinStartUp(i, numIter, numthreads) );
    }
    long headuptime = (System.nanoTime() - starttime)/1000;

    System.out.println("Startup time: " + headuptime + " microseconds with " + numthreads + " threads.");
  }
}
