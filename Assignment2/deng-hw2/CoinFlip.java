////////////////////////////////////////////////////////////////////////////////
// Detian Deng
////////////////////////////////////////////////////////////////////////////////
import java.util.Random;

class CoinFlip implements Runnable
{
  int thread_id;    // Variable containing specific id of this thread.
  long eachthread;
  long lastthread;
  int num_thread;
  Random toss = new Random();
  
  // Create some variables for testing.
  static long heads = 0L;

  public long getHeads()
  {
  	 return this.heads;
  } 

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
      this.heads = heads_inthread;
      //synchronized(CoinFlip.class){this.heads += heads_inthread;}
    }
    else
    {
      for ( int i=0; i<lastthread; i++ )
      {
        heads_inthread += this.toss.nextInt(2);
      }
      this.heads = heads_inthread;
      //synchronized(CoinFlip.class){this.heads += heads_inthread;}
    }
  }

  // Constructor: set thread id
  CoinFlip ( int id, long n_sample, int n_thread) 
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
      System.out.println ("Usage: CoinFlip #threads #iterations");
      return;
    } 

    long starttime = System.currentTimeMillis();
    // Get the number of threads we are going to run from the command line
    int numthreads = Integer.parseInt ( args[0] );
    long numIter = Long.parseLong ( args[1] );
    long totalHeads = 0;

    // Array to hold references to thread objects
    Thread[] threads = new Thread[numthreads];
    CoinFlip[] coinflips = new CoinFlip[numthreads];

    // create and start specified thread objects of class SynchronizedWorks
    for ( int i=0; i<numthreads; i++ )
    {
      coinflips[i] = new CoinFlip(i, numIter, numthreads);
      threads[i] = new Thread (coinflips[i]);
      threads[i].start();
    }

    // Await the completion of all threads
    for ( int i=0; i<numthreads; i++ )
    {
      try
      {
        threads[i].join();
        totalHeads += coinflips[i].getHeads();
      }
      catch (InterruptedException e)
      {
         System.out.println("Thread interrupted.  Exception: " + e.toString() +
                           " Message: " + e.getMessage()) ;
        return;
      }
    }

    long runtime = System.currentTimeMillis() - starttime;

    System.out.println(totalHeads + " heads in " + numIter + " coin tosses.");
    System.out.println("Elapsed time: " + runtime + "ms");
  }
}
