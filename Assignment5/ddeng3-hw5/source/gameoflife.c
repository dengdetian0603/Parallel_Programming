// gameoflife.c
// Name: Detian Deng
// JHED: ddeng3

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "mpi.h"

#define ITERATIONS 64
#define GRID_WIDTH  256
#define DIM  16     // assume a square grid

int main ( int argc, char** argv ) {

    // Dynamic data
    int global_grid[ 256 ] = 
       {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
        1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };

    // MPI Standard variable
    int num_procs;
    int ID, j;
    int iters = 0;

    // Messaging variables
    MPI_Status stat;

    // MPI Setup
    if ( MPI_Init( &argc, &argv ) != MPI_SUCCESS )
    {
        printf ( "MPI_Init error\n" );
    }

    MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );
    MPI_Comm_rank ( MPI_COMM_WORLD, &ID );

    assert ( DIM % num_procs == 0 );

    // TODO:
    int* shard;
    int shardLen = GRID_WIDTH / num_procs;
    // data for individual process
    shard = malloc(shardLen * sizeof(int));     
    int i;
    for (i = 0; i < shardLen; i++)
    {        
        shard[i] = global_grid[ID * shardLen + i];
    }

    // Serial code with 1 processor 
    if (num_procs == 1)
    {
        int row, col;
        // index list
        int lastrc[16]={15,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
        int nextRowc[16]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,0};
        // updated shard
        int newShard[GRID_WIDTH];
        int neighborAlive;  

        for (iters = 0; iters < ITERATIONS; iters++) 
        {
            for ( row = 0; row < DIM; row++){
                for ( col = 0; col < DIM; col++){
                    neighborAlive = 0;
                    neighborAlive += shard[lastrc[row] * DIM + lastrc[col]];       // + upper left value
                    neighborAlive += shard[lastrc[row] * DIM + col];               // + upper middle
                    neighborAlive += shard[lastrc[row] * DIM + nextRowc[col]];       // + upper right
                    neighborAlive += shard[row * DIM + lastrc[col]];               // + left 
                    neighborAlive += shard[row * DIM + nextRowc[col]];               // + right 
                    neighborAlive += shard[nextRowc[row] * DIM + lastrc[col]];       // + lower left
                    neighborAlive += shard[nextRowc[row] * DIM + col];               // + lower middle
                    neighborAlive += shard[nextRowc[row] * DIM + nextRowc[col]];       // + lower right

                    if ((shard[row * DIM + col] == 0)&&(neighborAlive == 3))
                    {   
                        newShard[row * DIM + col] = 1;                    // revive
                    }else if ((shard[row * DIM + col] == 1)&&((neighborAlive == 2)||(neighborAlive == 3)))
                    {
                        newShard[row * DIM + col] = 1;                    // survive
                    }else
                    {
                        newShard[row * DIM + col] = 0;                    // dead
                    }                   
                }
            }
            // update shard at each iteration
            for(row = 0; row < GRID_WIDTH; row++)
            {
                shard[row] = newShard[row];                                 
            }
            //printf("Serial Iteration %d : \n",iters + 1);
        }
        // final update global grid
        for(j = 0; j < GRID_WIDTH; j++)
        {
            global_grid[j] = shard[j];                                 
        }
    }
    // MPI code with 2, 4, 8, 16 processes
    else if((num_procs % 2 ==0) && (num_procs <= 16))
    {
        // Nerighbor process ID
        int last = (ID == 0) ? (num_procs - 1) : (ID - 1);
        int next = (ID + 1) % num_procs;
        for (iters = 0; iters < ITERATIONS; iters++) 
        {
            // Messages from neighbor process
            int FromLast[DIM];
            int FromNext[DIM];
            // Messages to neighbor process
            int ToLast[DIM];
            int ToNext[DIM];
            // Nearest neighbor exchange in a ring topology: e.g. 0 <--> 1 <--> 2 <--> 3 <--> 0
            for(i = 0; i < DIM; i++){
                ToLast[i] = shard[i];                   // first row of shard
                ToNext[i] = shard[shardLen - DIM + i];  // last row of shard
            }
            
            // Blocking Message Passing
            // tag = 1: pass down -->      tag = 2: pass up <--
            if( ID % 2 == 0){
                MPI_Send (&ToNext, DIM, MPI_INT, next, 1, MPI_COMM_WORLD);
                MPI_Recv (&FromLast,DIM,MPI_INT,last, 1, MPI_COMM_WORLD, &stat);
                MPI_Send (&ToLast, DIM, MPI_INT, last, 2, MPI_COMM_WORLD);
                MPI_Recv (&FromNext, DIM, MPI_INT, next, 2, MPI_COMM_WORLD, &stat);
            }
            else{
                MPI_Recv (&FromLast,DIM,MPI_INT,last, 1, MPI_COMM_WORLD, &stat);
                MPI_Send (&ToNext, DIM, MPI_INT, next, 1, MPI_COMM_WORLD);
                MPI_Recv (&FromNext, DIM, MPI_INT, next, 2, MPI_COMM_WORLD, &stat);
                MPI_Send (&ToLast, DIM, MPI_INT, last, 2, MPI_COMM_WORLD);
            }

            // Extended shard: [FromLast shard FromNext]
            int eShardLen = shardLen + 2 * DIM;
            int eShard[eShardLen];
            for(i = 0; i < DIM; i++)
            {
                eShard[i] = FromLast[i];
            }
            for(i = DIM; i < shardLen + DIM; i++)
            {
                eShard[i] = shard[i - DIM];
            }
            for(i = shardLen + DIM; i < eShardLen; i++)
            {
                eShard[i] = FromNext[i - DIM - shardLen];
            }

            // Update shard based on eShard at each iteration
            int row, col , lastRow, lastCol, nextRow, nextCol;
            int nAlive;
            for (i = 0; i < shardLen; i++)
            {
                row = i/DIM + 1;
                col = i%DIM;
                lastRow = row - 1;
                nextRow = row + 1;
                lastCol = (col==0) ? (DIM-1) : (col-1) ;
                nextCol = (col==(DIM-1)) ? (0) : (col + 1);

                nAlive = 0;
                nAlive += eShard[lastRow * DIM + lastCol] + 
                          eShard[lastRow * DIM + col] +
                          eShard[lastRow * DIM + nextCol] +
                          eShard[row * DIM + lastCol] +
                          eShard[row * DIM + nextCol] +
                          eShard[nextRow * DIM + lastCol] +
                          eShard[nextRow * DIM + col] +
                          eShard[nextRow * DIM + nextCol];

                if((shard[i] == 0) && (nAlive == 3))
                {
                    shard[i] = 1;
                }
                else if((shard[i] == 1) && ((nAlive == 2)||(nAlive ==3 )))
                {
                    shard[i] = 1;
                }
                else
                {
                    shard[i] = 0;
                }
            }
            //printf("MPI Iteration %d: \n", iters + 1);
        }
        // Combine all shards together
        if (ID != 0)
        {
            // all send to master process
            MPI_Send(shard, shardLen, MPI_INT, 0, 3, MPI_COMM_WORLD);
        }
        else
        {
            // master process receive all messeges
            int Received[num_procs-1][shardLen];
            for( i = 1; i < num_procs; i++)
            {
                MPI_Recv(&Received[i-1], shardLen, MPI_INT, i, 3, MPI_COMM_WORLD, &stat);
            }
            // update first shard
            for (i = 0; i < shardLen; i++)
            {
                global_grid[i] = shard[i];
            }
            // update the rest shards
            for (i = shardLen; i < GRID_WIDTH ; i++)
            {
                global_grid[i] = Received[(i - shardLen)/shardLen][i % shardLen];
            }
        }
    }



    // Output the updated grid state
    // FIXME: Feel free to print more iterations when you debug but only submit with the 64th iteration printing
    if ( ID == 0 && iters % ITERATIONS == 0 ) {
        printf ( "\nIteration %d: final grid:\n", iters );
        for ( j = 0; j < GRID_WIDTH; j++ ) {
            if ( j % DIM == 0 ) {
                printf( "\n" );
            }
            printf ( "%d  ", global_grid[j] );
        }
        printf( "\n" );
    }

    MPI_Finalize(); // finalize so I can exit
}