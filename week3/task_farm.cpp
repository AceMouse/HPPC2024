/*
  Assignment: Make an MPI task farm. A "task" is a randomly generated integer.
  To "execute" a task, the worker sleeps for the given number of milliseconds.
  The result of a task should be send back from the worker to the master. It
  contains the rank of the worker
*/

#include <iostream>
#include <random>
#include <chrono>
#include <thread>
#include <array>
#include <stack>

// To run an MPI program we always need to include the MPI headers
#include <mpi.h>

const int NTASKS= 100; //5000;  // number of tasks
const int RANDOM_SEED=1234;

void master (int nworker) {
    std::array<int, NTASKS> task, result;
    std::vector<int> workers(nworker+1, 0); // vector of curr. k assigned to worker rank = index


    // set up a random number generator
    std::random_device rd;
    //std::default_random_engine engine(rd());
    std::default_random_engine engine;
    engine.seed(RANDOM_SEED);
    // make a distribution of random integers in the interval [0:30]
    std::uniform_int_distribution<int> distribution(0, 1000);

    // initialize list of tasks
    for (int& t : task)   t = distribution(engine);  

    // initialize list of available workers
    std::stack<int> avail_ws;
    for (int w=1; w <= nworker; w++)  avail_ws.push(w);

    // loop through tasks and send out to available workers
    int w;
    int data;
    MPI_Status status;
        
    for (int tidx=0; tidx < NTASKS; tidx++) {
        
        if (!avail_ws.empty()) {
            w = avail_ws.top();   
            avail_ws.pop();
        }
        else {
            MPI_Recv(&data, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            w = data;    
            result[workers[w]] = data;
        }
        MPI_Send(&task[tidx], 1, MPI_INT, w, 0, MPI_COMM_WORLD);
        workers[w] = tidx;
    }

    // send termination signal to all workers
    for (int w_count=0; w_count < nworker; w_count++) {
        MPI_Recv(&data, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        w = data;
        result[workers[w]] = data;
        MPI_Send(&task[0], 1, MPI_INT, w, 1, MPI_COMM_WORLD);
    }
    
    // Print out a status on how many tasks were completed by each worker
    for (int worker=1; worker<=nworker; worker++) {
        int tasksdone = 0; int workdone = 0; 
        
        for (int itask=0; itask<NTASKS; itask++)
        if (result[itask]==worker) {
            tasksdone++;
            workdone += task[itask];
        }
        std::cout << "Master: Worker " << worker << " solved " << tasksdone << 
                    " tasks\n";    
    }
}

// call this function to complete the task. It sleeps for task milliseconds
void task_function(int task) {
    std::this_thread::sleep_for(std::chrono::milliseconds(task));
}

void worker (int rank) {
    
    int data = 0;
    MPI_Status status;
    while (true) {  
        MPI_Recv(&data, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == 1) break;
        task_function(data);
        MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
}

int main(int argc, char *argv[]) {
    int nrank, rank;

    MPI_Init(&argc, &argv);                // set up MPI
    MPI_Comm_size(MPI_COMM_WORLD, &nrank); // get the total number of ranks
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // get the rank of this process

    if (rank == 0)       // rank 0 is the master
        master(nrank-1); // there is nrank-1 worker processes
    else                 // ranks in [1:nrank] are workers
        worker(rank);

    MPI_Finalize();      // shutdown MPI
}
