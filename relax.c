#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>

double **masterArray;
double **dummyArray;
double **testingArray;

double arraySize;
int numThreads;
double tolerance;
double cellsPerThread;

int *visitorRecord;
int loop;
double globalMaxDiff;

pthread_cond_t condition;
pthread_mutex_t conditionMutex;
pthread_mutex_t maxDiffMutex;

typedef struct
{
    int secs;
    int usecs;
} TIMEDIFF;

// Timer function that returns time in seconds and microseconds between the 2 time given values
TIMEDIFF *diffTime(struct timeval *start, struct timeval *end)
{
    TIMEDIFF *diff = (TIMEDIFF *)malloc(sizeof(TIMEDIFF));

    if (start->tv_sec == end->tv_sec)
    {
        diff->secs = 0;
        diff->usecs = end->tv_usec - start->tv_usec;
    }
    else
    {
        diff->usecs = 1000000 - start->tv_usec;
        diff->secs = end->tv_sec - (start->tv_sec + 1);
        diff->usecs += end->tv_usec;
        if (diff->usecs >= 1000000)
        {
            diff->usecs -= 1000000;
            diff->secs += 1;
        }
    }

    return diff;
}

// Barrier implemented using condition variables as I'm working on an ARM64 Mac
int barrierWait(int threadNum)
{
    pthread_mutex_lock(&conditionMutex);

    // Adds calling thread's id to array of threads that are at barrier, then checks if all threads are at barier
    visitorRecord[threadNum] = threadNum;
    int allVisited = 1;
    for (int i = 0; i < numThreads; i++)
        if (visitorRecord[i] != i)
            allVisited = 0;

    // If all threads are at barrier, reset array to -1s and allow all threads to continue, else current thread enters busy wait
    if (allVisited == 1)
    {
        for (int i = 0; i < numThreads; i++)
            visitorRecord[i] = -1;
        pthread_cond_broadcast(&condition);
    }
    else
        pthread_cond_wait(&condition, &conditionMutex);

    pthread_mutex_unlock(&conditionMutex);
    return allVisited;
}

// Worker thread function
void avgAndCompare(int *threadNum)
{
    barrierWait(*threadNum);

    // Assigns each thread a first and last value of the array to work on
    int startCell = round(*threadNum * cellsPerThread);
    int endCell = round((*threadNum + 1) * cellsPerThread) - 1;
    double maxDiff;

    while (loop == 1)
    {
        maxDiff = 0;

        // Stores average of adjacent values for every cell in a temporary array
        for (int i = startCell; i <= endCell; i++)
        {
            int x = 1 + i / ((int)arraySize - 2);
            int y = 1 + i % ((int)arraySize - 2);

            dummyArray[x][y] = (masterArray[x - 1][y] +
                                masterArray[x + 1][y] +
                                masterArray[x][y - 1] +
                                masterArray[x][y + 1]) /
                               4;
        }

        barrierWait(*threadNum);

        // Inputs just caluclated values into master array and records largest difference between any 2 cells' previous and current value
        for (int i = startCell; i <= endCell; i++)
        {
            int x = 1 + i / ((int)arraySize - 2);
            int y = 1 + i % ((int)arraySize - 2);

            if (fabs(masterArray[x][y] - dummyArray[x][y]) > maxDiff)
            {
                maxDiff = fabs(masterArray[x][y] - dummyArray[x][y]);
            }
            masterArray[x][y] = dummyArray[x][y];
        }

        pthread_mutex_lock(&maxDiffMutex);
        if (maxDiff > globalMaxDiff)
            globalMaxDiff = maxDiff;
        pthread_mutex_unlock(&maxDiffMutex);

        // If calling thread is last one at barrier, check exit conditions and reset global max diff for next iteration
        int x = barrierWait(*threadNum);
        if (x == 1)
        {
            pthread_mutex_lock(&maxDiffMutex);
            if (globalMaxDiff < tolerance)
                loop = 0;
            globalMaxDiff = 0;
            pthread_mutex_unlock(&maxDiffMutex);
        }

        barrierWait(*threadNum);
    }

    pthread_exit(NULL);
}

// Function that takes running parameters, creates threads to perform the calculations, times the process, and outputs the resulting array
double run(double as, int nt, double tol, int verbose)
{
    // Update global variables
    arraySize = as;
    numThreads = nt;
    tolerance = tol;
    cellsPerThread = (as - 2) * (as - 2) / nt;

    globalMaxDiff = 0;
    loop = 1;

    // Dynamically allocate arrays
    masterArray = (double **)malloc(arraySize * sizeof(double *));
    dummyArray = (double **)malloc(arraySize * sizeof(double *));
    for (int i = 0; i < arraySize; i++)
    {
        masterArray[i] = (double *)malloc(arraySize * sizeof(double));
        dummyArray[i] = (double *)malloc(arraySize * sizeof(double));
    }

    visitorRecord = (int *)malloc(numThreads * sizeof(int));
    for (int i = 0; i < numThreads; i++)
        visitorRecord[i] = -1;

    // Initialise master array with values from testing array
    for (int i = 0; i < arraySize; i++)
    {
        for (int j = 0; j < arraySize; j++)
        {
            masterArray[i][j] = testingArray[i][j];
        }
    }

    // Begin timer
    struct timeval myTVstart, myTVend;
    TIMEDIFF *difference;
    gettimeofday(&myTVstart, NULL);

    // Initialise all threads and atomics
    pthread_t threads[numThreads];
    pthread_cond_init(&condition, NULL);
    pthread_mutex_init(&conditionMutex, NULL);
    pthread_mutex_init(&maxDiffMutex, NULL);

    // Create and join all threads
    int temp_arg[numThreads];
    for (int i = 0; i < numThreads; i++)
    {
        temp_arg[i] = i;
        pthread_create(&threads[i], NULL, (void *(*)(void *))avgAndCompare, &temp_arg[i]);
    }
    for (int i = 0; i < numThreads; i++)
        pthread_join(threads[i], NULL);

    // End timer
    gettimeofday(&myTVend, NULL);
    difference = diffTime(&myTVstart, &myTVend);

    if (verbose == 1)
    {
        printf("\nResult:\n");
        // Output result to console
        for (int i = 0; i < arraySize; i++)
        {
            for (int j = 0; j < arraySize; j++)
            {
                printf("%0.3f ", masterArray[i][j]);
            }
            printf("\n");
        }
    }

    // Free memory used for arrays
    for (int i = 0; i < arraySize; i++)
    {
        free(masterArray[i]);
        free(dummyArray[i]);
    }
    free(masterArray);
    free(dummyArray);
    free(visitorRecord);

    return difference->secs + ((double)difference->usecs / 1000000);
    ;
}

// Sequential version of algorithm for purpose of comparison
double runSeq(double as, double tol, int verbose)
{
    double **masterArraySeq;
    double **dummyArraySeq;
    double maxDiff = 9999;

    masterArraySeq = (double **)malloc(as * sizeof(double *));
    dummyArraySeq = (double **)malloc(as * sizeof(double *));

    for (int i = 0; i < as; i++)
    {
        masterArraySeq[i] = (double *)malloc(as * sizeof(double));
        dummyArraySeq[i] = (double *)malloc(as * sizeof(double));
    }

    for (int i = 0; i < as; i++)
    {
        for (int j = 0; j < as; j++)
        {
            masterArraySeq[i][j] = testingArray[i][j];
        }
    }

    // Begin timer
    struct timeval myTVstart, myTVend;
    TIMEDIFF *difference;
    gettimeofday(&myTVstart, NULL);

    while (maxDiff > tol)
    {
        maxDiff = 0;

        // Averaging 4 adjacent numbers for each and storing the result in a temporary array
        for (int i = 1; i < as - 1; i++)
        {
            for (int j = 1; j < as - 1; j++)
            {
                dummyArraySeq[i][j] = (masterArraySeq[i - 1][j] +
                                       masterArraySeq[i + 1][j] +
                                       masterArraySeq[i][j - 1] +
                                       masterArraySeq[i][j + 1]) /
                                      4;
            }
        }

        // Updating main array once all averages have been calculated
        // Also checking differences between all pairs of elements to see if within tolerance
        for (int i = 1; i < as - 1; i++)
        {
            for (int j = 1; j < as - 1; j++)
            {
                if (fabs(masterArraySeq[i][j] - dummyArraySeq[i][j]) > maxDiff)
                {
                    maxDiff = fabs(masterArraySeq[i][j] - dummyArraySeq[i][j]);
                }
                masterArraySeq[i][j] = dummyArraySeq[i][j];
            }
        }
    }

    // End timer
    gettimeofday(&myTVend, NULL);
    difference = diffTime(&myTVstart, &myTVend);

    if (verbose == 1)
    {
        printf("\nResult:\n");
        // Output result to console
        for (int i = 0; i < as; i++)
        {
            for (int j = 0; j < as; j++)
            {
                printf("%0.3f ", masterArraySeq[i][j]);
            }
            printf("\n");
        }
    }

    // Free memory used for arrays
    for (int i = 0; i < as; i++)
    {
        free(masterArraySeq[i]);
        free(dummyArraySeq[i]);
    }
    free(masterArraySeq);
    free(dummyArraySeq);

    return difference->secs + ((double)difference->usecs / 1000000);
}

int main()
{
    srand(time(0));
    rand();
    double t;

    //// Testing parameters ////
    double arraySize = 20.0;
    double tolerance = 0.01;
    int minThreads = 1;
    int maxThreads = 20;
    int numTests = 3;
    int verbose = 0;
    ////////////////////////////

    // Create and initialise testing array that all tests will use
    testingArray = (double **)malloc(arraySize * sizeof(double *));
    if (verbose == 1)
        printf("Initial array:\n");

    for (int i = 0; i < arraySize; i++)
    {
        testingArray[i] = (double *)malloc(arraySize * sizeof(double));
        for (int j = 0; j < arraySize; j++)
        {
            testingArray[i][j] = (double)rand() / (double)RAND_MAX;
            if (verbose == 1)
                printf("%0.3f ", testingArray[i][j]);
        }
        if (verbose == 1)
            printf("\n");
    }
    if (verbose == 1)
        printf("\n");

    // Run sequential algorithm on testing array
    t = 0;
    for (int i = 0; i < numTests; i++)
        t += runSeq(arraySize, tolerance, verbose);
    printf("seq %f\n", t / numTests);

    // Run parallel algorithm on testing array with a specified number of threads
    for (int i = minThreads; i <= maxThreads; i++)
    {
        t = 0;
        for (int j = 0; j < numTests; j++)
            t += run(arraySize, i, tolerance, verbose);
        printf("%d %f\n", i, t / numTests);
    }

    for (int i = 0; i < arraySize; i++)
        free(testingArray[i]);
    free(testingArray);
}