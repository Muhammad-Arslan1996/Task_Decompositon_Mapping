#include <iostream>
#include <iomanip>
#include <thread>
#include <stdlib.h>
#include "core/utils.h"
#include "core/graph.h"

#ifdef USE_INT
#define INIT_PAGE_RANK 100000
#define EPSILON 1000
#define PAGE_RANK(x) (15000 + (5 * x) / 6)
#define CHANGE_IN_PAGE_RANK(x, y) std::abs(x - y)
typedef int64_t PageRankType;
#else
#define INIT_PAGE_RANK 1.0
#define EPSILON 0.01
#define DAMPING 0.85
#define PAGE_RANK(x) (1 - DAMPING + DAMPING * x)
#define CHANGE_IN_PAGE_RANK(x, y) std::fabs(x - y)
typedef float PageRankType;
#endif

std::mutex getVertexLock;
void printThreadStatistics(uint n_workers, std::vector<long> vertexCount, std::vector<long> edgeCount, std::vector<double> b1,
  std::vector<double> b2, std::vector<double> nVTime, std::vector<double> thread_time_taken){

  std::cout << "thread_id, num_vertices, num_edges, barrier1_time, barrier2_time, getNextVertex_time, total_time\n";
  for(uint i = 0; i < n_workers; i++){
    std::cout <<i<< ", "<<vertexCount[i]<< ", "<<edgeCount[i]<< ", "<<b1[i]<< ", "<<b2[i]<< ", "<<nVTime[i]<< ", "
    << std::setprecision(TIME_PRECISION) << thread_time_taken[i]<<"\n";
  }
}

void printStatistics(double sum_of_page_ranks, double total_time_taken, double partitioning_time){

  std::cout << "Sum of page rank : " << sum_of_page_ranks << "\n";
  std::cout <<"Partitioning time (in seconds) : "<< std::setprecision(TIME_PRECISION) << partitioning_time<<"\n";
  std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION)<< total_time_taken << "\n";
}

void pageRankSerial(Graph &g, int max_iters)
{
    uintV n = g.n_;

    PageRankType *pr_curr = new PageRankType[n];
    PageRankType *pr_next = new PageRankType[n];

    for (uintV i = 0; i < n; i++)
    {
        pr_curr[i] = INIT_PAGE_RANK;
        pr_next[i] = 0.0;
    }

    // Push based pagerank
    timer t1;
    double time_taken = 0.0;
    // Create threads and distribute the work across T threads
    // -------------------------------------------------------------------
    t1.start();
    for (int iter = 0; iter < max_iters; iter++)
    {
        // for each vertex 'u', process all its outNeighbors 'v'
        for (uintV u = 0; u < n; u++)
        {
            uintE out_degree = g.vertices_[u].getOutDegree();
            for (uintE i = 0; i < out_degree; i++)
            {
                uintV v = g.vertices_[u].getOutNeighbor(i);
                pr_next[v] += (pr_curr[u] / out_degree);
            }
        }
        for (uintV v = 0; v < n; v++)
        {
            pr_next[v] = PAGE_RANK(pr_next[v]);

            // reset pr_curr for the next iteration
            pr_curr[v] = pr_next[v];
            pr_next[v] = 0.0;
        }
    }
    time_taken = t1.stop();
    // -------------------------------------------------------------------

    PageRankType sum_of_page_ranks = 0;
    for (uintV u = 0; u < n; u++)
    {
        sum_of_page_ranks += pr_curr[u];
    }
    std::cout << "Sum of page rank : " << sum_of_page_ranks << "\n";
    std::cout << "Time taken (in seconds) : " << time_taken << "\n";
    delete[] pr_curr;
    delete[] pr_next;
}

void pageRankProcess(Graph *g, uint iterStart, uint iterEnd, uint tid, uint n,
  double *thread_time, std::atomic<PageRankType>* pr_next, PageRankType *pr_curr, std::mutex* l1,
   CustomBarrier* my_barrier, double *b1, double *b2, long* edgeCountThread,
   long* vertexCountThread, int max_iters){

  timer thread_timer, b1T, b2T;
  long edgeCount = 0;
  long vertexCount = 0;
  PageRankType* temp_next = new PageRankType [n];
  thread_timer.start();
  for(int i = 0; i < max_iters; i++){
    for (uintV u = iterStart; u < iterEnd; u++)
    {
        uintE out_degree = g->vertices_[u].getOutDegree();
        edgeCount += out_degree;
        for (uintE i = 0; i < out_degree; i++)
        {
          uintV v = g->vertices_[u].getOutNeighbor(i);
          temp_next[v] = pr_next[v];
          while (!pr_next[v].compare_exchange_weak(temp_next[v], pr_next[v] + (pr_curr[u] / out_degree)));
        }
    }
    *edgeCountThread = edgeCount;
    b1T.start();
    my_barrier->wait();
    *b1 += b1T.stop();
    for (uintV v = iterStart; v < iterEnd; v++)
    {
        vertexCount++;
        temp_next[v] = pr_next[v];
        while (!pr_next[v].compare_exchange_weak(temp_next[v], PAGE_RANK(pr_next[v])));
        // reset pr_curr for the next iteration
        pr_curr[v] = pr_next[v];
        pr_next[v] = 0.0;
    }
    *vertexCountThread = vertexCount;
    *thread_time = thread_timer.stop();
    b2T.start();
    my_barrier->wait();
    *b2 += b2T.stop();
  }
}

void getStartEndVertices(uintV *start, uintV *end, uintV n, uint n_workers){

  for(uint i = 0; i < n_workers - 1; i++){
    start[i] = (n/n_workers)*i;
    end[i] = ((i+1)*(n/n_workers));
  }
  start[n_workers-1] = n/n_workers*(n_workers-1);
  end[n_workers-1] = n;

}

void getStarEndEdges(Graph *g, uintV *start, uintV *end,  uintV n, uintE m, uint n_workers){
  int j = 0;
  uintE limit = 0;
  uintE out_degree;
  start[j] = 0;
  for(int u = 0; u < n; u++){
    out_degree = g->vertices_[u].getOutDegree();
    limit += out_degree;
    if(limit >= (m/n_workers)*(j+1)){
      end[j] = u;
      j++;
      start[j] = u;
      }
    if(j == n_workers-1){
      break;
    }
  }
  end[j] = n;
}

void pageRankParallel(Graph &g, int max_iters, uint n_workers, CustomBarrier &my_barrier, uint strategy)
{
    uintV n = g.n_;
    uintE m = g.m_;
    PageRankType *pr_curr = new PageRankType[n];
    std::atomic<PageRankType>* pr_next = new std::atomic <PageRankType>[n];
    std::mutex* l1 = new std::mutex[n];
    std::vector<double> thread_time_taken(n_workers,0);
    std::vector<long> vertexCount(n_workers,0);
    std::vector<long> edgeCount(n_workers, 0);
    std::vector<double> barrier1_time(n_workers, 0);
    std::vector<double> barrier2_time(n_workers, 0);
    std::vector<double> getVertex_time(n_workers, 0);
    uintV start[n_workers];
    uintV end[n_workers];
    timer partitionTimer;
    timer t1;
    double total_time_taken = 0.0;
    double partitioning_time = 0.0;
    for (uintV i = 0; i < n; i++)
    {
        pr_curr[i] = INIT_PAGE_RANK;
        pr_next[i] = 0.0;
    }
    if(strategy == 1){
      partitionTimer.start();
      getStartEndVertices(start, end, n, n_workers);
      partitioning_time += partitionTimer.stop();
    }
    else if(strategy == 2){
      partitionTimer.start();
      getStarEndEdges(&g , start, end, n, m, n_workers);
      partitioning_time += partitionTimer.stop();
    }
    std::thread t[n_workers];
    t1.start();
    for(uint i = 0; i < n_workers; i++){
      t[i] = std::thread(pageRankProcess, &g, start[i], end[i], i, n, &thread_time_taken[i],
        std::ref(pr_next), std::ref(pr_curr), std::ref(l1), &my_barrier, &barrier1_time[i],
        &barrier2_time[i], &edgeCount[i], &vertexCount[i], max_iters);
      }

    for(uint i = 0; i < n_workers; i++){
      t[i].join();
    }
    total_time_taken = t1.stop();

    printThreadStatistics(n_workers, vertexCount, edgeCount, barrier1_time, barrier2_time, getVertex_time, thread_time_taken);
    PageRankType sum_of_page_ranks = 0;
    for (uintV u = 0; u < n; u++){
        sum_of_page_ranks += pr_curr[u];
    }
    printStatistics(sum_of_page_ranks, total_time_taken, partitioning_time);

    delete[] pr_curr;
    delete[] pr_next;
    delete[] l1;
}

uintV getNextVertexToBeProcessed(std::atomic <long> &sharedCurr, uintV n, uint granularity){

  if(sharedCurr >= n){
    return -1;
  }
  return sharedCurr.fetch_add(granularity);
}

void pageRankDynamic(Graph *g, uint tid, uintV n, double *thread_time, std::atomic<PageRankType>* pr_next,
  PageRankType *pr_curr, std::mutex* l1, CustomBarrier* my_barrier, double *b1,
  double *b2, long* edgeCountThread, long* vertexCountThread, double* getvertexTime,
  std::atomic <long>* getVertexShared1,std::atomic <long>* getVertexShared2, uint granularity){

     timer thread_timer, b1T, gvTimer, b2T;
     long edgeCount = 0;
     long vertexCount = 0;
     double getVertecTimer = 0.0;
     PageRankType* temp_next = new PageRankType [n];
     thread_timer.start();
     while(true){
       gvTimer.start();
       uintV vertices = getNextVertexToBeProcessed(*getVertexShared1, n, granularity);
       *getvertexTime += gvTimer.stop();
       if(vertices == -1)
         break;
       for(int u = vertices; u  < vertices+granularity && u < n; u++){
         uintE out_degree = g->vertices_[u].getOutDegree();
         for (uintE i = 0; i < out_degree; i++)
         {
           uintV v = g->vertices_[u].getOutNeighbor(i);
           edgeCount++;
           temp_next[v] = pr_next[v];
           while (!pr_next[v].compare_exchange_weak(temp_next[v], pr_next[v] + (pr_curr[u] / out_degree)));
         }
       }
     }
     *edgeCountThread += edgeCount;
     b1T.start();
     my_barrier->wait();
     *b1 += b1T.stop();

     while(true){
       gvTimer.start();
       uintV vertices = getNextVertexToBeProcessed(*getVertexShared2, n, granularity);
       *getvertexTime += gvTimer.stop();
       if(vertices == -1)
         break;
      for(int v = vertices; v < vertices+granularity && v < n; v++){
        vertexCount++;
        temp_next[v] = pr_next[v];
        while (!pr_next[v].compare_exchange_weak(temp_next[v], PAGE_RANK(pr_next[v])));
        // reset pr_curr for the next iteration
        pr_curr[v] = pr_next[v];
        pr_next[v] = 0.0;
      }
     }
     *vertexCountThread += vertexCount;
     *thread_time += thread_timer.stop();
     b2T.start();
     my_barrier->wait();
     *b2 += b2T.stop();
}

void pageRankParallelDynamic(Graph &g, int max_iters, uint n_workers, CustomBarrier &my_barrier, uint granularity){

  uintV n = g.n_;
  PageRankType *pr_curr = new PageRankType[n];
  std::atomic<PageRankType>* pr_next = new std::atomic <PageRankType>[n];
  std::mutex* l1 = new std::mutex[n];
  std::vector<double> thread_time_taken(n_workers,0);
  std::vector<long> vertexCount(n_workers,0);
  std::vector<long> edgeCount(n_workers, 0);
  std::vector<double> barrier1_time(n_workers, 0);
  std::vector<double> barrier2_time(n_workers, 0);
  std::vector<double> getVertex_time(n_workers, 0);
  timer t1;
  double total_time_taken = 0.0;
  double partitioning_time = 0.0;
  std::atomic <long> getVertexShared1{0};
  std::atomic <long> getVertexShared2{0};
  for (uintV i = 0; i < n; i++)
  {
      pr_curr[i] = INIT_PAGE_RANK;
      pr_next[i] = 0.0;
  }
  std::thread t[n_workers];
  t1.start();
  for (int iter = 0; iter < max_iters; iter++)
  {
    getVertexShared1 = 0;
    getVertexShared2 = 0;
    for(uint i = 0; i < n_workers; i++){
      t[i] = std::thread(pageRankDynamic, &g, i, n, &thread_time_taken[i], std::ref(pr_next),
      std::ref(pr_curr), std::ref(l1), &my_barrier, &barrier1_time[i], &barrier2_time[i],
      &edgeCount[i], &vertexCount[i], &getVertex_time[i], &getVertexShared1, &getVertexShared2, granularity);
    }

    for(uint i = 0; i < n_workers; i++){
      t[i].join();
    }
  }
  total_time_taken = t1.stop();

  printThreadStatistics(n_workers, vertexCount, edgeCount, barrier1_time, barrier2_time, getVertex_time, thread_time_taken);
  PageRankType sum_of_page_ranks = 0;
  for (uintV u = 0; u < n; u++){
      sum_of_page_ranks += pr_curr[u];
  }
  for(int i =0; i<n_workers; i++){
    partitioning_time += getVertex_time[i];
  }
  printStatistics(sum_of_page_ranks, total_time_taken, partitioning_time);

  delete[] pr_curr;
  delete[] pr_next;
  delete[] l1;
}

int main(int argc, char *argv[])
{
    cxxopts::Options options("page_rank_push", "Calculate page_rank using serial and parallel execution");
    options.add_options("", {
                                {"nWorkers", "Number of workers", cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_WORKERS)},
                                {"nIterations", "Maximum number of iterations", cxxopts::value<uint>()->default_value(DEFAULT_MAX_ITER)},
                                {"strategy", "Strategy to be used", cxxopts::value<uint>()->default_value(DEFAULT_STRATEGY)},
                                {"granularity", "Granularity to be used", cxxopts::value<uint>()->default_value(DEFAULT_GRANULARITY)},
                                {"inputFile", "Input graph file path", cxxopts::value<std::string>()->default_value("/scratch/assignment1/input_graphs/roadNet-CA")},
                            });

    auto cl_options = options.parse(argc, argv);
    uint n_workers = cl_options["nWorkers"].as<uint>();
    uint strategy = cl_options["strategy"].as<uint>();
    uint granularity = cl_options["granularity"].as<uint>();
    uint max_iterations = cl_options["nIterations"].as<uint>();
    std::string input_file_path = cl_options["inputFile"].as<std::string>();

#ifdef USE_INT
    std::cout << "Using INT\n";
#else
    std::cout << "Using FLOAT\n";
#endif
    std::cout << std::fixed;
    std::cout << "Number of workers : " << n_workers << "\n";
    std::cout << "Task decomposition strategy : " << strategy << "\n";
    std::cout << "Granularity : " << granularity << "\n";
    std::cout << "Iterations : " << max_iterations << "\n";

    Graph g;
    std::cout << "Reading graph\n";
    g.readGraphFromBinary<int>(input_file_path);
    CustomBarrier my_barrier(n_workers);
    std::cout << "Created graph\n";
    switch (strategy)
    {
    case 0:
        std::cout << "\nSerial\n";
        pageRankSerial(g, max_iterations);
        break;
    case 1:
        std::cout << "\nVertex-based work partitioning\n";
        pageRankParallel(g, max_iterations, n_workers, my_barrier, strategy);
        break;
    case 2:
        std::cout << "\nEdge-based work partitioning\n";
        pageRankParallel(g, max_iterations, n_workers, my_barrier, strategy);
        break;
    case 3:
        std::cout << "\nDynamic work partitioning\n";
        pageRankParallelDynamic(g, max_iterations, n_workers, my_barrier, granularity);
        break;
    default:
        break;
    }

    return 0;
}
