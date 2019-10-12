#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <thread>
#include <atomic>
#include <vector>
#include "core/utils.h"
#include "core/graph.h"

std::atomic <long> triangle_count;
std::vector<uintV> vertexVec;

uintV countTriangles(uintV *array1, uintE len1, uintV *array2, uintE len2, uintV u, uintV v)
{

    uintE i = 0, j = 0; // indexes for array1 and array2
    uintV count = 0;
    while ((i < len1) && (j < len2))
    {
        if (array1[i] == array2[j])
        {
            if ((array1[i] != u) && (array1[j] != v))
            {
                count++;
            }
            else
            {
                // triangle with self-referential edge -> ignore
            }
            i++;
            j++;
        }
        else if (array1[i] < array2[j])
        {
            i++;
        }
        else
        {
            j++;
        }
    }
    return count;
}

void printThreadStatistics(uint n_workers, std::thread *t, long *triangleCountEachThread,
  std::vector<double> &thread_time_taken, std::vector<long> &vertexCount, std::vector<long> &edgeCount){

  std::cout << "thread_id, num_vertices, num_edges, triangle_count, time_taken\n";
  for(uint i = 0; i < n_workers; i++){
    t[i].join();
    std::cout <<i<< ", " << vertexCount[i] << ", "<< edgeCount[i] << ", "<<triangleCountEachThread[i] << ", "
    << std::setprecision(TIME_PRECISION) << thread_time_taken[i]<<"\n";
  }
}

void printStatistics(double partitioning_time, double total_time_taken){
  // Print the overall statistics
  std::cout << "Number of triangles : " << triangle_count << "\n";
  std::cout << "Number of unique triangles : " << triangle_count / 3 << "\n";
  std::cout << "Partitioning time (in seconds) : " << std::setprecision(TIME_PRECISION) << partitioning_time << "\n";
  std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION) << total_time_taken << "\n";
}
void processEdge(Graph *g, uintE start, uintE end, uint tid, double* thread_time,
  long* triangleCountEachThread, std::vector<std::pair<uintV,uintV>> edgesVector){

  timer thread_timer;
  thread_timer.start();
  long local_triangle_count = 0;
  long count = 0;
  for (uintE i = start; i < end; i++){
    uintE u = edgesVector[i].first;
    uintE v = edgesVector[i].second;
    local_triangle_count += countTriangles(g->vertices_[u].getInNeighbors(),
                                           g->vertices_[u].getInDegree(),
                                           g->vertices_[v].getOutNeighbors(),
                                           g->vertices_[v].getOutDegree(),
                                           u,
                                           v);
  }
  *triangleCountEachThread = local_triangle_count;
  triangle_count.fetch_add(local_triangle_count);
  *thread_time = thread_timer.stop();

}

void getEdges(Graph *g, std::vector<std::pair<uintV,uintV>> &edgesVector, uintV n){

  for(int u = 0; u < n; u++){
  uintE out_degree = g->vertices_[u].getOutDegree();
    for (uintE i = 0; i < out_degree; i++){
      uintV v = g->vertices_[u].getOutNeighbor(i);
      edgesVector.push_back(std::make_pair(u,v));
    }
  }
}
void processVertex(Graph *g, uintV start, uintV end, uint tid, double* thread_time,
  long* triangleCountEachThread, long* edgeCountThread, long* vertexCountThread){
  timer thread_timer;
  thread_timer.start();
  long edgeCount = 0;
  long vertexCount = 0;
  long local_triangle_count = 0;
  //std::cout << "edgeCountThread1: " << *edgeCountThread << "\n";
  for (uintV u = start; u < end; u++)
  {
      // For each outNeighbor v, find the intersection of inNeighbor(u) and outNeighbor(v)
      vertexCount++;
      uintE out_degree = g->vertices_[u].getOutDegree();
      for (uintE i = 0; i < out_degree; i++)
      {
        uintV v = g->vertices_[u].getOutNeighbor(i);
        edgeCount++;
        local_triangle_count += countTriangles(g->vertices_[u].getInNeighbors(),
                                              g->vertices_[u].getInDegree(),
                                              g->vertices_[v].getOutNeighbors(),
                                              g->vertices_[v].getOutDegree(),
                                              u,
                                              v);
      }
  }
  *triangleCountEachThread = local_triangle_count;
  *edgeCountThread = edgeCount;
  *vertexCountThread = vertexCount;
  triangle_count.fetch_add(local_triangle_count);
  *thread_time = thread_timer.stop();

}

void getStartEndVertices(uintV *start, uintV *end, uintV n, uint n_workers){

  for(uint i = 0; i < n_workers - 1; i++){
    start[i] = (n/n_workers)*i;
    end[i] = ((i+1)*(n/n_workers));
  }
  start[n_workers-1] = n/n_workers*(n_workers-1);
  end[n_workers-1] = n;

}

void triangleCountingParallel(Graph &g, uint n_workers, int strategy)
{
  uintV n = g.n_;
  double total_time_taken = 0.0;
  timer t1;
  timer partitionTimer;
  long triangleCountEachThread[n_workers];
  std::vector<double> thread_time_taken(n_workers,0);
  std::vector<long> vertexCount(n_workers,0);
  std::vector<long> edgeCount(n_workers, 0);
  double partitionTime = 0.0;
  uintV start[n_workers];
  uintV end[n_workers];
  std::thread t[n_workers];

  t1.start();
  if(strategy == 1){
    partitionTimer.start();
    getStartEndVertices(start, end, n, n_workers);
    partitionTime = partitionTimer.stop();
    for(uint i = 0; i < n_workers; i++){
      t[i] = std::thread(processVertex, &g, start[i], end[i], i, &thread_time_taken[i], &triangleCountEachThread[i], &edgeCount[i], &vertexCount[i]); // i is tid
    }
  }

  if(strategy == 2){
    std::vector<std::pair<uintV,uintV>> edgesVector;
    partitionTimer.start();
    getEdges(&g ,edgesVector, n);
    uintE totalEdges = edgesVector.size();
    getStartEndVertices(start, end, totalEdges, n_workers);
    partitionTime = partitionTimer.stop();
    for(uint i = 0; i < n_workers; i++){
      edgeCount[i] = end[i] - start[i];
      t[i] = std::thread(processEdge, &g, start[i], end[i], i, &thread_time_taken[i], &triangleCountEachThread[i], edgesVector);
    }
  }

  printThreadStatistics(n_workers, t, triangleCountEachThread, thread_time_taken, vertexCount, edgeCount);
  total_time_taken = t1.stop();
  printStatistics(partitionTime, total_time_taken);
}

uintV getNextVertexToBeProcessed(std::atomic <long> &sharedCurr, uintV n){

  if(sharedCurr >= n){
    return -1;
  }
  return sharedCurr.fetch_add(1);
}

void dynamicCounting(Graph *g, uint tid, uintV n, double* thread_time,
  long* triangleCountEachThread, long* vertexCountThread, long* edgeCountThread,
  std::atomic <long>* sharedCurr, double* partitionTime){

  timer thread_timer;
  thread_timer.start();
  timer partitionTimer;
  long local_triangle_count = 0;
  long edgeCount = 0;
  long vertexCount = 0;
  while(true){

    partitionTimer.start();
    uintV u = getNextVertexToBeProcessed(*sharedCurr, n);
    *partitionTime += partitionTimer.stop();
    if(u == -1)
      break;
    vertexCount ++;
    uintE out_degree = g->vertices_[u].getOutDegree();
    for (uintE i = 0; i < out_degree; i++)
    {
      uintV v = g->vertices_[u].getOutNeighbor(i);
      edgeCount++;
      local_triangle_count += countTriangles(g->vertices_[u].getInNeighbors(),
                                             g->vertices_[u].getInDegree(),
                                             g->vertices_[v].getOutNeighbors(),
                                             g->vertices_[v].getOutDegree(),
                                             u,
                                             v);
    }
  }
  *vertexCountThread = vertexCount;
  *edgeCountThread = edgeCount;
  *triangleCountEachThread = local_triangle_count;
  triangle_count.fetch_add(local_triangle_count);
  *thread_time = thread_timer.stop();
}

void parallelStrategy3(Graph &g, uint n_workers)
{
  uintV n = g.n_;
  double partitionTime = 0.0;
  double total_time_taken = 0.0;
  std::atomic <long> sharedCurr{0};
  timer t1;
  std::vector<double> thread_time_taken(n_workers,0);
  std::vector<long> vertexCount(n_workers,0);
  std::vector<long> edgeCount(n_workers, 0);
  long triangleCountEachThread[n_workers];

  t1.start();
  std::thread t[n_workers];
  for(uint i = 0; i < n_workers; i++){
    t[i] = std::thread(dynamicCounting, &g, i, n, &thread_time_taken[i], &triangleCountEachThread[i], &vertexCount[i], &edgeCount[i], &sharedCurr, &partitionTime); // i is tid
  }
  printThreadStatistics(n_workers, t, triangleCountEachThread, thread_time_taken, vertexCount, edgeCount);
  total_time_taken = t1.stop();
  printStatistics(partitionTime, total_time_taken);
}

void triangleCountSerial(Graph &g)
{
    uintV n = g.n_;
    long triangle_count = 0;
    double time_taken = 0.0;
    timer t1;
    t1.start();
    for (uintV u = 0; u < n; u++)
    {
        uintE out_degree = g.vertices_[u].getOutDegree();
        for (uintE i = 0; i < out_degree; i++)
        {
            uintV v = g.vertices_[u].getOutNeighbor(i);
            triangle_count += countTriangles(g.vertices_[u].getInNeighbors(),
                                             g.vertices_[u].getInDegree(),
                                             g.vertices_[v].getOutNeighbors(),
                                             g.vertices_[v].getOutDegree(),
                                             u,
                                             v);
        }
    }
    time_taken = t1.stop();
    std::cout << "Number of triangles : " << triangle_count << "\n";
    std::cout << "Number of unique triangles : " << triangle_count / 3 << "\n";
    std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION) << time_taken << "\n";
}


int main(int argc, char *argv[])
{
    cxxopts::Options options("triangle_counting_serial", "Count the number of triangles using serial and parallel execution");
    options.add_options("custom", {
                                      {"nWorkers", "Number of workers", cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_WORKERS)},
                                      {"strategy", "Strategy to be used", cxxopts::value<uint>()->default_value(DEFAULT_STRATEGY)},
                                      {"inputFile", "Input graph file path", cxxopts::value<std::string>()->default_value("/scratch/assignment1/input_graphs/roadNet-CA")},
                                  });

    auto cl_options = options.parse(argc, argv);
    uint n_workers = cl_options["nWorkers"].as<uint>();
    uint strategy = cl_options["strategy"].as<uint>();
    std::string input_file_path = cl_options["inputFile"].as<std::string>();
    std::cout << std::fixed;
    std::cout << "Number of workers : " << n_workers << "\n";
    std::cout << "Task decomposition strategy : " << strategy << "\n";

    Graph g;
    std::cout << "Reading graph\n";
    g.readGraphFromBinary<int>(input_file_path);
    std::cout << "Created graph\n";

    switch (strategy)
    {
    case 0:
        std::cout << "\nSerial\n";
        triangleCountSerial(g);
        break;
    case 1:
        std::cout << "\nVertex-based work partitioning\n";
        triangleCountingParallel(g, n_workers, 1);
        break;
    case 2:
        std::cout << "\nEdge-based work partitioning\n";
        triangleCountingParallel(g, n_workers, 2);
        break;
    case 3:
        std::cout << "\nDynamic work partitioning\n";
        parallelStrategy3(g, n_workers);
        break;
    default:
        break;
    }

    return 0;
}
