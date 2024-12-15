#define _POSIX_C_SOURCE 199309L

#include<stdio.h>
#include<pthread.h>
#include<stdlib.h>
#include<string.h>
#include<float.h>
#include<math.h>
#include<omp.h>

#include <iostream>
#include <vector>

/* OSX timer includes */
#ifdef __MACH__
#include <mach/mach.h>
#include <mach/mach_time.h>
#endif

/**
 * @brief Return the number of seconds since an unspecified tim
 e (e.g., Unix
 *        epoch). This is accomplished with a high-resolution m
 onotonic timer,
 *        suitable for performance timing.
 *
 * @return The number of seconds.
 */
static inline double monotonic_seconds()
{
#ifdef __MACH__
  /* OSX */
  static mach_timebase_info_data_t info;
  static double seconds_per_unit;
  if(seconds_per_unit == 0) {
    mach_timebase_info(&info);
    seconds_per_unit = (info.numer / info.denom) / 1e9;
  }
  return seconds_per_unit * mach_absolute_time();
#else
  /* Linux systems */
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec * 1e-9;
#endif 
}

#define BUFFER_SIZE 1000000
#define ITERATIONS 20

struct cluster_assg_args {
  double **row_pointer;
  int thread_id; 
  int rows;
};


double **cluster_medoids = NULL;
int *cluster_assignment_map = NULL;
int converged = 1;
int rows_test = 0;
int dims_test = 0;
int clusters = 0;
int p = 0;


static void print_time(double const seconds) {
  printf("k-medoids clustering time: %0.04fs\n", seconds);
}

inline double distance(double *v1, double *v2) {

  int j;
  double sq_dist = 0.0;
  for(j = 0; j < dims_test; j++) {
    double cost = v1[j] - v2[j];		
    sq_dist += cost * cost;
  }
  return sqrt(sq_dist);
}

inline void sum(double *v1, double *v2) {

  int j;
  for(j = 0; j < dims_test; j++) {
    v1[j] += v2[j];
  }
}

int main(int argc, char *argv[]) {

  if(argc <= 3) {
    printf("Please enter command line params in this order: filename number_of_clusters number_of_threads. \n");
    exit(0);
  }

  int i, j, k;
  const char *delim = " \0";
  const char *file = argv[1];
  FILE *input_file = fopen(file, "r");
  char buffer[BUFFER_SIZE]; 
  char *token;
  int rows = 0, dims = 0;
  clusters = strtol(argv[2], NULL, 10);
  p = strtol(argv[3], NULL, 10);

  double **mat = NULL;

  int dims_set = 0;
  int first_line = 1;
  if(input_file == NULL) {
    printf("File open failed!");
    return(1);
  } else {
    while(fgets(buffer, BUFFER_SIZE, input_file) != NULL) {
      token = strtok(buffer, delim);

      if(first_line) {

        rows_test = (int) strtol(token, NULL, 10);

        while(token != NULL) {
          dims_test = (int) strtol(token, NULL, 10);
          token = strtok(NULL, delim);
        }

        mat = (double **) malloc(rows_test * sizeof(double *));
        int row;
        for(row = 0; row < rows_test; row++) {
          mat[row] = (double *) malloc(dims_test * sizeof(double));
        }

        first_line = 0;
        continue;
      }

      while(token != NULL) {
        mat[rows][dims] = strtod(token, NULL);
        token = strtok(NULL, delim);
        dims++;
      }

      dims_set = 1;
      rows++;
      dims = 0;
    }

    if(ferror(input_file)) {
      perror("Error occured: ");
    }

    fclose(input_file);
  }

  //std::cout << "Done Reading.. " << std::endl; 

  double start = monotonic_seconds();

  //Final global medoids - allocation and init.
  cluster_medoids = (double **) malloc(clusters * sizeof(double *));

  for (k = 0; k < clusters; k++) {
    cluster_medoids[k] = (double *) malloc(dims_test * sizeof(double));
  }

  for (k = 0; k < clusters; k++) {
    for (j = 0; j < dims_test; j++) {
      // Using the first k datapoints as first clusters 
      cluster_medoids[k][j] = mat[k][j];
    }	
  }

  //Cluster assignment map.
  cluster_assignment_map = (int *) malloc(rows_test * sizeof(int));

  int *cluster_xadj, *cluster_xadj_cp, *cluster_indices;
  cluster_xadj = (int *) malloc((clusters + 1) * sizeof(int));
  cluster_xadj_cp = (int *) malloc((clusters + 1) * sizeof(int));
  cluster_indices = (int *) malloc(rows_test * sizeof(int));

  int *cluster_min_ids;
  double *cluster_min_dist;
  cluster_min_ids = (int *) malloc(clusters * sizeof(int));
  cluster_min_dist = (double *) malloc(clusters * sizeof(double));

  int iter_count = 0;
  int myid;
  omp_set_num_threads(p);

  while (iter_count < ITERATIONS) {

    std::cout << "Iter count: " << iter_count << std::endl; 

    double part1_time = monotonic_seconds(); 

    int chunk_size = rows_test/p + 1;
#pragma omp parallel default(none) private(k, j) \
    shared(chunk_size, rows_test, dims_test, mat, clusters, converged, cluster_medoids, cluster_assignment_map)
    {

      #pragma omp for schedule(dynamic, chunk_size)
      for(i = 0; i < rows_test; i++) {

        double min = DBL_MAX;
        int cluster_id = 0;
        for (k = 0; k < clusters; k++) {
          double dist = distance(mat[i], cluster_medoids[k]);
          if (dist < min) {
            min = dist;
            cluster_id = k;
          }
        } // End of inner iter over cluster centers.

        if (cluster_assignment_map[i] != cluster_id) {
          cluster_assignment_map[i] = cluster_id;
          converged = 0;
        }

      } //end of pragma for

    } // Ending parallel section here..

    part1_time = monotonic_seconds() - part1_time; 
    printf("part1_time: %0.4f s \n", part1_time);

    // Initializing CSR structure.
    for (int i=0; i<=clusters; i++) { 
      cluster_xadj[i] = 0;
      cluster_xadj_cp[i] = 0;
    }

    for (i=0; i<rows_test; i++) { 
      cluster_xadj_cp[cluster_assignment_map[i]]++;
    }

    for (int i=0; i<clusters; i++)
      cluster_xadj_cp[i+1] += cluster_xadj_cp[i];

    for (int i=clusters; i>0; i--)
      cluster_xadj_cp[i] = cluster_xadj_cp[i-1];
    cluster_xadj_cp[0] = 0;

    for (int i=0; i<=clusters; i++)
      cluster_xadj[i] = cluster_xadj_cp[i];

    int cluster, pos; 
    for (int i=0; i<rows_test; i++) {
      cluster = cluster_assignment_map[i]; 
      pos = cluster_xadj_cp[cluster]++;
      cluster_indices[pos] = i;
    }

    for (int i=0; i<clusters; i++)
      cluster_min_dist[i] = DBL_MAX;

    double part2_time;
    part2_time = monotonic_seconds();

    // Restarting a new parallel section here..
    int thread_id;
    int cluster_chunk_size = clusters/p + 1;
#pragma omp parallel default(none) private(thread_id, k, j) \
    shared(chunk_size, cluster_chunk_size, rows_test, dims_test, mat, clusters, converged, cluster_medoids, cluster_assignment_map, cluster_xadj, cluster_indices, cluster_min_dist, cluster_min_ids)
    {

      #pragma omp for schedule(dynamic, 100)
      for (i=0; i<rows_test; i++) { // Checking if i should be the new medoid of its cluster

        thread_id = omp_get_thread_num();

        // Perform the all pairs comparison
        // between each data point in the cluster i
        int cluster_id = cluster_assignment_map[i];

        int idx_j;
        double avg_dist=0;
        for (j=cluster_xadj[cluster_id]; j<cluster_xadj[cluster_id+1]; j++) {
          idx_j = cluster_indices[j];
          double dist = distance(mat[i], mat[idx_j]);
          avg_dist += dist;
        }
        avg_dist = avg_dist / (cluster_xadj[cluster_id+1] - cluster_xadj[cluster_id]);

        // serial access -- critical section
        #pragma omp critical
        {
          // Code that should only be executed by one thread at a time
          // Each data point needs to check if it 
            // resulted in the min. avg. dist.
          if (avg_dist < cluster_min_dist[cluster_id]) {
            cluster_min_dist[cluster_id] = avg_dist;
            cluster_min_ids[cluster_id] = i; 
          }
        }

      } // Done iterating over each cluster

    } //end of parallel section

    part2_time = monotonic_seconds() - part2_time; 
    printf("part2_time: %0.4f s \n", part2_time);

    // Update the new medoid dimensions. 
    int min_id; 
    for(i=0; i<clusters; i++) { 
      min_id = cluster_min_ids[i];
      for (j=0; j<dims_test; j++)
        cluster_medoids[i][j] = mat[min_id][j];
    }

    iter_count++;
/*
    if (converged && iter_count)
      break; 

    converged = 1;
*/
  } //end while - Indent required.
  printf("it stop at %d\n",iter_count);
  double time = monotonic_seconds() - start;
  print_time(time);

  FILE *f = fopen("clusters.txt", "w");
  if (f == NULL) {
    printf("Error opening file!\n");
    exit(1);
  }

  for(i = 0; i < rows_test; i++) {
    fprintf(f, "%d\n", cluster_assignment_map[i]);
  }

  fclose(f);

  f = fopen("medoids.txt", "w");
  if (f == NULL) {
    printf("Error opening file!\n");
    exit(1);
  }

  for(i = 0; i < clusters; i++) {
    for(j = 0; j < dims_test; j++) {
      fprintf(f, "%e ", cluster_medoids[i][j]);
    }
    fprintf(f, "\n");
  }

  fclose(f);

  //REMEMBER TO FREE POINTERS - THERE IS NO GARBAGE COLLECTION! 
  for(k = 0; k < clusters; k++) {
    free(cluster_medoids[k]);
  }
  free(cluster_medoids);
  free(cluster_assignment_map);

  return 0;
}
