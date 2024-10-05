#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <omp.h>
#include <limits>
#include <cmath>

/* Gives us high-resolution timers. */
#define _POSIX_C_SOURCE 199309L
#include <time.h>

/* OSX timer includes */
#ifdef __MACH__
  #include <mach/mach.h>
  #include <mach/mach_time.h>
#endif

/**
* @brief Return the number of seconds since an unspecified time (e.g., Unix
*        epoch). This is accomplished with a high-resolution monotonic timer,
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

/**
* @brief Output the seconds elapsed while clustering.
*
* @param seconds Seconds spent on k-medoids clustering, excluding IO.
*/
static void print_time(double const seconds)
{
  std::cout << "k-medoids clustering time: " << std::fixed << std::setprecision(4) << seconds << "s" << std::endl;
}

int main(int argc, char* argv[])
{
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <filename> <num_clusters> <num_threads>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    int N_c = std::stoi(argv[2]);
    int num_threads = std::stoi(argv[3]);

    omp_set_num_threads(num_threads);
    omp_set_nested(1);
    omp_set_max_active_levels(2);

    // Read data from file
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open " + filename);
    }

    // Read first line to get N and D
    int N, D;
    file >> N >> D;

    // Create and resize the 2D vector
    std::vector<std::vector<double>> data(N, std::vector<double>(D));

    // Read the rest of the file
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < D; ++j) {
            if (!(file >> data[i][j])) {
                throw std::runtime_error("Error reading data from file");
            }
        }
    }

    // k-medoids clustering
    const int n_iters = 20;
    std::vector<int> clusters(N, 0);
    std::vector<std::vector<double>> medoids(N_c, std::vector<double>(D));

    double start_time = monotonic_seconds();
    
    // Initialize medoids with the first N_c-1 points
    std::vector<int> medoid_indices(N_c);
    for (int i = 0; i < N_c; ++i) {
        medoid_indices[i] = i;
        medoids[i] = data[i];
    }

    // Pre-compute squared distances between all points (upper triangular only)
    std::vector<std::vector<double>> distances(N);
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < N; ++i) {
        distances[i].resize(N - i);
        for (int j = i; j < N; ++j) {
            double dist = 0.0;
            for (int k = 0; k < D; ++k) {
                double diff = data[i][k] - data[j][k];
                dist += diff * diff;
            }
            distances[i][j - i] = dist; // Store squared distance
        }
    }

    // Add this before the main iteration loop
    std::vector<std::vector<int>> cluster_points(N_c);

    for (int iter = 0; iter < n_iters; ++iter) {
        // Assign points to nearest medoids
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < N; ++i) {
            double min_dist = std::numeric_limits<double>::max();
            int nearest_medoid = 0;
            for (int j = 0; j < N_c; ++j) {
                int m = medoid_indices[j];
                double dist = (i <= m) ? distances[i][m - i] : distances[m][i - m];
                if (dist < min_dist) {
                    min_dist = dist;
                    nearest_medoid = j;
                }
            }
            clusters[i] = nearest_medoid;
        }

        // Clear and rebuild cluster_points
        #pragma omp parallel
        {
            std::vector<std::vector<int>> local_cluster_points(N_c);

            #pragma omp for schedule(static)
            for (int i = 0; i < N; ++i) {
                local_cluster_points[clusters[i]].push_back(i);
            }

            #pragma omp critical
            {
                for (int i = 0; i < N_c; ++i) {
                    cluster_points[i].insert(cluster_points[i].end(), 
                                             local_cluster_points[i].begin(), 
                                             local_cluster_points[i].end());
                }
            }
        }

        // Update medoids
        bool medoids_changed = false;

        #pragma omp parallel for schedule(dynamic) reduction(|:medoids_changed)
        for (int i = 0; i < N_c; ++i) {
            double min_dist = std::numeric_limits<double>::max();
            int best_medoid = medoid_indices[i];
            
            #pragma omp parallel for schedule(dynamic) reduction(min:min_dist) \
                                     shared(best_medoid)
            for (size_t j = 0; j < cluster_points[i].size(); ++j) {
                int point_j = cluster_points[i][j];
                double total_dist = 0.0;
                
                for (int k : cluster_points[i]) {
                    int min_idx = std::min(point_j, k);
                    int max_idx = std::max(point_j, k);
                    total_dist += distances[min_idx][max_idx - min_idx];
                }
                
                if (total_dist < min_dist) {
                    min_dist = total_dist;
                    #pragma omp critical
                    {
                        if (total_dist < min_dist) {
                            min_dist = total_dist;
                            best_medoid = point_j;
                        }
                    }
                }
            }
            
            if (best_medoid != medoid_indices[i]) {
                medoids[i] = data[best_medoid];
                medoid_indices[i] = best_medoid;
                medoids_changed = true;
            }
        }

        if (!medoids_changed) break;
    }

    double end_time = monotonic_seconds();

    // Print time out
    double elapsed_time = end_time - start_time;
    print_time(elapsed_time);

    // Write elapsed time to a file
    std::ostringstream time_filename;
    time_filename << "elapsed_time_" << N_c << "_" << num_threads << ".txt";
    std::ofstream time_file(time_filename.str());
    if (!time_file.is_open()) {
        throw std::runtime_error("Unable to open " + time_filename.str() + " for writing");
    }
    time_file << std::fixed << std::setprecision(4) << elapsed_time << " " << N_c << " " << num_threads << std::endl;
    time_file.close();

    // Write final cluster assignments to clusters_openmp_Nclusters_Nthreads.txt
    std::ostringstream clusters_filename;
    clusters_filename << "clusters_openmp_" << N_c << "_" << num_threads << ".txt";
    std::ofstream clusters_file(clusters_filename.str());
    if (!clusters_file.is_open()) {
        throw std::runtime_error("Unable to open " + clusters_filename.str() + " for writing");
    }
    for (int i = 0; i < N; ++i) {
        clusters_file << clusters[i] << std::endl;
    }
    clusters_file.close();

    // Write final medoids to medoids_openmp_Nclusters_Nthreads.txt
    std::ostringstream medoids_filename;
    medoids_filename << "medoids_openmp_" << N_c << "_" << num_threads << ".txt";
    std::ofstream medoids_file(medoids_filename.str());
    if (!medoids_file.is_open()) {
        throw std::runtime_error("Unable to open " + medoids_filename.str() + " for writing");
    }
    for (int i = 0; i < N_c; ++i) {
        for (int j = 0; j < D; ++j) {
            medoids_file << medoids[i][j] << " ";
        }
        medoids_file << std::endl;
    }
    medoids_file.close();

    // end
    return 0;
}