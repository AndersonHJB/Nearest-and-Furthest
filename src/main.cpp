#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <limits>
#include <omp.h>
#include <random>
#include <chrono>
#include <string>
#include <algorithm>

// 定义点结构
struct Point {
    double x, y;
};

// 随机生成点
std::vector<Point> generateRandomPoints(size_t numPoints) {
    std::vector<Point> points;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    for (size_t i = 0; i < numPoints; ++i) {
        points.push_back({dis(gen), dis(gen)});
    }
    return points;
}

// 从 CSV 文件读取点
std::vector<Point> readPointsFromCSV(const std::string& filename) {
    std::vector<Point> points;
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string xStr, yStr;
        if (std::getline(ss, xStr, ',') && std::getline(ss, yStr)) {
            points.push_back({std::stod(xStr), std::stod(yStr)});
        }
    }
    return points;
}

// 计算普通距离
double standardDistance(const Point& a, const Point& b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return std::sqrt(dx * dx + dy * dy);
}

// 计算环绕距离
double wraparoundDistance(const Point& a, const Point& b) {
    double dx = std::abs(a.x - b.x);
    double dy = std::abs(a.y - b.y);
    dx = std::min(dx, 1.0 - dx);
    dy = std::min(dy, 1.0 - dy);
    return std::sqrt(dx * dx + dy * dy);
}

// 串行计算最近和最远距离
void computeDistancesSerial(const std::vector<Point>& points, bool useWraparound,
                            std::vector<double>& nearestDistances, std::vector<double>& furthestDistances,
                            double& avgNearest, double& avgFurthest, bool useNaiveAlgorithm) {
    size_t n = points.size();
    nearestDistances.resize(n, std::numeric_limits<double>::max());
    furthestDistances.resize(n, 0.0);
    avgNearest = 0.0;
    avgFurthest = 0.0;

    if (useNaiveAlgorithm) {
        // 朴素算法：距离计算两次
        for (size_t i = 0; i < n; ++i) {
            double nearest = std::numeric_limits<double>::max();
            double furthest = 0.0;
            for (size_t j = 0; j < n; ++j) {
                if (i == j) continue;
                double dist = useWraparound ? wraparoundDistance(points[i], points[j])
                                            : standardDistance(points[i], points[j]);
                if (dist < nearest) nearest = dist;
                if (dist > furthest) furthest = dist;
            }
            nearestDistances[i] = nearest;
            furthestDistances[i] = furthest;
            avgNearest += nearest;
            avgFurthest += furthest;
        }
        // 计算平均值
        avgNearest /= n;
        avgFurthest /= n;
    } else {
        // 优化算法：距离计算一次
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                double dist = useWraparound ? wraparoundDistance(points[i], points[j])
                                            : standardDistance(points[i], points[j]);
                if (dist < nearestDistances[i]) nearestDistances[i] = dist;
                if (dist > furthestDistances[i]) furthestDistances[i] = dist;
                if (dist < nearestDistances[j]) nearestDistances[j] = dist;
                if (dist > furthestDistances[j]) furthestDistances[j] = dist;
            }
            avgNearest += nearestDistances[i];
            avgFurthest += furthestDistances[i];
        }
        avgNearest /= n;
        avgFurthest /= n;
    }
}

// 并行计算最近和最远距离
void computeDistancesParallel(const std::vector<Point>& points, bool useWraparound,
                              std::vector<double>& nearestDistances, std::vector<double>& furthestDistances,
                              double& avgNearest, double& avgFurthest, bool useNaiveAlgorithm,
                              const std::string& scheduleType, int chunkSize) {
    size_t n = points.size();
    nearestDistances.resize(n, std::numeric_limits<double>::max());
    furthestDistances.resize(n, 0.0);

    // 设置OpenMP调度策略
    omp_sched_t ompSchedule;
    if (scheduleType == "static") {
        ompSchedule = omp_sched_static;
    } else if (scheduleType == "dynamic") {
        ompSchedule = omp_sched_dynamic;
    } else if (scheduleType == "guided") {
        ompSchedule = omp_sched_guided;
    } else {
        ompSchedule = omp_sched_auto;
    }
    omp_set_schedule(ompSchedule, chunkSize);

    avgNearest = 0.0;
    avgFurthest = 0.0;

    if (useNaiveAlgorithm) {
        // 朴素算法：距离计算两次
#pragma omp parallel for schedule(runtime) reduction(+:avgNearest,avgFurthest)
        for (size_t i = 0; i < n; ++i) {
            double nearest = std::numeric_limits<double>::max();
            double furthest = 0.0;
            for (size_t j = 0; j < n; ++j) {
                if (i == j) continue;
                double dist = useWraparound ? wraparoundDistance(points[i], points[j])
                                            : standardDistance(points[i], points[j]);
                if (dist < nearest) nearest = dist;
                if (dist > furthest) furthest = dist;
            }
            nearestDistances[i] = nearest;
            furthestDistances[i] = furthest;
            avgNearest += nearest;
            avgFurthest += furthest;
        }
        // **修正：计算平均值**
        avgNearest /= n;
        avgFurthest /= n;
    } else {
        // 优化算法：距离计算一次
        std::vector<double> localNearestDistances(n, std::numeric_limits<double>::max());
        std::vector<double> localFurthestDistances(n, 0.0);

#pragma omp parallel
        {
            std::vector<double> privateNearestDistances(n, std::numeric_limits<double>::max());
            std::vector<double> privateFurthestDistances(n, 0.0);

#pragma omp for schedule(runtime)
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = i + 1; j < n; ++j) {
                    double dist = useWraparound ? wraparoundDistance(points[i], points[j])
                                                : standardDistance(points[i], points[j]);
                    if (dist < privateNearestDistances[i]) privateNearestDistances[i] = dist;
                    if (dist > privateFurthestDistances[i]) privateFurthestDistances[i] = dist;
                    if (dist < privateNearestDistances[j]) privateNearestDistances[j] = dist;
                    if (dist > privateFurthestDistances[j]) privateFurthestDistances[j] = dist;
                }
            }

#pragma omp critical
            {
                for (size_t i = 0; i < n; ++i) {
                    if (privateNearestDistances[i] < nearestDistances[i])
                        nearestDistances[i] = privateNearestDistances[i];
                    if (privateFurthestDistances[i] > furthestDistances[i])
                        furthestDistances[i] = privateFurthestDistances[i];
                }
            }
        }

        // 计算平均值
        avgNearest = 0.0;
        avgFurthest = 0.0;
#pragma omp parallel for reduction(+:avgNearest,avgFurthest)
        for (size_t i = 0; i < n; ++i) {
            avgNearest += nearestDistances[i];
            avgFurthest += furthestDistances[i];
        }
        avgNearest /= n;
        avgFurthest /= n;
    }
}

// 输出结果到文件
void writeResults(const std::string& filename, const std::vector<double>& data) {
    std::ofstream file(filename);
    for (double value : data) {
        file << value << "\n";
    }
    file.close();
}

// 处理单个数据集
void processDataset(const std::string& description, const std::vector<Point>& points, const std::string& outputPrefix,
                    int numThreads, const std::string& scheduleType, int chunkSize, bool useNaiveAlgorithm) {
    std::cout << "Processing dataset: " << description << "\n";

    omp_set_num_threads(numThreads);

    std::vector<double> nearestDistances, furthestDistances;
    double avgNearest = 0.0, avgFurthest = 0.0;

    // 串行计算：朴素算法
    auto start = std::chrono::high_resolution_clock::now();
    computeDistancesSerial(points, false, nearestDistances, furthestDistances, avgNearest, avgFurthest, true);
    auto end = std::chrono::high_resolution_clock::now();
    double durationSerialNaive = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Serial Naive (standard geometry) completed in " << durationSerialNaive << " seconds\n";
    std::cout << "Average nearest distance: " << avgNearest << "\n";
    std::cout << "Average furthest distance: " << avgFurthest << "\n\n";

    // 并行计算
    start = std::chrono::high_resolution_clock::now();
    computeDistancesParallel(points, false, nearestDistances, furthestDistances, avgNearest, avgFurthest,
                             useNaiveAlgorithm, scheduleType, chunkSize);
    end = std::chrono::high_resolution_clock::now();
    double durationParallel = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Parallel (" << (useNaiveAlgorithm ? "Naive" : "Optimized") << ", standard geometry) completed in "
              << durationParallel << " seconds\n";
    std::cout << "Average nearest distance: " << avgNearest << "\n";
    std::cout << "Average furthest distance: " << avgFurthest << "\n\n";

    writeResults(outputPrefix + "_nearest_standard.txt", nearestDistances);
    writeResults(outputPrefix + "_furthest_standard.txt", furthestDistances);

    // 重复以上步骤，针对环绕几何
    // 串行计算：朴素算法
    start = std::chrono::high_resolution_clock::now();
    computeDistancesSerial(points, true, nearestDistances, furthestDistances, avgNearest, avgFurthest, true);
    end = std::chrono::high_resolution_clock::now();
    durationSerialNaive = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Serial Naive (wraparound geometry) completed in " << durationSerialNaive << " seconds\n";
    std::cout << "Average nearest distance: " << avgNearest << "\n";
    std::cout << "Average furthest distance: " << avgFurthest << "\n\n";

    // 并行计算
    start = std::chrono::high_resolution_clock::now();
    computeDistancesParallel(points, true, nearestDistances, furthestDistances, avgNearest, avgFurthest,
                             useNaiveAlgorithm, scheduleType, chunkSize);
    end = std::chrono::high_resolution_clock::now();
    durationParallel = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Parallel (" << (useNaiveAlgorithm ? "Naive" : "Optimized") << ", wraparound geometry) completed in "
              << durationParallel << " seconds\n";
    std::cout << "Average nearest distance: " << avgNearest << "\n";
    std::cout << "Average furthest distance: " << avgFurthest << "\n\n";

    writeResults(outputPrefix + "_nearest_wraparound.txt", nearestDistances);
    writeResults(outputPrefix + "_furthest_wraparound.txt", furthestDistances);
}

// 显示程序使用方法
void printUsage(const char* programName) {
    std::cout << "Usage: " << programName << " [options]\n"
              << "Options:\n"
              << "  --threads <num>         Number of threads (default: max available)\n"
              << "  --schedule <type>       OpenMP schedule type: static, dynamic, guided, auto (default: dynamic)\n"
              << "  --chunk_size <size>     Chunk size for scheduling (default: 1)\n"
              << "  --algorithm <type>      Algorithm type: naive, optimized (default: naive)\n"
              << "  --help                  Show this help message\n";
}

int main(int argc, char* argv[]) {
    int numThreads = omp_get_max_threads();
    std::string scheduleType = "dynamic";
    int chunkSize = 1;
    bool useNaiveAlgorithm = true;

    // 解析命令行参数
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--threads" && i + 1 < argc) {
            numThreads = std::stoi(argv[++i]);
        } else if (arg == "--schedule" && i + 1 < argc) {
            scheduleType = argv[++i];
        } else if (arg == "--chunk_size" && i + 1 < argc) {
            chunkSize = std::stoi(argv[++i]);
        } else if (arg == "--algorithm" && i + 1 < argc) {
            std::string algo = argv[++i];
            if (algo == "naive") {
                useNaiveAlgorithm = true;
            } else if (algo == "optimized") {
                useNaiveAlgorithm = false;
            } else {
                std::cerr << "Unknown algorithm type: " << algo << "\n";
                return 1;
            }
        } else if (arg == "--help") {
            printUsage(argv[0]);
            return 0;
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            printUsage(argv[0]);
            return 1;
        }
    }

    // 处理随机生成的数据集
    size_t numPoints = 100000;
    std::vector<Point> randomPoints = generateRandomPoints(numPoints);
    processDataset("Random points (100,000)", randomPoints, "output_random", numThreads, scheduleType, chunkSize,
                   useNaiveAlgorithm);

    // 处理 CSV 数据集
    std::vector<Point> csvPoints = readPointsFromCSV("data/100000 locations.csv");
    processDataset("CSV points (100,000)", csvPoints, "100000_output_csv", numThreads, scheduleType, chunkSize,
                   useNaiveAlgorithm);

    std::vector<Point> csvPoints2 = readPointsFromCSV("data/200000 locations.csv");
    processDataset("CSV points (200,000)", csvPoints2, "200000_output_csv", numThreads, scheduleType, chunkSize,
                   useNaiveAlgorithm);

    return 0;
}

