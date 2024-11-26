#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <limits>
#include <omp.h>
#include <random>
#include <chrono>

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
    return std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2));
}

// 计算环绕距离
double wraparoundDistance(const Point& a, const Point& b) {
    double dx = std::min(std::abs(a.x - b.x), 1.0 - std::abs(a.x - b.x));
    double dy = std::min(std::abs(a.y - b.y), 1.0 - std::abs(a.y - b.y));
    return std::sqrt(dx * dx + dy * dy);
}

// 串行计算最近和最远距离
void computeDistancesSerial(const std::vector<Point>& points, bool useWraparound,
                            std::vector<double>& nearestDistances, std::vector<double>& furthestDistances,
                            double& avgNearest, double& avgFurthest) {
    size_t n = points.size();
    nearestDistances.resize(n, std::numeric_limits<double>::max());
    furthestDistances.resize(n, 0.0);
    avgNearest = 0.0;
    avgFurthest = 0.0;

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            double dist = useWraparound ? wraparoundDistance(points[i], points[j])
                                        : standardDistance(points[i], points[j]);
            nearestDistances[i] = std::min(nearestDistances[i], dist);
            furthestDistances[i] = std::max(furthestDistances[i], dist);
            nearestDistances[j] = std::min(nearestDistances[j], dist);
            furthestDistances[j] = std::max(furthestDistances[j], dist);
        }
        avgNearest += nearestDistances[i];
        avgFurthest += furthestDistances[i];
    }

    // 计算平均值
    avgNearest /= n;
    avgFurthest /= n;
}

// 并行计算最近和最远距离
void computeDistancesParallel(const std::vector<Point>& points, bool useWraparound,
                              std::vector<double>& nearestDistances, std::vector<double>& furthestDistances,
                              double& avgNearest, double& avgFurthest) {
    size_t n = points.size();
    nearestDistances.resize(n, std::numeric_limits<double>::max());
    furthestDistances.resize(n, 0.0);
    std::vector<double> nearestSums(n, 0.0);
    std::vector<double> furthestSums(n, 0.0);

#pragma omp parallel
    {
        std::vector<double> localNearestDistances(n, std::numeric_limits<double>::max());
        std::vector<double> localFurthestDistances(n, 0.0);

#pragma omp for schedule(dynamic)
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                double dist = useWraparound ? wraparoundDistance(points[i], points[j])
                                            : standardDistance(points[i], points[j]);

                // 更新i和j的最近和最远距离
                if (dist < localNearestDistances[i]) localNearestDistances[i] = dist;
                if (dist > localFurthestDistances[i]) localFurthestDistances[i] = dist;
                if (dist < localNearestDistances[j]) localNearestDistances[j] = dist;
                if (dist > localFurthestDistances[j]) localFurthestDistances[j] = dist;
            }
        }

        // 合并结果
#pragma omp critical
        {
            for (size_t i = 0; i < n; ++i) {
                if (localNearestDistances[i] < nearestDistances[i]) nearestDistances[i] = localNearestDistances[i];
                if (localFurthestDistances[i] > furthestDistances[i]) furthestDistances[i] = localFurthestDistances[i];
            }
        }
    }

    // 计算平均值
    avgNearest = 0.0;
    avgFurthest = 0.0;
    for (size_t i = 0; i < n; ++i) {
        avgNearest += nearestDistances[i];
        avgFurthest += furthestDistances[i];
    }
    avgNearest /= n;
    avgFurthest /= n;
}

// 输出结果到文件
void writeResults(const std::string& filename, const std::vector<double>& data) {
    std::ofstream file(filename);
    for (double value : data) {
        file << value << "\n";
    }
    file.close();
}

// 测试单个数据集
void processDataset(const std::string& description, const std::vector<Point>& points, const std::string& outputPrefix) {
    std::cout << "Processing dataset: " << description << "\n";

    std::vector<double> nearestDistances, furthestDistances;
    double avgNearest = 0.0, avgFurthest = 0.0;

    // 串行：普通几何
    auto start = std::chrono::high_resolution_clock::now();
    computeDistancesSerial(points, false, nearestDistances, furthestDistances, avgNearest, avgFurthest);
    auto end = std::chrono::high_resolution_clock::now();
    double duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Serial (standard geometry) completed in " << duration << " seconds\n";
    std::cout << "Average nearest distance: " << avgNearest << "\n";
    std::cout << "Average furthest distance: " << avgFurthest << "\n\n";

    writeResults(outputPrefix + "_nearest_standard_serial.txt", nearestDistances);
    writeResults(outputPrefix + "_furthest_standard_serial.txt", furthestDistances);

    // 并行：普通几何
    start = std::chrono::high_resolution_clock::now();
    computeDistancesParallel(points, false, nearestDistances, furthestDistances, avgNearest, avgFurthest);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Parallel (standard geometry) completed in " << duration << " seconds\n";
    std::cout << "Average nearest distance: " << avgNearest << "\n";
    std::cout << "Average furthest distance: " << avgFurthest << "\n\n";

    writeResults(outputPrefix + "_nearest_standard_parallel.txt", nearestDistances);
    writeResults(outputPrefix + "_furthest_standard_parallel.txt", furthestDistances);

    // 串行：环绕几何
    start = std::chrono::high_resolution_clock::now();
    computeDistancesSerial(points, true, nearestDistances, furthestDistances, avgNearest, avgFurthest);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Serial (wraparound geometry) completed in " << duration << " seconds\n";
    std::cout << "Average nearest distance: " << avgNearest << "\n";
    std::cout << "Average furthest distance: " << avgFurthest << "\n\n";

    writeResults(outputPrefix + "_nearest_wraparound_serial.txt", nearestDistances);
    writeResults(outputPrefix + "_furthest_wraparound_serial.txt", furthestDistances);

    // 并行：环绕几何
    start = std::chrono::high_resolution_clock::now();
    computeDistancesParallel(points, true, nearestDistances, furthestDistances, avgNearest, avgFurthest);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Parallel (wraparound geometry) completed in " << duration << " seconds\n";
    std::cout << "Average nearest distance: " << avgNearest << "\n";
    std::cout << "Average furthest distance: " << avgFurthest << "\n\n";

    writeResults(outputPrefix + "_nearest_wraparound_parallel.txt", nearestDistances);
    writeResults(outputPrefix + "_furthest_wraparound_parallel.txt", furthestDistances);
}

int main() {
    // 处理随机生成的数据集
    size_t numPoints = 100000;
    std::vector<Point> randomPoints = generateRandomPoints(numPoints);
    processDataset("Random points (100,000)", randomPoints, "output_random");

    // 处理 CSV 数据集
    std::vector<Point> csvPoints = readPointsFromCSV("data/100000 locations.csv");
    processDataset("CSV points (100,000)", csvPoints, "output_csv");

    return 0;
}