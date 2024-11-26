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

// 计算最近和最远距离（优化避免重复计算）
void computeDistances(const std::vector<Point>& points, bool useWraparound,
                      std::vector<double>& nearestDistances, std::vector<double>& furthestDistances) {
    size_t n = points.size();
    nearestDistances.resize(n, std::numeric_limits<double>::max());
    furthestDistances.resize(n, 0.0);

    // 并行计算距离
#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            double dist = useWraparound ? wraparoundDistance(points[i], points[j])
                                        : standardDistance(points[i], points[j]);
            // 更新最近和最远距离
#pragma omp critical
            {
                nearestDistances[i] = std::min(nearestDistances[i], dist);
                furthestDistances[i] = std::max(furthestDistances[i], dist);
                nearestDistances[j] = std::min(nearestDistances[j], dist);
                furthestDistances[j] = std::max(furthestDistances[j], dist);
            }
        }
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

// 测试单个数据集
void processDataset(const std::string& description, const std::vector<Point>& points, const std::string& outputPrefix) {
    std::cout << "Processing dataset: " << description << "\n";

    std::vector<double> nearestDistances, furthestDistances;

    // 普通几何距离
    auto start = std::chrono::high_resolution_clock::now();
    computeDistances(points, false, nearestDistances, furthestDistances);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Standard geometry completed in "<< duration.count() << " seconds\n";

    writeResults(outputPrefix + "_nearest_standard.txt", nearestDistances);
    writeResults(outputPrefix + "_furthest_standard.txt", furthestDistances);

    // 环绕几何距离
    start = std::chrono::high_resolution_clock::now();
    computeDistances(points, true, nearestDistances, furthestDistances);
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "Wraparound geometry completed in "
              << duration.count() << " seconds\n";

    writeResults(outputPrefix + "_nearest_wraparound.txt", nearestDistances);
    writeResults(outputPrefix + "_furthest_wraparound.txt", furthestDistances);

    // 打印平均距离
    double nearestSum = 0.0, furthestSum = 0.0;
    for (double d : nearestDistances) nearestSum += d;
    for (double d : furthestDistances) furthestSum += d;

    std::cout << "Average nearest distance: " << (nearestSum / points.size()) << "\n";
    std::cout << "Average furthest distance: " << (furthestSum / points.size()) << "\n";
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