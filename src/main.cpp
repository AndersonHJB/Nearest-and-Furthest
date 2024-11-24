#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <limits>
#include <omp.h>

struct Point {
    double x, y;
};

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

// 计算最近和最远距离
void computeDistances(const std::vector<Point>& points, bool useWraparound,
                      std::vector<double>& nearestDistances, std::vector<double>& furthestDistances) {
    size_t n = points.size();
    nearestDistances.resize(n, std::numeric_limits<double>::max());
    furthestDistances.resize(n, 0.0);

#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (i != j) {
                double dist = useWraparound ? wraparoundDistance(points[i], points[j])
                                            : standardDistance(points[i], points[j]);
                nearestDistances[i] = std::min(nearestDistances[i], dist);
                furthestDistances[i] = std::max(furthestDistances[i], dist);
            }
        }
    }
}

// 处理单个文件
void processFile(const std::string& filename, const std::string& outputPrefix) {
    std::vector<Point> points = readPointsFromCSV(filename);
    std::cout << points.size() << std::endl;
    std::vector<double> nearestDistances, furthestDistances;

    // 普通几何
    computeDistances(points, false, nearestDistances, furthestDistances);

    // 输出结果
    std::ofstream nearestOut(outputPrefix + "_nearest_standard.txt");
    std::ofstream furthestOut(outputPrefix + "_furthest_standard.txt");
    for (size_t i = 0; i < points.size(); ++i) {
        nearestOut << nearestDistances[i] << "\n";
        furthestOut << furthestDistances[i] << "\n";
    }
    nearestOut.close();
    furthestOut.close();

    // 环绕几何
    computeDistances(points, true, nearestDistances, furthestDistances);

    nearestOut.open(outputPrefix + "_nearest_wraparound.txt");
    furthestOut.open(outputPrefix + "_furthest_wraparound.txt");
    for (size_t i = 0; i < points.size(); ++i) {
        nearestOut << nearestDistances[i] << "\n";
        furthestOut << furthestDistances[i] << "\n";
    }

    // 打印平均距离
    double nearestSum = 0.0, furthestSum = 0.0;
    for (double d : nearestDistances) nearestSum += d;
    for (double d : furthestDistances) furthestSum += d;

    std::cout << "File: " << filename << "\n";
    std::cout << "Average nearest distance: " << (nearestSum / points.size()) << "\n";
    std::cout << "Average furthest distance: " << (furthestSum / points.size()) << "\n";
}

int main(int argc, char *argv[]) {
    // 处理文件 1: 100000 locations.csv
    processFile("data/100000 locations.csv", "output_100000");
    processFile("data/200000 locations.csv", "output_200000");

    return 0;
    // std::vector<Point> points = readPointsFromCSV("data/100000 locations.csv");
    // // std::cout << points.size() << std::endl;
    // // std::cout << points.data() << std::endl;
    // std::cout << "Loaded points from CSV:" << std::endl;
    // for (size_t i = 0; i < points.size(); ++i) {
    //     std::cout << "Point " << i + 1 << ": (" << points[i].x << ", " << points[i].y << ")" << std::endl;
    // }
    // std::vector<Point> points = readPointsFromCSV("data/100000 locations.csv");
    // if (points.empty()) {
    //     std::cerr << "Error: No points were read from the file!" << std::endl;
    // } else {
    //     std::cout << "Successfully read " << points.size() << " points from the file." << std::endl;
    //     for (size_t i = 0; i < std::min(points.size(), size_t(10)); ++i) {
    //         std::cout << "Point " << i << ": (" << points[i].x << ", " << points[i].y << ")" << std::endl;
    //     }
    // }

}
