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



int main(int argc, char *argv[]) {

    std::vector<Point> points = readPointsFromCSV("data/100000 locations.csv");
    // std::cout << points.size() << std::endl;
    // std::cout << points.data() << std::endl;
    std::cout << "Loaded points from CSV:" << std::endl;
    for (size_t i = 0; i < points.size(); ++i) {
        std::cout << "Point " << i + 1 << ": (" << points[i].x << ", " << points[i].y << ")" << std::endl;
    }
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
