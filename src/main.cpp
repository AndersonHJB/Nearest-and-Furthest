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
std::vector<Point> read_csv(const std::string& filename) {
    std::vector<Point> points;
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string x_str, y_str;
        if (std::getline(ss, x_str, ',') && std::getline(ss, y_str)) {
            points.push_back({std::stod(x_str), std::stod(y_str)});
        }
    }
    return points;
}

// 计算标准几何距离
double standard_distance(const Point& p1, const Point& p2) {
    return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
}

int main(int argc, char *argv[]) {

    std::vector<Point> points = read_csv("data/100000 locations.csv");
    // std::cout << points.size() << std::endl;
    // std::cout << points.data() << std::endl;
    // std::cout << "Loaded points from CSV:" << std::endl;
    // for (size_t i = 0; i < points.size(); ++i) {
    //     std::cout << "Point " << i + 1 << ": (" << points[i].x << ", " << points[i].y << ")" << std::endl;
    // }
}
