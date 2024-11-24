import matplotlib.pyplot as plt
import numpy as np

# 读取距离文件
def read_distances(file_path):
    with open(file_path, "r") as file:
        distances = [float(line.strip()) for line in file]
    return distances

# 绘制直方图
def plot_histogram(data, title, xlabel, output_file):
    plt.figure(figsize=(10, 6))
    plt.hist(data, bins=50, edgecolor='black', alpha=0.7)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel("Frequency")
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print("Histogram plot saved to {}".format(output_file))

# 处理单个文件集(普通和环绕几何）
def process_results(file_prefix):
    # 文件路径
    nearest_standard_file = f'../{file_prefix}_nearest_standard.txt'
    furthest_standard_file = f"../{file_prefix}_furthest_standard.txt"
    nearest_wraparound_file = f"../{file_prefix}_nearest_wraparound.txt"
    furthest_wraparound_file = f"../{file_prefix}_furthest_wraparound.txt"

    # 读取数据
    nearest_standard = read_distances(nearest_standard_file)
    furthest_standard = read_distances(furthest_standard_file)
    nearest_wraparound = read_distances(nearest_wraparound_file)
    furthest_wraparound = read_distances(furthest_wraparound_file)

    # 绘制直方图
    plot_histogram(nearest_standard, "Nearest Distance (Standard Geometry)", "Distance",
                   f"{file_prefix}_nearest_standard_hist.png")
    plot_histogram(furthest_standard, "Furthest Distance (Standard Geometry)", "Distance",
                   f"{file_prefix}_furthest_standard_hist.png")
    plot_histogram(nearest_wraparound, "Nearest Distance (Wraparound Geometry)", "Distance",
                   f"{file_prefix}_nearest_wraparound_hist.png")
    plot_histogram(furthest_wraparound, "Furthest Distance (Wraparound Geometry)", "Distance",
                   f"{file_prefix}_furthest_wraparound_hist.png")


# 主程序
if __name__ == "__main__":
    # 处理文件 1: 100000 locations
    process_results("output_100000")

    # 处理文件 2: 200000 locations
    process_results("output_200000")



