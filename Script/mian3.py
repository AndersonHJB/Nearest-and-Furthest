import numpy as np
import matplotlib.pyplot as plt
import os

def read_distances(filename):
    """读取距离数据文件"""
    with open(filename, 'r') as f:
        distances = [float(line.strip()) for line in f]
    return np.array(distances)

def plot_distribution(distances, title, xlabel, ylabel, output_filename):
    """绘制距离分布的直方图"""
    plt.figure(figsize=(10, 6))
    plt.hist(distances, bins=100, density=True, alpha=0.7, color='blue')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.savefig(output_filename)
    plt.close()
    print(f"图像已保存为 {output_filename}")

def process_dataset(prefix, description):
    """处理单个数据集，绘制最近和最远距离的分布图"""
    datasets = [
        (f"../{prefix}_nearest_standard.txt", "标准几何 - 最近距离", f"{prefix}_nearest_standard.png"),
        (f"../{prefix}_furthest_standard.txt", "标准几何 - 最远距离", f"{prefix}_furthest_standard.png"),
        (f"../{prefix}_nearest_wraparound.txt", "环绕几何 - 最近距离", f"{prefix}_nearest_wraparound.png"),
        (f"../{prefix}_furthest_wraparound.txt", "环绕几何 - 最远距离", f"{prefix}_furthest_wraparound.png"),
    ]

    for filename, title_suffix, output_filename in datasets:
        if not os.path.exists(filename):
            print(f"文件 {filename} 不存在，跳过。")
            continue
        distances = read_distances(filename)
        title = f"{description} - {title_suffix}"
        plot_distribution(distances, title, "距离", "频率", output_filename)

def main():
    # 处理随机生成的数据集
    process_dataset("output_random", "随机点 (100,000)")

    # 处理CSV数据集 (100,000)
    process_dataset("100000_output_csv", "CSV点 (100,000)")

    # 处理CSV数据集 (200,000)
    process_dataset("200000_output_csv", "CSV点 (200,000)")

if __name__ == "__main__":
    main()
