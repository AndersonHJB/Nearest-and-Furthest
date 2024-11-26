import numpy as np
import matplotlib.pyplot as plt

def read_distances(filename):
    """读取距离数据文件，返回距离列表"""
    with open(filename, 'r') as f:
        distances = [float(line.strip()) for line in f]
    return distances

def plot_distribution(distances, title, xlabel, ylabel, output_file):
    """绘制距离分布的直方图"""
    plt.figure(figsize=(10, 6))
    plt.hist(distances, bins=100, edgecolor='black', alpha=0.7)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.savefig(output_file)
    plt.close()

def process_and_plot(prefix, dataset_name):
    """处理指定前缀的文件，绘制并保存图表"""
    geometries = ['standard', 'wraparound']
    versions = ['serial', 'parallel']
    distance_types = ['nearest', 'furthest']

    for geometry in geometries:
        for version in versions:
            for distance_type in distance_types:
                filename = f"../{prefix}_{distance_type}_{geometry}_{version}.txt"
                distances = read_distances(filename)
                title = f"{dataset_name} - {geometry.capitalize()} Geometry ({version.capitalize()})\n{distance_type.capitalize()} Distance Distribution"
                xlabel = f"{distance_type.capitalize()} Distance"
                ylabel = "Frequency"
                output_file = f"{prefix}_{distance_type}_{geometry}_{version}.png"
                plot_distribution(distances, title, xlabel, ylabel, output_file)
                print(f"Generated plot: {output_file}")

if __name__ == "__main__":
    # 处理随机生成的数据集
    process_and_plot('output_random', 'Random Points (100,000)')

    # 处理CSV文件的数据集
    process_and_plot('output_csv', 'CSV Points (100,000)')
