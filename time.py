import argparse
import numpy as np
import os
import subprocess
import sys
import matplotlib.pyplot as plt

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run the analysis on given datasets.')
    parser.add_argument('--datasets', nargs='+', default=['HIV', 'SARS-CoV-2', 'Swine H1N1'],
                        help='Names of the datasets (e.g., HIV SARS-CoV-2 Swine H1N1)')
    parser.add_argument('--num_of_seed', type=int, default=10, help='Number of seeds (e.g., 5, 10, 20)')

    # Ignore any arguments that are not recognized
    args, unknown = parser.parse_known_args()

    running_times = {'SparseNJ': [], 'FNJ': [], 'RNJ': []}

    # Process each dataset
    for dataset in args.datasets:
        # Construct the full path to the dataset file
        dataset_path = os.path.join('data', f'{dataset}.npy')

        # Load the dataset
        data = np.load(dataset_path, allow_pickle=False)
        print(f"Dataset {dataset} loaded with shape {data.shape}")

        for algorithm in ['sparseNJ', 'FNJ', 'RNJ']:
            result = subprocess.run([sys.executable, f'{algorithm}.py', '--dataset', dataset_path],
                                    capture_output=True, text=True)
            print(f"{algorithm} stdout: {result.stdout}")
            print(f"{algorithm} stderr: {result.stderr}")
            output = result.stdout.strip().split('\n')

            if algorithm == 'sparseNJ':
                times = []
                for line in output:
                    if "SNJ clustering time=" in line:  # Changed from f"{algorithm} clustering time="
                        time_str = line.split('=')[1].strip().strip('[]')
                        times = [float(t) for t in time_str.split()]
                        break
                if times:
                    average_time = sum(times) / len(times)
                    running_times['SparseNJ'].append(average_time)  # Changed from algorithm to 'SparseNJ'

            else:
                time = None
                for line in output:
                    if f"{algorithm} clustering time=" in line:
                        time = float(line.split('=')[1].strip())
                        break
                if time is not None:
                    running_times[algorithm].append(time)

    # Plotting the results
    plt.figure(figsize=(12, 6))
    x = np.arange(len(args.datasets))
    width = 0.2

    offsets = [-width, 0, width]
    colors = ['#2ca02c', '#d62728', '#1f77b4']
    labels = ['SparseNJ', 'FNJ', 'RNJ']

    print("Running times:", running_times)

    for i, (label, offset, color) in enumerate(zip(labels, offsets, colors)):
        if running_times[label]:
            plt.bar(x + offset, running_times[label], width, label=label, color=color, edgecolor='black')
        else:
            print(f"No running times for {label}")

    plt.xlabel('Dataset', fontsize=16, labelpad=12)
    plt.ylabel('Elapsed Time (seconds)', fontsize=16, labelpad=12)
    plt.title('Elapsed Time Comparison', fontsize=18, pad=22)
    plt.xticks(x, args.datasets, rotation=45, ha="right")
    plt.yticks(fontsize=14)
    plt.legend(fontsize=14, frameon=True, shadow=True)
    plt.tight_layout()
    plt.savefig('elapsed_time_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
