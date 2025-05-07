import time
import random
import subprocess
import matplotlib.pyplot as plt

def generate_random_sequence(length):
    return ''.join(random.choices('ACGT', k=length))

def measure_execution_time_results(script_path, query, target):
    start = time.monotonic_ns()
    process = subprocess.run(
        ['python', script_path, '--A', query, '--B', target],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    stop = time.monotonic_ns()
    if process.returncode != 0:
        print(f"Error running {script_path}: {process.stderr.decode()}")
    return stop - start

def main():
    query_length = 42
    target_sizes = [100, 200, 500, 1000, 2000, 5000, 10000, 50000, 100000]
    query = generate_random_sequence(query_length)

    ryan_script = 'sw_ryan_layer.py'
    vincent_script = 'sw_vincent_bowen.py'
    biopython_script = 'sw_biopython.py'
    sci_kit_script = 'sw_skbio.py'

    ryan_times = []
    vincent_times = []
    biopython_times = []
    sci_kit_times = []

    for target_size in target_sizes:
        target = generate_random_sequence(target_size)

        ryan_time = measure_execution_time_results(ryan_script, query, target)
        vincent_time = measure_execution_time_results(vincent_script, query, target)
        biopython_time = measure_execution_time_results(biopython_script, query, target)
        sci_kit_time = measure_execution_time_results(sci_kit_script, query, target)

        ryan_times.append(ryan_time)
        vincent_times.append(vincent_time)
        biopython_times.append(biopython_time)
        sci_kit_times.append(sci_kit_time)

    plt.figure(figsize=(10, 6))
    plt.plot(target_sizes, ryan_times, label="Ryan Layer's Provided Smith-Waterman Algorithm")
    plt.plot(target_sizes, vincent_times, label="Vincent Bowen's Rudimentary Smith-Waterman Algorithm")
    plt.plot(target_sizes, biopython_times, label="Biopython Pairwise Sequence Alignment")
    plt.plot(target_sizes, sci_kit_times, label="SciKit-Bio Striped Smith-Waterman")
    plt.xlabel('Target Size')
    plt.ylabel('Execution Time (nanoseconds)')
    plt.title('Sequence Alignment Algorithm Speed Comparison')
    plt.legend()
    plt.grid(True)
    plt.savefig('results.png')

if __name__ == "__main__":
    main()