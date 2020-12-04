from pathlib import Path
import argparse

# Read WT and proposed residue abbreviations
parser = argparse.ArgumentParser(description='print post-execution states for each job in job array')
parser.add_argument('input_file', type=str, help='path to post-execution states file')
args = parser.parse_args()

# Read and parse file
with open(args.input_file, 'r') as f:
	lines = [line for line in f.readlines()]
	for i, line in enumerate(lines):
		if line[:3] == 'Job':
			print(line.split(',')[0])
		if "RUNLIMIT" in line:
			print(lines[i-2])