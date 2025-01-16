# fasterq-dump SRR14814229 --skip-technical --split-files --progress --threads 4
# python3 script-sra-call.py SRR14814229

import argparse
import subprocess

parser = argparse.ArgumentParser(description = 'Download and convert SRA files to FASTQ format')
# group = parser.add_mutually_exclusive_group(required=True)
parser.add_argument(
    '-s', '--sra', nargs = '+', help = 'One or more SRA numbers to process'
)
# group.add_argument(
#     '-f', '--file', type = str,
#     help = 'Path to CSV file containing SRA numbers in the first column'
# )
parser.add_argument(
    '-o', '--output', type = str, default = './sra_data',
    help = 'Output directory for downloaded files (default: ./sra_data)'
)

subprocess.run(
    ['fasterq-dump', sra_number, 
        '--outdir', output_dir,
        '--skip-technical',
        '--split-files',  
        '--progress',
        '--threads', '4' 
    ],
    check = True,
    capture_output = True,
    text = True
)
