#!/home/jiahuih/anaconda3/bin/python3

import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

# define the path
script_path = "./pdb_parser_modified.py"
pdb_directory = "/home/jiahuih/data/jiahuih/SLC/Resolution/SLC6_HotSpots/AFMv2_3_output_dir/SURF4_SLC6A4"


# define the name for the input
pdb_files = os.listdir(pdb_directory)
pdb_files = [file for file in pdb_files if file.startswith("ranked_") 
            and file.endswith(".pdb")]

# Define chain name needed in the PDB_parser.py
chain1 = input('First chain: ')
chain2 = input('Second chain: ')

# Loop through each filtered .pdb file and parallelize the process
def run_script(pdb_file):
    pdb_path = os.path.join(pdb_directory, pdb_file)
    command = f"python {script_path} {pdb_path} {chain1} {chain2}"
    subprocess.run(command, shell=True)
with ThreadPoolExecutor() as executor:
    executor.map(run_script, pdb_files)




