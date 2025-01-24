#! /home/jiahuih/anaconda3/bin/python3

import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime

def parse_pdb(pdbfile, chain):
    fpdb = open(pdbfile)
    dict_chain = {}
    for line in fpdb:
        if line[0:4] == 'ATOM':
            ch=line[21]
            if ch == chain:
                x=float(line[30:38])
                y=float(line[38:46])
                z=float(line[46:54])
                at=line[12:16].strip()
                res=line[17:20]
                n=line[22:26].strip()
                coord=[x,y,z]
                dict_chain[(n, res)] = dict_chain.get((n, res),{})
                dict_chain[(n, res)][at] = coord
    return dict_chain

def distance(chain1, chain2):
    final = {}
    for res1 in chain1.keys():
        for at1 in chain1[res1]:
            min_d = 3.5
            for res2 in chain2.keys():
                for at2 in chain2[res2]:
                    d = dist(chain1[res1][at1], chain2[res2][at2])
                    if d < min_d:
                        final[res1] = final.get(res1, {})
                        final[res1][at1] = final[res1].get(at1, {})
                        final[res1][at1][res2] = d , at2
                        min_d = d
    return final

def dist(at1,at2):
    d=np.sqrt((at1[0]-at2[0])**2+(at1[1]-at2[1])**2 + (at1[2]-at2[2])**2)
    return d

def print_interacting_atoms(pdbfile, di):
    n1, res1, n2, res2, at1, at2, dist = (list() for x in range(7))
    for x in di.keys():
        for y in di[x].keys():
            for z in di[x][y]:
                n1.append(x[0])
                res1.append(x[1])
                n2.append(z[0])
                res2.append(z[1])
                at1.append(y)
                at2.append(di[x][y][z][1])                    
                dist.append(di[x][y][z][0])
    tr = {'N-1' : n1, 'RES-1' : res1, 'Atom-1' : at1, 'N-2' : n2, 'RES-2' : res2, 'Atom-2' : at2}
    index = pd.MultiIndex.from_arrays([n1, res1, at1, n2, res2, at2], names=(tr.keys()))
    df = pd.DataFrame(data = {'Distance' : dist}, index = index)
    file_name = os.path.splitext(pdbfile)[0]
    output_file = f'interaction_csv/{file_name}.csv'
    df.to_csv(output_file)
    current_time = datetime.now()
    print(f'{file_name} saved at {current_time.strftime("%H:%M:%S")}')

if __name__ == '__main__':
    if len(sys.argv) >= 1:
        pdbfile = sys.argv[1]
        chain1 = sys.argv[2]
        chain2 = sys.argv[3]
        dict_chain1 = parse_pdb(pdbfile, chain1)
        dict_chain2 = parse_pdb(pdbfile, chain2)
        di = distance(dict_chain1, dict_chain2)
        print_interacting_atoms(pdbfile, di)
    else: print('parsepdb.py pdbfile')
