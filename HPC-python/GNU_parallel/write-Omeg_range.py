import numpy as n
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("Omeg_range_str")

args = parser.parse_args()
X = args.Omeg_range_str.split(',')
start = float(X[0])
stop  = float(X[1])
inc = float(X[2])
with open("Omeg_range.txt", "w") as db :
  for Omeg in n.arange(start,stop,inc): #5.01,7.02,0.1
    db.writelines(f"{Omeg}\n")