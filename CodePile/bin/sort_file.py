#!/usr/bin/python
import sys
dir = sys.argv[1]
topology = sys.argv[2]

if (dir[-1] != '/'):
    dir = dir + '/'
file_in = open(dir+topology+"relations.rlt", "r")
file_out = open(dir+topology+"relations_sort.rlt", "w")
data = []
for line in file_in:
    if (line.split()[0] == '#'):
        file_out.write(line)
        continue
    data.append(line.split())
def cmp1 (a):
    return int(a[0])
data.sort(key=cmp1)
for i in data:
    file_out.write('\t'.join(i))
    file_out.write("\n")
file_in.close()
file_out.close()

file_in = open(dir+topology+"states.states", "r")
file_out = open(dir+topology+"states_sort.states", "w")
data = []
for line in file_in:
    if (line.split()[0] == '#'):
        file_out.write(line)
        continue
    data.append(line.split())
def cmp1 (a):
    return int(a[0])
data.sort(key=cmp1)
for i in data:
    file_out.write('\t'.join(i))
    file_out.write("\n")
file_in.close()
file_out.close()