#!/usr/bin/python

import sys
from SuperSet import *

usage = """
usage: %s input_fn
""" % sys.argv[0]

if len(sys.argv) != 2 :
  print usage
  sys.exit(-1)

#=================================================================
def getKey(seq) :
   mm = {}
   pos = 0
   while True :
     k = seq.find('[', pos)
     if k == -1 :
       break
     mm[str(k)] = '0'
     pos = k+4
   m = []
   for x in mm :
     m.append(x)
   m.sort()
   key = ' '.join(m)
   return key

#=================================================================
h = {}
a = open(sys.argv[1]).readlines()
#out2 = open(sys.argv[1] + "_mm_data", 'w')
out = open(sys.argv[1] + "_supersets", 'w')

# build hash table: keys are the query IDs,
# values are lists of the metadata for all hits for the query
last_line = ""
j=0
while j < len(a) :
  line = a[j]
  j+= 1
  if len(line) > 1 and line[0] != '#' and not line[0].isdigit():
    t = []
    key = ""
    if line[0] != " " :
      t = last_line.split()
      key = t[5]
      seq = line[:-1] #strip of the '\n'
      while True :
        c = a[j][0]
        if c == 'A' or c == 'G' or c == 'T' or c == 'C' or c == ']' or c == '[' :
          seq += a[j][:-1]
          j += 1
        else :
          break
      t.append(seq)
      if not h.has_key(key) :
        h[key] = []
      h[key].append(t)
  last_line = line

print 'num queries:', len(h.keys())
ct = 0
for x in h.values() :
  ct += len(x)
print 'num hits (sum over hits for the queries):', ct


supersets = []

ss = False

names = {}

for query in h.keys() :

 sets = {}
 for z in h[query] :
   if not names.has_key(z[1]) :
     names[z[1]] = 0
   names[z[1]] += 1

   key = getKey(z[8])
   if not sets.has_key(key) :
     ss = SS()
     ss.setKmer(z[8])
     sets[key] = ss
     supersets.append(ss)

 for z in h[query] :
   key = getKey(z[8])

   gid = ''
   for c in range(len(z[1])-1, -1, -1) :
     if z[1][c].isdigit() :
       gid += z[1][c]
   gid = gid[::-1]
   gid = int(gid)

   length = z[0]
   loc = z[2]

   #next, find mismatch positions
   seq = z[8]
   mm = {}
   pos = 0
   used = 0

   while True :
     k = seq.find('[', pos)
     if k == -1 : break
     pos = k+4
     mm[k-3*used] = seq[k+1]
     used += 1
   hit = Hit(gid, length, loc, mm)
   sets[key].insertHit(hit)

#ss = mergeSuperSets(supersets)
print "writing supersets; superset count: ", len(supersets)
for ss in supersets :
  out.write(ss.__str__())

out.close()
