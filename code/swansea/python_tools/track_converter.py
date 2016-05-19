#!/usr/bin/python
import csv
points = {}
with open('t3d_tide-jmj_real_gen.slf.tracks.txt', 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter="\t")
    for row in reader:
        if len(row) == 2:
            continue
        try:
            points[row[1]].append(row)
        except KeyError:
            points[row[1]] = []
            points[row[1]].append(row)

for i in points.keys():
    w = open('t3d_tide-jmj_real_gen.slf.tracks.t%d.csv' % int(i), 'wb+')
    writer = csv.writer(w, delimiter=",")
    for x in points[i]:
        writer.writerow(x)
    w.close()

