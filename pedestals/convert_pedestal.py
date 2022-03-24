import csv

reader = csv.reader(open("Pedestal_run_000011_highgain.txt"), delimiter=" ")

layer=-1
chip=-1
channel=-1

for i, line in enumerate(reader):
	
	if i < 2:
		continue
	elif i < 10:
		layer   = line[0]
		chip    = line[1]
		channel = line[2]



