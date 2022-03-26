import csv
import numpy as np
import ROOT

reader = csv.reader(open("fev11_cob_rotate_chip_channel_x_y_mapping.txt"), delimiter=" ")
writer = csv.writer(open('cob_tmp2.txt', 'w'), delimiter=" ")

chip=-1
x0=9999
y0=9999
# channel=-1
x=9999
y=9999

mapxy_chip = ROOT.TH2F("mapxy_chip","map-xy; x; y",32,-90,90,32,-90,90)

x_tmp = [None]*64
y_tmp = [None]*64

channel = np.zeros((15,8,8))

for ichip, ix, iy in np.ndindex(channel.shape):
	channel[ichip, ix, iy] = -1

for i, line in enumerate(reader):

	if i == 0:
		continue

	input_string = [None]*5

	chip = int(line[0])
	x0   = float(line[1])
	y0   = float(line[2])
	x    = float(line[3])
	y    = float(line[4])




'''
for i, line in enumerate(reader):

	input_string = [None]*6
	
	if i < 1:
		continue

	chip    = int(line[0])
	x0      = float(line[1])
	y0      = float(line[2])
	# channel = int(line[3])
	x       = float(line[4])
	y       = float(line[5])

	if chip == 11:
		x_tmp[channel] = x
		y_tmp[channel] = y

	if chip == 13:
		# print(x_tmp[channel])
		x = x_tmp[channel] + 46.1
		# x = -x
		x = "{:.1f}".format(x)
		y = y_tmp[channel]

	input_string[0] = str(chip)
	input_string[1] = str(x0)
	input_string[2] = str(y0)
	input_string[3] = str(channel)
	input_string[4] = str(x)
	input_string[5] = str(y)

	writer.writerow(input_string)
'''

