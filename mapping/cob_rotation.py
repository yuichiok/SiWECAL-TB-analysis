import csv
import numpy as np
# import matplotlib.pyplot as plt
import ROOT

reader = csv.reader(open("backup/fev11_cob_chip_channel_x_y_mapping.txt"), delimiter=" ")
writer = csv.writer(open('cob_tmp2.txt', 'w'), delimiter=" ")

chip=-1
x0=9999
y0=9999
# channel=-1
x=9999
y=9999

mapxy_chip   = ROOT.TH2F("mapxy_chip","map-xy; x; y",32,-90,90,32,-90,90)

mapxy_13     = ROOT.TH2F("mapxy_13","map-xy; x; y",8,-45,0,8,-90,-45)
mapxy_13_rot = ROOT.TH2F("mapxy_13_rot","map-xy; x; y",8,-45,0,8,-90,-45)

mapxy_15     = ROOT.TH2F("mapxy_15","map-xy; x; y",8,-90,-45,8,-90,-45)
mapxy_15_rot = ROOT.TH2F("mapxy_15_rot","map-xy; x; y",8,-90,-45,8,-90,-45)

channel_13 = np.zeros((8,8))
channel_15 = np.zeros((8,8))

for ix, iy in np.ndindex(channel_13.shape):
	channel_13[ix, iy] = -1
	channel_15[ix, iy] = -1

for i, line in enumerate(reader):

	if i == 0:
		continue

	input_string = [None]*6

	chip = int(line[0])
	x0   = float(line[1])
	y0   = float(line[2])
	chnl = int(line[3])
	x    = float(line[4])
	y    = float(line[5])

	if chip == 13:
		mapxy_chip.Fill(-x,-y,chnl)
		mapxy_13.Fill(-x,-y,chnl)
	elif chip == 15:
		mapxy_chip.Fill(-x,-y,chnl)
		mapxy_15.Fill(-x,-y,chnl)		


for ix, iy in np.ndindex(channel_13.shape):
	ch = int(mapxy_13.GetBinContent(ix+1,iy+1))
	channel_13[ix,iy] = ch

	ch = int(mapxy_15.GetBinContent(ix+1,iy+1))
	channel_15[ix,iy] = ch

print(channel_13)
print(channel_15)

channel_13_rot180 = np.rot90(channel_13,2)
channel_15_rot180 = np.rot90(channel_15,2)

print(channel_13_rot180)
print(channel_15_rot180)

for ix, iy in np.ndindex(channel_13.shape):
	cx = mapxy_13_rot.GetXaxis().GetBinCenter(ix+1)
	cy = mapxy_13_rot.GetYaxis().GetBinCenter(iy+1)
	ch = int(channel_13_rot180[ix,iy])
	mapxy_13_rot.Fill(cx,cy,ch)

	cx2 = mapxy_15_rot.GetXaxis().GetBinCenter(ix+1)
	cy2 = mapxy_15_rot.GetYaxis().GetBinCenter(iy+1)
	ch2 = int(channel_15_rot180[ix,iy])
	mapxy_15_rot.Fill(cx2,cy2,ch2)


reader = csv.reader(open("backup/fev11_cob_chip_channel_x_y_mapping.txt"), delimiter=" ")

for i, line in enumerate(reader):

	if i == 0:
		writer.writerow(line)
	else:
		input_string = [None]*6

		chip = int(line[0])
		x0   = float(line[1])
		y0   = float(line[2])
		chnl = int(line[3])
		x    = float(line[4])
		y    = float(line[5])

		input_string[0] = str(chip)
		input_string[1] = str(x0)
		input_string[2] = str(y0)
		# input_string[3] = str(channel)

		if chip == 13:
			binx = mapxy_13_rot.GetXaxis().FindBin(-x)
			biny = mapxy_13_rot.GetYaxis().FindBin(-y)
			ch = int(mapxy_13_rot.GetBinContent(binx,biny))
			input_string[3] = str(ch)
		elif chip == 15:
			binx = mapxy_15_rot.GetXaxis().FindBin(-x)
			biny = mapxy_15_rot.GetYaxis().FindBin(-y)
			ch = int(mapxy_15_rot.GetBinContent(binx,biny))
			print(binx,biny,ch)
			input_string[3] = str(ch)
		else:
			input_string[3] = str(chnl)

		input_string[4] = str(x)
		input_string[5] = str(y)

		writer.writerow(input_string)




# c0 = ROOT.TCanvas("c0","c0",500,500)
# mapxy_chip.Draw("col,text")
# c0.Print("plots/map_chip.png")

# c1 = ROOT.TCanvas("c1","c1",500,500)
# mapxy_13.Draw("col,text")
# c1.Print("plots/map_13.png")

# c2 = ROOT.TCanvas("c2","c2",500,500)
# mapxy_15.Draw("col,text")
# c2.Print("plots/map_15.png")

# c3 = ROOT.TCanvas("c3","c3",500,500)
# mapxy_13_rot.Draw("col,text")
# c3.Print("plots/map_rot_13.png")

# c4 = ROOT.TCanvas("c4","c4",500,500)
# mapxy_15_rot.Draw("col,text")
# c4.Print("plots/map_rot_15.png")




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

