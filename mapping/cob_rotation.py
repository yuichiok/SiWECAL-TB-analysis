import csv
import numpy as np
# import matplotlib.pyplot as plt
import ROOT

reader = csv.reader(open("backup/fev11_cob_chip_channel_x_y_mapping.txt"), delimiter=" ")

MyFile = ROOT.TFile("output.root","RECREATE");

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

channel_13_rot180 = np.rot90(channel_13,2)
channel_15_rot180 = np.rot90(channel_15,2)

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

# writer_tmp = csv.writer(open('cob_tmp1.txt', 'w'), delimiter=" ")
with open('cob_tmp1.txt', 'w') as outfile1:

	writer_tmp = csv.writer(outfile1, delimiter=" ")

	for i, line in enumerate(reader):

		if i == 0:
			writer_tmp.writerow(line)
		else:
			input_string = [None]*6

			chip = int(line[0])
			x0   = float(line[1])
			y0   = float(line[2])
			chnl = int(line[3])
			x    = float(line[4])
			y    = float(line[5])

			input_string[0] = str(chip).zfill(2)
			input_string[1] = str(x0)
			input_string[2] = str(y0)

			if chip == 13:
				binx = mapxy_13_rot.GetXaxis().FindBin(-x)
				biny = mapxy_13_rot.GetYaxis().FindBin(-y)
				ch = int(mapxy_13_rot.GetBinContent(binx,biny))
				input_string[3] = str(ch).zfill(2)
			elif chip == 15:
				binx = mapxy_15_rot.GetXaxis().FindBin(-x)
				biny = mapxy_15_rot.GetYaxis().FindBin(-y)
				ch = int(mapxy_15_rot.GetBinContent(binx,biny))
				input_string[3] = str(ch).zfill(2)
			else:
				input_string[3] = str(chnl).zfill(2)

			input_string[4] = str(x)
			input_string[5] = str(y)

			writer_tmp.writerow(input_string)


txt = open("cob_tmp1.txt").readlines()
# open("cob_tmp2.txt","w").write("\n".join(sorted(txt.split("\n"))))
with open("cob_tmp2.txt","w") as f:
	f.write("".join(sorted(txt)))

reader = csv.reader(open("cob_tmp2.txt"), delimiter=" ")
# writer = csv.writer(open('cob_tmp3.txt', 'w'), delimiter=" ")

with open('cob_tmp3.txt', 'w') as outfile2:

	writer = csv.writer(outfile2, delimiter=" ")

	line0 = "chip x0 y0 channel x y"
	writer.writerow(line0.split(" "))

	for i, line in enumerate(reader):

		if line[0]=='chip':
			continue
		else:
			input_string = [None]*6

			# print(line)
			chip = int(line[0])
			x0   = float(line[1])
			y0   = float(line[2])
			chnl = int(line[3])
			x    = float(line[4])
			y    = float(line[5])

			input_string[0] = str(int(chip))
			input_string[1] = str(x0)
			input_string[2] = str(y0)
			input_string[3] = str(int(chnl))
			input_string[4] = str(x)
			input_string[5] = str(y)

			writer.writerow(input_string)

ROOT.gStyle.SetOptStat(0)
MyFile.cd()
mapxy_chip.Write()
mapxy_13.Write()
mapxy_13_rot.Write()
mapxy_15.Write()
mapxy_15_rot.Write()

