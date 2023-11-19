import numpy as np
import math
import matplotlib.pyplot as plt

# bit SNR in dB
X = [1, 1.5, 2, 2.5, 3, 3.5]

# my SCL: each with 500 errorblock
L2 = [[0.3348962, 0.1537988, 0.0820345, 0.0282072, 0.0107582, 0.0036057],
      [0.3190810, 0.1567398, 0.0760803, 0.0297106, 0.0106687, 0.0038410],
      [0.3006615, 0.1557632, 0.0792267, 0.0299904, 0.0108439, 0.0035432]]

L4 = [[0.2596054, 0.1372495, 0.0630199, 0.0262826, 0.0104351, 0.0033085],
      [0.2624672, 0.1336898, 0.0589831, 0.0229684, 0.0096658, 0.0031515],
      [0.2574665, 0.1303441, 0.0646245, 0.0239120, 0.0092415, 0.0033088]]

L8 = [[0.2358491, 0.1220405, 0.0589067, 0.0232515, 0.0087504, 0.0031934],
      [0.2385496, 0.1214182, 0.0519049, 0.0229442, 0.0087428, 0.0032571],
      [0.2338634, 0.1276813, 0.0546388, 0.0259565, 0.0090719, 0.0035773]]

L16 = [[0.2361833, 0.1177302, 0.0538851, 0.0223574, 0.0090181, 0.0034127],
       [0.2264493, 0.1191895, 0.0508440, 0.0233252, 0.0086050, 0.0033119],
       [0.2427184, 0.1284357, 0.0589414, 0.0239923, 0.0092297, 0.0030990]]

L32 = [[0.2248201, 0.1305824, 0.0541594, 0.0243025, 0.0099065, 0.0034298],
       [0.2188184, 0.1194743, 0.0569087, 0.0235206, 0.0093139, 0.0035924],
       [0.2295684, 0.1223092, 0.0562810, 0.0218924, 0.0090705, 0.0033591]]

# Take average
l2 = []
for j in range(len(L2[0])):
    ave = 0
    for i in range(len(L2)):
        ave += L2[i][j]
    ave = ave / len(L2)
    l2.append(ave)
l4 = []
for j in range(len(L4[0])):
    ave = 0
    for i in range(len(L4)):
        ave += L4[i][j]
    ave = ave / len(L4)
    l4.append(ave)
l8 = []
for j in range(len(L8[0])):
    ave = 0
    for i in range(len(L8)):
        ave += L8[i][j]
    ave = ave / len(L8)
    l8.append(ave)
l16 = []
for j in range(len(L16[0])):
    ave = 0
    for i in range(len(L16)):
        ave += L16[i][j]
    ave = ave / len(L16)
    l16.append(ave)
l32 = []
for j in range(len(L32[0])):
    ave = 0
    for i in range(len(L32)):
        ave += L32[i][j]
    ave = ave / len(L32)
    l32.append(ave)
# SCL ref from Po-Chung (50 error block)
y2 = [0.314, 0.189, 0.0789, 0.0269, 0.0115, 0.00397]
y4 = [0.286, 0.133, 0.0571, 0.0232, 0.00936, 0.00394]
y8 = [0.266, 0.121, 0.0546, 0.0227, 0.00993, 0.00394]
y16 = [0.251, 0.118, 0.0546, 0.0219, 0.00993, 0.00394]
y32 = [0.251, 0.118, 0.0546, 0.0219, 0.00993, 0.00394]

plt.figure()
line0, = plt.plot(X, l2, 'ko-', markersize = 3, linewidth = 1)
line1, = plt.plot(X, l4, 'y^-', markersize = 4, linewidth = 1)
line2, = plt.plot(X, l8, 'bx-', markersize = 4, linewidth = 1)
line3, = plt.plot(X, l16, 'gv-', markersize = 4, linewidth = 1)
line4, = plt.plot(X, l32, 'r+-', markersize = 4, linewidth = 1)
plt.xlabel("BSNR (dB)")
plt.ylabel("BLER")
#plt.title("N = 128  K = 64    Error block = 500")
plt.grid('--', which='both', axis='both')
plt.yscale("log")
#plt.xticks([1, 2, 3, 4])
#plt.xlim([1, 4])
plt.legend((line0, line1, line2, line3, line4),
    ('L = 2', 'L = 4', 'L = 8', 'L = 16', 'L = 32'))
plt.show()

