# This python code is used to draw a ternary diagrams

import matplotlib.pyplot as plt
import numpy as np
import math
filename = 'Directionality.txt'  # read Directionality.txt file
size=np.loadtxt(filename, delimiter=',',dtype=float)

coor_x=[-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,-2,-1.5,-1,-0.5,0,0.5,1,1.5,-1,-0.5,0,0.5,1,1.5,2,0,0.5,1,1.5,2,2.5,1,1.5,2,2.5,3,2,2.5,3,3.5,3,3.5,4,4,4.5,5]
coor_y=[0,0.866025404,1.732050808,2.598076211,3.464101615,4.330127019,5.196152423,6.062177826,6.92820323,7.794228634,0,0.866025404,1.732050808,2.598076211,3.464101615,4.330127019,5.196152423,6.062177826,6.92820323,0,0.866025404,1.732050808,2.598076211,3.464101615,4.330127019,5.196152423,6.062177826,0,0.866025404,1.732050808,2.598076211,3.464101615,4.330127019,5.196152423,0,0.866025404,1.732050808,2.598076211,3.464101615,4.330127019,0,0.866025404,1.732050808,2.598076211,3.464101615,0,0.866025404,1.732050808,2.598076211,0,0.866025404,1.732050808,0,0.866025404,0]
plt.figure(figsize=(9.11, 11.35), dpi=600)
 
plt.scatter(coor_x, coor_y, c=size, s=3400, marker='h',cmap='coolwarm')
cb=plt.colorbar(orientation="horizontal",shrink=0.7)
cb.ax.tick_params(labelsize=16)
font = {'family' : 'Times New Roman',
	'color'  : 'black',
	'weight' : 'normal',
	'size'   : 22,
	}
cb.set_label('Directionality',fontdict=font)  # set figure labels
plt.yticks([])
plt.axis('off')
plt.savefig("Directionality")   # set the figure file name and save it
