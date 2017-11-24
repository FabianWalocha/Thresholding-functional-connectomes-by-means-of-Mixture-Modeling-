"""
Circular graph used to represent connections between different Brain ROIs
Produces: Figure 4, figure 5, figure 9, figure 10, figure 11
Made using python 3.5.3.1Qt5
Author: Fabian Walocha (2017), fawalocha@gmail.com
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import matplotlib.cm as cm
from mne.viz import circular_layout, plot_connectivity_circle


# Create colormap for the difference map (blue to black to red)
cdict1 = {'red':   ((0.0, 0.0, 0.0),
           (0.5, 0.0, 0.1),
           (1.0, 1.0, 1.0)),

 'green': ((0.0, 0.0, 0.0),
           (1.0, 0.0, 0.0)),

 'blue':  ((0.0, 0.0, 1.0),
           (0.5, 0.1, 0.0),
           (1.0, 0.0, 0.0))
}
plt.register_cmap(name='BlueRed1', data=cdict1)

names = [['Empirical precision',
          'Empirical precision',
          'Empirical precision',
          'Empirical precision'],
         ['Ledoit-Wolf',
          'Ledoit-Wolf',
          'Ledoit-Wolf',
          'Ledoit-Wolf'],
         ['Permutation testing',
          'Permutation testing',
          'Permutation testing',
          'Permutation testing'],
         ['Mixture modeling',
          'Mixture modeling',
          'Mixture modeling',
          'Mixture modeling'],
          ['Proportional thresholding, 5%',
           'Proportional thresholding, 5%',
           'Proportional thresholding, 5%',
           'Proportional thresholding, 5%'],
           ['Proportional thresholding, 10%',
           'Proportional thresholding, 10%',
           'Proportional thresholding, 10%',
           'Proportional thresholding, 10%']]
names2 = ['EP','ledW','permT','MM','N1','N2']


# Load labels
label_names = sio.loadmat("data/rois.mat")
label_names = np.squeeze(label_names['rois'])
label_names_tmp = []
for ind3 in range(36):
    label_names_tmp.append(label_names[ind3][0])
label_names = label_names_tmp
lh_labels = label_names[0:18]
rh_labels = label_names[18:36]

for num,idx1 in enumerate(['EP','ledW','permT','MM','N1','N2']):
    for num2,idx2 in enumerate(['rest','task','diff','diff3']):

        # load data
        con_res = sio.loadmat("results/res_"+idx1+idx2+".mat")
        if num2 == 0 or num2 == 1:
            con_res = con_res['res_'+idx1+'3']
        else:
            con_res = con_res['res_'+idx1]

        # Order the nodes to be symmetrical
        node_order = list()
        node_order.extend(lh_labels[::-1])  # reverse the order
        node_order.extend(rh_labels)
        
        # Order color palette to reflect symmetry
        colors = iter(iter(cm.rainbow(np.linspace(0, 1, len(lh_labels)))))
        node_colors = []
        for idx3 in range(len(lh_labels)):
            node_colors.append(next(colors));
        node_colors.extend(node_colors);
        
        node_angles = circular_layout(label_names, node_order, start_pos=90,
                                      group_boundaries=[0, len(label_names) / 2])
        

        # Difference map
        if  idx2 == 'diff' or idx2 == 'diff3':
            [fig1,_] = plot_connectivity_circle(con_res, label_names, n_lines=None,
                                     node_angles=node_angles, node_colors=node_colors,
                                     colormap = 'BlueRed1',vmin = -99,vmax = 99,
                                     title=names[num][num2])
        
        # Heatmaps
        else:
            
            [fig1,_] = plot_connectivity_circle(con_res, label_names, n_lines=None,
                                     node_angles=node_angles, node_colors=node_colors,
                                     colormap = 'afmhot', vmin = 0, vmax = 207,
                                     title=names[num][num2])

        fig1.savefig('plots/circular/circle'+names2[num]+idx2+'.png', dpi=300, facecolor='black',bbox_inches='tight')
        
        #
        fig1.show()
