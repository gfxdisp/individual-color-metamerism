import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import loadmat
import sys
import os
sys.path.append('./ColorVideoVDP')
from pycvvdp.dolby_ictcp import ictcp
import torch
import json
# Read data from CSV files
all_data = pd.read_csv('ObserverMetamerism/ObserverMetamerism/Data/YC_AllData.csv')
ciexyz31_1 = pd.read_csv('ObserverMetamerism/ObserverMetamerism/Data/AuxData/ciexyz31_1.csv', header=None).to_numpy()
obs_1000_cmf_struct = loadmat("ObserverMetamerism/ObserverMetamerism/Data/AuxData/IndividualObs_2deg.mat")
obs_1000_cmfs = obs_1000_cmf_struct['xyz_CMFs']
display_1_spd = pd.read_csv('ObserverMetamerism/ObserverMetamerism/Data/Spectra/C2_Spectra.csv', header=None).to_numpy()
display_2_spd = pd.read_csv('ObserverMetamerism/ObserverMetamerism/Data/Spectra/X310_Spectra.csv', header=None).to_numpy()
display_3_spd = pd.read_csv('ObserverMetamerism/ObserverMetamerism/Data/Spectra/Projector_Spectra.csv', header=None).to_numpy()
display_4_spd = pd.read_csv('ObserverMetamerism/ObserverMetamerism/Data/Spectra/VG246_Spectra.csv', header=None).to_numpy()
display_spds = [display_1_spd, display_2_spd, display_3_spd, display_4_spd]

exp_num = len(all_data)
color_rgb_s = pd.read_csv("ObserverMetamerism/ObserverMetamerism/Data/all_11_colors_rgb.csv", header=None).to_numpy() / 256.0
display_patterns = np.array([[1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]])
display_configurations = np.array([[-1, 1, 2, 3], [-1, -1, 4, 5], [-1, -1, -1, 6], [-1, -1, -1, -1]])
num_obs = obs_1000_cmfs.shape[2]
e_cell = np.zeros((6, 11), dtype=object)
score_cell = np.empty((6, 11), dtype=object)

deltaE_ITP_class = ictcp()
# Calculate E_Cell
for display_pattern_index in range(6):
    display_1 = display_patterns[display_pattern_index, 0]
    display_2 = display_patterns[display_pattern_index, 1]
    for color_index in range(11):
        display_spd_1 = display_spds[display_1 - 1]
        display_spd_2 = display_spds[display_2 - 1]
        E_set = []
        for obs in range(num_obs):
            X1, X2, Y1, Y2, Z1, Z2 = 0, 0, 0, 0, 0, 0
            for lamda in range(390,785,5):
                X1 = X1 + display_spd_1[color_index, lamda - 380] * obs_1000_cmfs[int((lamda - 390) / 5), 0, obs] * 5
                X2 = X2 + display_spd_2[color_index, lamda - 380] * obs_1000_cmfs[int((lamda - 390) / 5), 0, obs] * 5
                Y1 = Y1 + display_spd_1[color_index, lamda - 380] * obs_1000_cmfs[int((lamda - 390) / 5), 1, obs] * 5
                Y2 = Y2 + display_spd_2[color_index, lamda - 380] * obs_1000_cmfs[int((lamda - 390) / 5), 1, obs] * 5
                Z1 = Z1 + display_spd_1[color_index, lamda - 380] * obs_1000_cmfs[int((lamda - 390) / 5), 2, obs] * 5
                Z2 = Z2 + display_spd_2[color_index, lamda - 380] * obs_1000_cmfs[int((lamda - 390) / 5), 2, obs] * 5
            XYZ1 = torch.tensor([X1, Y1, Z1])[None,:,None,None,None]
            XYZ2 = torch.tensor([X2, Y2, Z2])[None,:,None,None,None]
            delta_E_ITP = deltaE_ITP_class.predict_xyz(XYZ1, XYZ2)
            E_set.append(delta_E_ITP)
        e_cell[display_pattern_index, color_index] = E_set

# with open(r'e_cell.json', 'r') as fp:
#     data = json.load(fp)
#     e_cell = np.array(data)

# Initialize variables
for i in range(6):
    for j in range(11):
        score_cell[i, j] = []
# Populate score_Cell
for i in range(exp_num):
    display_1 = int(all_data.loc[i, "Display1"])
    display_2 = int(all_data.loc[i, "Display2"])
    color_index = int(all_data.loc[i, "ColorIndex"])
    display_pattern_index = display_configurations[display_1 - 1, display_2 - 1]
    score = all_data.loc[i, "MatchingScore"]
    score_cell[display_pattern_index - 1, color_index - 1].append(score)

# Calculate means and std
e_means = np.zeros((6, 11))
e_std = np.zeros((6, 11))
score_means = np.zeros((6, 11))
score_std = np.zeros((6, 11))
point_color = np.zeros((6, 11, 3))
point_display_pattern = np.zeros((6, 11))

for i in range(6):
    for j in range(11):
        if all(score == 0 for score in score_cell[i, j]):
            e_means[i, j] = np.nan
            e_std[i, j] = np.nan
            score_means[i, j] = np.nan
            score_std[i, j] = np.nan
        else:
            e_means[i, j] = np.mean(e_cell[i, j])
            e_std[i, j] = np.std(e_cell[i, j])
            score_means[i, j] = np.mean(score_cell[i, j])
            score_std[i, j] = np.std(score_cell[i, j])

        point_color[i, j, :] = color_rgb_s[j, :]
        point_display_pattern[i, j] = i

# Flatten the arrays for plotting
e_means_flat = e_means.flatten()
e_std_flat = e_std.flatten()
score_means_flat = score_means.flatten()
score_std_flat = score_std.flatten()
point_color_flat = point_color.reshape((-1, 3))
point_display_pattern_flat = point_display_pattern.flatten()

valid_indices = ~np.isnan(e_means_flat) & ~np.isnan(score_means_flat)

# Plot the scatter plot with error bars
plt.figure()
plt.scatter(e_std_flat[valid_indices], score_std_flat[valid_indices], c=point_color_flat[valid_indices], edgecolors=point_color_flat[valid_indices])
for indice in range(valid_indices.size):
    if not valid_indices[indice]:
        continue
    plt.text(e_std_flat[indice], score_std_flat[indice] - 0.02,
             "{}, {}".format(display_patterns[int(point_display_pattern_flat[indice]), 0],
                             display_patterns[int(point_display_pattern_flat[indice]), 1]),
             horizontalalignment='center')

for i in range(11):
    x = 1.4
    y = - (i - 9) * 0.02 + 0.95
    plt.scatter(x + 0.27, y, s=60, c=[color_rgb_s[i, :]], edgecolors=[color_rgb_s[i, :]])
    plt.text(x , y, 'Color {}'.format(i + 1), horizontalalignment='left', verticalalignment='center')

plt.text(1.4, 0.89, 'Display 1 - LG C2', horizontalalignment='left', verticalalignment='center')
plt.text(1.4, 0.87, 'Display 2 - Sony X310', horizontalalignment='left', verticalalignment='center')
plt.text(1.4, 0.85, 'Display 3 - Samsung Laser Projector', horizontalalignment='left', verticalalignment='center')
plt.text(1.4, 0.83, 'Display 4 - ASUS VG246', horizontalalignment='left', verticalalignment='center')

plt.xlim([0, 2])
plt.xlabel('deltaE 2000 (Standard Deviation of 1000 values)')
plt.ylabel('Score (Standard Deviation of 74 values)')
plt.title('74 Subject Scores and 1000 CMFs deltaE 2000 Standard deviation points')
# plt.savefig('E:/About_Cambridge/All Research Projects/Color/ObserverMetamerism/ObserverMetamerism/Data/deltaE_subjective_pure_std_no_zero.png')
plt.show()
