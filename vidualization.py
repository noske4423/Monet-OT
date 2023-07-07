import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pylab as pl
import anndata as ad
import math
import matplotlib.colors as mcolors
import ot
from matplotlib import ticker, cm, colors
import pickle
import random
from statistics import stdev, variance, median
from PIL import Image
import scipy.stats
import define_function as d_f
import csv
from datetime import datetime
import os
import gc
import multiprocessing as mp
import time
import psutil

organ_list = ['Blood', 'BM', 'Brain', 'Liver', 'LN', 'SmallInt', 'Spleen']
disease_list = ['WT', 'AD']
color_list = ['dodgerblue', 'orange', 'limegreen', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan', 'black']

celltype_ref = {'Blood': 'MCA_broad_cell_type',
                'BM': 'Tabula.Muris_cell_type',
                'Brain': 'Azimuth_predicted.subclass',
                'Liver': 'Tabula.Muris_cell_type',
                'LN': 'MCA_Blood_broad_cell_type',
                'SmallInt': 'MCA_broad_cell_type',
                'Spleen': 'Tabula.Muris_cell_type'
                }
ms_path = '/Users/kataokayuunosuke/MS'
image_path = f'{ms_path}/image'

i=0
organ_position_list = [0]
line_position_list = [-0.5]
organ_celltype_list = []
for organ in organ_list:
    # load celltype list
    celltype_list_path = f'{ms_path}/{organ}_cell_type_sample_list.txt'
    with open(celltype_list_path, mode='r') as f:
        celltype_list = f.read().splitlines()
    for celltype in celltype_list:
        organ_celltype_list.append(celltype)
        i += 1

    organ_position = (line_position_list[-1] + i -0.5)/2
    if (organ_position*2)%2 == 0:
        organ_position = int(organ_position)

    organ_position_list.append(organ_position)
    line_position_list.append(i - 0.5)

num_array = np.arange(1, 53, 1)
position_list = []
label_name_list = []
for i in num_array:
    position_list.append(i-1)
    if int(i-1) in organ_position_list and i != 1:
        organ_index = organ_position_list.index(int(i-1))
        label_name_list.append(organ_list[organ_index-1]+'   '+str(int(i)))

    else:
        label_name_list.append(str(int(i)))

for organ_position in organ_position_list:
    if type(organ_position) == float:
        organ_index = organ_position_list.index(organ_position)
        label_name_list.append(organ_list[organ_index-1]+'       ')
        position_list.append(organ_position)

# save position_list
with open(f'{ms_path}/organomix/position_list.txt', mode='w') as f:
    f.write('\n'.join(map(str, position_list)))

# load position_list
with open(f'{ms_path}/organomix/position_list.txt', mode='r') as f:
    position_list = f.read().splitlines()
    position_list = list(map(float, position_list))

# save label_name_list
with open(f'{ms_path}/organomix/label_name_list.txt', mode='w') as f:
    f.write('\n'.join(map(str, label_name_list)))

# load label_name_list
with open(f'{ms_path}/organomix/label_name_list.txt', mode='r') as f:
    label_name_list = f.read().splitlines()

# save line_position_list
with open(f'{ms_path}/organomix/line_position_list.txt', mode='w') as f:
    f.write('\n'.join(map(str, line_position_list)))

# load line_position_list
with open(f'{ms_path}/organomix/line_position_list.txt', mode='r') as f:
    line_position_list = f.read().splitlines()
    line_position_list = list(map(float, line_position_list))





# データの生成
data = np.random.random((52, 52))

# フィギュアの生成
fig, ax1 = plt.subplots(figsize=(10, 10))

# データをプロット
img = ax1.imshow(data)

# 線を描く（6x6のブロックに分割）
for i in range(6):
    ax1.plot([-0.5, 51.5], [line_position_list[i+1], line_position_list[i+1]], color='black', linewidth=1)
    ax1.plot([line_position_list[i+1], line_position_list[i+1]], [-0.5, 51.5], color='black', linewidth=1)

# 元の軸のラベルを設定
xticks = position_list
yticks = position_list

ax1.set_xticks(xticks)
ax1.set_yticks(yticks)

# ラベルを設定
xlabels = label_name_list
ylabels = label_name_list

ax1.set_xticklabels(xlabels, rotation=-90)
ax1.set_yticklabels(ylabels)

# x軸のメモリを上側に表示
ax1.xaxis.tick_top()

# メモリの線を消す
ax1.tick_params(axis='both', length=0)

plt.show()

# 保存
fig.savefig(f'{ms_path}/organomix/test.png', dpi=300, bbox_inches='tight', pad_inches=0.05)