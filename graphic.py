#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd

from main import N, SIZE


def clear_list(data):
    return [x for x in data if x != '']

def transform_data(data):
    data = clear_list(data)
    frames = np.zeros(N)
    for frame in data:
        frame = clear_list(frame.split('\n'))
        frames = np.vstack([frames, frame])
    frames = np.delete(frames, 0, 0)
    return frames

def read_data(file):
    with open(file, 'r', encoding='utf-8') as reader:
        lines = reader.read()
    data_splited = lines.split('frame\n')
    data_splited = transform_data(data_splited)
    return data_splited

data = read_data('coords.txt')
print(np.shape(data))

frames = list(data.copy())
for i, frame in enumerate(frames):
    frames[i] = [r.split() for r in frame]
frames = np.array(frames, dtype=np.float32)
#frames = fix_values(frames)

s = (1, N*len(frames))
frames_counter = np.reshape(frames[:,:,0], s)[0]
x_counter = np.reshape(frames[:,:,1], s)[0]
y_counter = np.reshape(frames[:,:,2], s)[0]
z_counter = np.reshape(frames[:,:,3], s)[0]

fs = pd.Series(frames_counter)
xs = pd.Series(x_counter)
ys = pd.Series(y_counter)
zs = pd.Series(z_counter)

newdf = {'frame' : fs, 'x' : xs, 'y' : ys, 'z' : zs}

df_data = pd.DataFrame(newdf)
df_data.to_csv('dfdata.csv')


import plotly.express as px

df = df_data
fig = px.scatter_3d(df, x='x', y='y', z='z', animation_frame='frame')
fig.update_layout(xaxis=dict(range=[0,SIZE]), yaxis=dict(range=[0,SIZE]))
fig.show()

#atom1 = frames[:,0,:]
#dfa = pd.DataFrame(atom1, columns=['frame', 'x', 'y', 'z'])
#print(dfa.info())
#dfa.to_csv('atom1.csv')

#fig = px.scatter_3d(dfa, x='x', y='y', z='z', animation_frame='frame')
#fig.update_layout(xaxis_range=[0,SIZE], yaxis_range=[0,SIZE], zaxis_range=[0,SIZE])

#fig = px.scatter(dfa, x='x', y='y')
#fig.update_layout(xaxis_range=[0,SIZE], yaxis_range=[0,SIZE])

#fig.show()
