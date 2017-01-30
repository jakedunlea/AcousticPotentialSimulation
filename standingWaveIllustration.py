#! Python3

import numpy as np
import matplotlib.pyplot as plt

f = 20
wl = 1
w = 2 * np.pi * f
k = 2*np.pi / wl
x = np.linspace(0, 4, 500)
t = 0
phs = [-np.pi, np.pi]


def wave_func(amp, wavenum, dist, ang_freq, time, phase):
    u = amp * np.sin(wavenum * dist - ang_freq * time + phase)
    return u

waves = [wave_func(1, k, x, w, t, ph) for ph in phs]
waves.append(waves[0] + waves[1])

fig, ax = plt.subplots(3, 1, sharex=True)
color = plt.cm.gist_rainbow(np.linspace(0, 1, 4))

for idx, a in enumerate(ax):
    a.plot(x, waves[idx], color=color[idx])
    a.plot(x, np.zeros_like(x), color='black')
    a.set_axis_bgcolor('white')
    a.set_ylim([-2.25, 2.25])
    a.xaxis.set_ticklabels([])
    a.tick_params(axis='both', colors='black')
    y = a.yaxis
    labels = [-3, -2, -1, 0, 1, 2, 3]
    print(labels)
    y.set_ticklabels(labels, color='black')

    for spine in a.spines.values():
        spine.set_color('black')

plt.show()
