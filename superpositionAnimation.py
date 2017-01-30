#! python3

"""
Matplotlib Animation Example

author: Jake Vanderplas
email: vanderplas@astro.washington.edu
website: http://jakevdp.github.com
license: BSD
Please feel free to use and modify this, but keep the above information. Thanks!

Jake Dunlea:
This code was just used (badly) as a personal tool to visualise the superposition of opposing waves and the effect
of frequency diff and phase diff on this superposition
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib.widgets import Slider


fig2, ax2 = plt.subplots(3, 1, sharex=True)

for ax in ax2:
    ax.set_xlim([0, 1])

line2, = ax2[0].plot([], [], lw=2, color='red')
line3, = ax2[1].plot([], [], lw=2, color='green')
line4, = ax2[2].plot([], [], lw=2, color='yellow')
h_line2, = ax2[0].plot([], [], lw=0.5, color='white')
h_line3, = ax2[1].plot([], [], lw=0.5, color='white')
h_line4, = ax2[2].plot([], [], lw=0.5, color='white')

plt.subplots_adjust(bottom=0.2)

axcolor = 'grey'
axfreq = plt.axes([0.2, 0.05, 0.65, 0.03], axisbg=axcolor)
axph = plt.axes([0.2, 0.1, 0.65, 0.03], axisbg=axcolor)

f1 = 10
freqdiff = 0
f2 = f1 + freqdiff
sfreq = Slider(axfreq, 'Freq Diff', 0, 30.0, valinit=freqdiff, color='red')
sfreq.label.set_color('white')

phfac = 0.05
phdiff = 0
ph1 = np.pi * phfac
ph2 = np.pi * (phfac + phdiff)
sph = Slider(axph, 'Phase Diff', 0, 5.0, valinit=phdiff, color='red')
sph.label.set_color('white')


def init2():
    line2.set_data([], [])
    line3.set_data([], [])
    line4.set_data([], [])
    h_line2.set_data([], [])
    h_line3.set_data([], [])
    h_line4.set_data([], [])
    return line2, line3, line4, h_line2, h_line3, h_line4,


def animate2(i):
    t = np.linspace(0, 1, 500)

    y1 = np.sin(2 * np.pi * f1 * t - ph1 * i)
    y2 = np.sin(2 * np.pi * f2 * t + ph2 * i)
    line2.set_data(t, y1)
    line3.set_data(t, y2)
    line4.set_data(t, y1+y2)
    h_line2.set_data(t, 0)
    h_line3.set_data(t, 0)
    h_line4.set_data(t, 0)
    return line2, line3, line4, h_line2, h_line3, h_line4,

anim2 = animation.FuncAnimation(fig2, animate2, init_func=init2,
                                frames=200, interval=50, blit=True)

transFigure = fig2.transFigure.inverted()

coord1 = transFigure.transform(ax2[0].transData.transform([0.1, 0.5]))
coord2 = transFigure.transform(ax2[2].transData.transform([0.1, -0.5]))

# line = matplotlib.lines.Line2D((coord1[0], coord2[0]), (coord1[1], coord2[1]),
#                                transform=fig2.transFigure, color='white')

# fig2.lines = line,
fig2.patch.set_facecolor('black')

for a in ax2:
    a.set_axis_bgcolor('black')
    a.set_ylim([-2.25, 2.25])
    a.tick_params(axis='x', colors='white')
    a.tick_params(axis='y', colors='white')
    a.yaxis.label.set_color('white')
    a.xaxis.label.set_color('white')

    for spine in a.spines.values():
        spine.set_color('white')

plt.show()
