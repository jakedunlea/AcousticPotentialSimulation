#! python3

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import numpy as np
import GUI
import runSimulation


def write_pressure_animation(variables=None):
    """
    Loads, animates and writes to file the components of a pressure wave (pd, pr1, pr2, pr3)
    as well as the pressure wav itself. The location to which the animations are written
    depends on the parameters with which the pressure wave was simulated.
    :param variables: a dictionary of keywords with the parameters with which the pressure wave is
    :return:
    """
    inp = input("\n>>> This will override current pressure animation files, continue? [y/n]")
    var_dict = variables
    if var_dict is None:
        var_dict = GUI.VariableSelectionBox().variables

    arr_loc = 'saved_arrays\\2-D\\pressure\\{shape}\\{h}m height\\{s} step\\Pressure Arrays.npz'.format(
        shape=var_dict['ref_shape'],
        h=round(var_dict['h'], 3),
        s=var_dict['ss']
    )

    # Tries to load a dictionary of pressure arrays for the location outlined by the variable selected
    pressure_dict = None

    try:
        pressure_dict = np.load(arr_loc)
    except FileNotFoundError:
        inp = input("\n>>> File not found, run simulation with these parameters now? [y/n]"
                    "\n >> ")
        if inp == 'y':
            print('\n>>> Running simulation with specified parameters...')
            runSimulation.simulate(variables=var_dict)
            try:
                pressure_dict = np.load(arr_loc)
            except FileNotFoundError:
                print("\n>>> Issue occurred.")
                raise FileNotFoundError
        else:
            raise FileNotFoundError

    p_titles = ['pd', 'pr1', 'pr2', 'pr3']
    p_arr = [pressure_dict[p] for p in p_titles]

    p_arr.append(np.sum(p_arr, 0))

    anim_titles = ['Direct Pressure Wave', 'First Reflected Pressure Wave',
                   'Second Reflected Pressure Wave', 'Third Reflected Pressure Wave',
                   'Superimposed Pressure Wave']

    flat_arr = np.array(p_arr).flatten()
    all_min = np.min(np.real(flat_arr))
    all_max = np.max(np.real(flat_arr))

    freq = var_dict['freq']
    period = 1 / freq
    w = 2 * np.pi * freq

    def cycle_sound(time, p_field):
        a = 1j * np.array(np.exp(1j * w * time))
        im_arr = p_field * a
        return np.real(im_arr)

    for idx, p in enumerate(p_arr):
        fig = plt.figure()

        t = 0
        im = plt.imshow(cycle_sound(t, p),
                        origin='lower',
                        vmin=all_min, vmax=all_max,
                        animated=True)

        def updateim(*args):
                global t
                t += period / 50
                im.set_array(cycle_sound(t, p))
                return im,

        frame = plt.gca()

        frame.axes.xaxis.set_ticklabels([])
        frame.axes.yaxis.set_ticklabels([])

        meta = dict(title=anim_titles[idx], artist='Jake Dunlea')

        ani = animation.FuncAnimation(fig, updateim, frames=500, interval=50, blit=True)

        print('\n>>> Writing {p}'.format(p=anim_titles[idx]))
        my_writer = animation.FFMpegWriter(fps=20, metadata=meta, bitrate=720)
        ani.save('animations\\c\\{p}.mp4'.format(p=anim_titles[idx]), writer=my_writer)


def write_acoustic_potential_animation(acoustic_potential, variables=None):
    """
    Loads, animates and writes to file the changes in acoustic potential per period.
    The location to which the animations are written depends on the parameters with which
    the pressure wave was simulated.
    :param acoustic_potential: Acoustic Levitation Potential array (complex)
    :param variables: a dictionary of keywords with the parameters with which the pressure wave
                      was simulated.
    :return:
    """
    var_dict = variables
    if var_dict is None:
        var_dict = GUI.VariableSelectionBox().variables

    freq = var_dict['freq']
    period = 1 / freq
    w = 2 * np.pi * freq

    def cycle_sound(time, ap):
        a = 1j * np.array(np.exp(1j * w * time))
        im_arr = ap * a
        return np.real(im_arr)

    fig = plt.figure()

    t = 0
    im = plt.imshow(cycle_sound(t, acoustic_potential),
                    origin='lower',
                    animated=True)

    def updateim(*args):
            global t
            t += period / 50
            im.set_array(cycle_sound(t, acoustic_potential))
            return im,

    frame = plt.gca()

    frame.axes.xaxis.set_ticklabels([])
    frame.axes.yaxis.set_ticklabels([])

    meta = dict(title='Acoustic Levitation Potential Simulation',
                artist='Jake Dunlea')

    ani = animation.FuncAnimation(fig, updateim, frames=500, interval=50, blit=True)

    print('\n>>> Writing {p}'.format(p='Acoustic Levitation Potential Simulation'))
    my_writer = animation.FFMpegWriter(fps=20, metadata=meta, bitrate=720)
    ani.save('animations\\{shape}\\{p}.mp4'.format(shape=var_dict['ref_shape'],
                                                   p='Acoustic Levitation Potential Simulation'
                                                   ),
             writer=my_writer)

    plt.show()

if __name__ == '__main__':
    import calculateAcousticPotential
    ap = calculateAcousticPotential.calculate_acoustic_potential()
    t = 0
    write_acoustic_potential_animation(ap)
