#!python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import GUI
import AcousticPotential
import runSimulation
import runAnimation
import os


def calculate_acoustic_potential(variables=None):
    """
    Calculates and saves Acoustic Levitation Potential of a sound field
    :param variables: Dictionary of sound field parameters to use in simulation
    :return: Acoustic Levitation Potential array
    """

    var_dict = variables
    if var_dict is None:
        var_dict = GUI.VariableSelectionBox().variables

    p_arr_loc = 'saved_arrays\\2-D\\pressure\\{shape}\\{h}m height\\{s} step\\Pressure Arrays.npz'.format(
        shape=var_dict['ref_shape'],
        h=round(var_dict['h'], 3),
        s=var_dict['ss']
    )

    # Tries to load a dictionary of pressure arrays for the location outlined by the variable selected

    try:
        pressure_dict = np.load(p_arr_loc)
    except FileNotFoundError:
        inp = input("\n>>> File not found, run simulation with these parameters now? [y/n]"
                    "\n >> ")
        if inp == 'y':
            print('\n>>> Running simulation with specified parameters...')
            runSimulation.simulate(variables=var_dict, save=True)
            try:
                pressure_dict = np.load(p_arr_loc)
            except FileNotFoundError:
                print("\n>>> Issue occurred.")
                raise FileNotFoundError
        else:
            raise FileNotFoundError

    p_titles = pressure_dict.keys()
    p_arr = [pressure_dict[p] for p in p_titles]
    total_p_field = np.array(p_arr).sum(axis=0)
    p_arr.append(total_p_field)

    # Get sim_obj and use functions in AcousticPotential to get Acoustic Potential and plot and save
    sim_obj = AcousticPotential.Simulate(var_dict)

    ap_save_loc = 'saved_arrays\\2-D\\pressure\\{shape}\\{h}m height\\{s} step\\Acoustic Potential.npz'.format(
        shape=var_dict['ref_shape'],
        h=round(var_dict['h'], 3),
        s=var_dict['ss']
    )

    try:
        ap_file = np.load(ap_save_loc)
        acoustic_potential = ap_file['ap']
    except FileNotFoundError:
        print("\n>>> File not found, calculating Acoustic Levitation Potential now.")

        if not os.path.exists(os.path.dirname(ap_save_loc)):
            os.makedirs(os.path.dirname(ap_save_loc))

        acoustic_potential = sim_obj.levitation_potential(total_p_field)

        np.savez(ap_save_loc, ap=acoustic_potential)

    print('\n>>> AP min =', np.real(acoustic_potential).min(),
          '\n>>> AP Max =', np.real(acoustic_potential).max())

    plt.imshow(np.real(acoustic_potential),
               origin='lower')
    plt.colorbar()
    plt.show()

    return acoustic_potential

if __name__ == '__main__':
    ap = calculate_acoustic_potential()
