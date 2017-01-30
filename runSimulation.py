#! python3

import GUI
import AcousticPotential
import numpy as np
import os


def simulate(variables=None, save=True):
    """
    Runs simulation of the direct pressure and the reflected
    pressures to the third reflection.
    :param variables: Uses a dictionary of variables to run simulation;
                      If None, opens GUI box to select variables
    :param save: If True, saves to relevant directory based on params.
    :return: Dictionary of pressure arrays.
    """

    var_dict = variables
    if var_dict is None:
        var_dict = GUI.VariableSelectionBox().variables

    sim_obj = AcousticPotential.Simulate(var_dict)

    pressure_arrs = [sim_obj.pd(), sim_obj.pr1(),
                     sim_obj.pr2()]  # , sim_obj.pr3()]

    pressure_arr_titles = ['pd', 'pr1', 'pr2']  # , 'pr3']

    title_arr_dict = dict(zip(pressure_arr_titles, pressure_arrs))

    if save:
        save_loc = 'saved_arrays\\{d}-D\\pressure\\{rs}\\{h}m height\\{ss} step\\Pressure Arrays'.format(
                    rs=sim_obj.ref_shape,
                    h=round(sim_obj.h, 3),
                    ss=sim_obj.ss,
                    d=sim_obj.obs_dims
                )

        if not os.path.exists(os.path.dirname(save_loc)):
            os.makedirs(os.path.dirname(save_loc))

        np.savez(save_loc, **title_arr_dict)

    return dict(zip(pressure_arr_titles, pressure_arrs))

if __name__ == '__main__':
    simulate()
