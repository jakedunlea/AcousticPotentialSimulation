import numpy as np
import matplotlib.pyplot as plt
import AcousticPotential
import calculateAcousticPotential
import GUI
from scipy.interpolate import griddata
import os


def calculate_acoustophoretic_force(variables=None):
    """

    :param variables:
    :return:
    """

    var_dict = variables
    if var_dict is None:
        var_dict = GUI.VariableSelectionBox().variables

    ap_loc = 'saved_arrays\\2-D\\pressure\\{shape}\\{h}m height\\{s} step\\Acoustic Potential.npz'.format(
        shape=var_dict['ref_shape'],
        h=round(var_dict['h'], 3),
        s=var_dict['ss']
    )

    try:
        ap = np.load(ap_loc)['ap']
    except FileNotFoundError:
        inp = input("\n>>> File not found, run simulation with these parameters now? [y/n]"
                    "\n >> ")
        if inp == 'y':
            print('\n>>> Running simulation with specified parameters...')
            calculateAcousticPotential.calculate_acoustic_potential(variables=var_dict)
            try:
                ap = np.load(ap_loc)['ap']
            except FileNotFoundError:
                print("\n>>> Issue occurred.")
                raise FileNotFoundError
        else:
            raise FileNotFoundError

    sim_obj = AcousticPotential.Simulate(var_dict)

    acoustophoretic_force = sim_obj.acoustophoretic_force(ap)

    af_save_loc = 'saved_arrays\\2-D\\pressure\\{shape}\\{h}m height\\{s} step\\Acoustophoretic Force.npz'.format(
        shape=var_dict['ref_shape'],
        h=round(var_dict['h'], 3),
        s=var_dict['ss']
    )

    np.savez(af_save_loc, af=acoustophoretic_force)

    return acoustophoretic_force


def af_compare():
    """

    :return:
    """
    cvd = GUI.VariableSelectionBox().variables
    fvd = cvd.copy()
    cvd['ref_shape'] = 'c'
    fvd['ref_shape'] = 'f'
    c_af = calculate_acoustophoretic_force(cvd)
    f_af = calculate_acoustophoretic_force(fvd)

    def do_interpolate(vdict, af):

        sim_obj = AcousticPotential.Simulate(vdict)

        x = np.linspace(-sim_obj.ref_w/2, sim_obj.ref_w/2 + sim_obj.ss, 300)
        y = np.linspace(0, sim_obj.h, 300)

        xx, yy = np.meshgrid(x, y)

        # noinspection PyTypeChecker
        inted = griddata((xx, yy), af, (sim_obj.xi, sim_obj.yi), method='cubic')

        ys = AcousticPotential.pos_on_semicircle(x, sim_obj.R, sim_obj.c_xy)
        ymask = yy < ys
        return inted

    fig, ax = plt.subplots(1, 2, True, True)

    mx = np.array([c_af, f_af]).max()
    mn = np.array([c_af, f_af]).min()

    im = ax[0].imshow(c_af[0], origin='lower', vmax=mx, vmin=mn)
    ax[1].imshow(f_af[0], origin='lower', vmax=mx, vmin=mn)

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax)

    ax[0].set_title('Curved Reflector')
    ax[1].set_title('Planar Reflector')

    # ax[0].quiver(c_af[1][::5, ::5], c_af[0][::5, ::5])
    # ax[1].quiver(f_af[1][::5, ::5], f_af[0][::5, ::5])

    plt.show()

if __name__ == '__main__':
    af_compare()
