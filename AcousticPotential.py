import numpy as np
from scipy.interpolate import griddata
import progressbar
from tkinter import filedialog
import tkinter as tk
import os


def distance(xy1, xy2):
    """
    Calculates the distance between two xy coordinates.
    :param xy1: XY coordinates of point 1 (x1, y1)
    :param xy2: XY coordinates of point 2 (x2, y2)
    :return: Float of distance
    """
    x_dist = xy2[0] - xy1[0]
    y_dist = xy2[1] - xy1[1]
    dist = np.sqrt(x_dist ** 2 + y_dist ** 2)
    return dist


def pos_on_semicircle(x, r, cxy):
    """
    Calculates the corresponding y-coordinate of the point at the top of a circle given an x-coordinate
    :param x: x--coordinate to calculate corresponding y-coordinate at top of circle
    :param r: radius of circle
    :param cxy: coordinates of the centre of the circle (cx, cy)
    :return: float of y-coordinate
    """
    pos = np.sqrt(r ** 2 - (x - cxy[0]) ** 2) + cxy[1]

    return pos


def save_array_as(arr, loc=None):
    """
    Saves array as .npy file
    If optional param loc is not provided a dialog box opens to determine location
    :param arr: array to be saved
    :param loc: optional location to save array
    :return: array saved at loc
    """
    if loc is None:
        root = tk.Tk()
        root.loc = filedialog.asksaveasfilename(initialdir="/",
                                                title="Save as",
                                                filetypes=(("npy files", "*.npy"),
                                                           ("all files", "*.*"))
                                                )
        np.save(root.loc, arr)
        root.destroy()
    else:
        np.save(loc, arr)


def load_array(loc=None):
    """
    Loads NumPy array from a .npy or .npz file location
    :param loc: Optional location to NumPy array
                If not provided, dialog opens to find location
    :return: NumPy array
    """
    if loc is None:
        root = tk.Tk()
        root.loc = filedialog.askopenfilename(initialdir="/",
                                              title="Save as",
                                              filetypes=(("npy files", "*.npy"),
                                                         ("npz files", "*.npz"),
                                                         ("all files", "*.*"))
                                              )
        ret = np.load(root.loc)
        root.destroy()
    else:
        ret = np.load(loc)

    return ret


class Simulate(object):
    def __init__(self, value_dict, obs_dims=2, save=False):

        if value_dict is not None and value_dict != {}:
            self.save = save
            self.obs_dims = obs_dims
            self.value_dict = value_dict

            self.ref_shape = value_dict['ref_shape']
            self.freq = value_dict['freq']
            self.period = 1 / self.freq
            self.ang_freq = value_dict['ang_freq']
            self.c_air = value_dict['c_air']
            self.dens = value_dict['dens']
            self.wl = value_dict['wl']
            self.k = value_dict['k']
            self.h = value_dict['h']
            self.R = value_dict['R']
            self.c_xy = value_dict['c_xy']
            self.trans_dia = value_dict['trans_dia']
            self.ref_w = value_dict['ref_w']
            self.ss = value_dict['ss']

            self.N = int(self.trans_dia / 0.002)
            self.I = int(self.ref_w / 0.002)

            self.N_coords = self.n_coords()['N_coords']
            self.I_coords = self.i_coords()['I_coords']

            self.M_coords = self.m_coords(obs_dims=self.obs_dims)
            self.M = len(self.M_coords)
            self.M_x, self.M_y = list(zip(*self.M_coords))
            xi, yi = np.linspace(min(self.M_x), max(self.M_x), 300), np.linspace(min(self.M_y), max(self.M_y), 300)
            self.xi, self.yi = np.meshgrid(xi, yi)

            self.NM_dist_mat = self.nm_dist_mat()
            self.NI_dist_mat = self.ni_dist_mat()
            self.IM_dist_mat = self.im_dist_mat()

            self.A_n = self.n_coords()['A_n']
            self.A_i = self.i_coords()['A_i']

        else:
            print('wtf')

    def nm_dist_mat(self):
        """
        Produces a array of the distance between the points on the transducer face (n) and
        the points in the space between the transducer and the reflector (m).
        :return: N x M array of floats
        """
        mat = np.zeros([self.N, self.M])
        for n in range(self.N):
            for m in range(self.M):
                mat[n, m] = distance(self.N_coords[n], self.M_coords[m])
        return mat

    def ni_dist_mat(self):
        """
        Produces a array of the distance between the points on the transducer face (n) and
        the points on the reflector face (i).
        :return: N x I array of floats
        """
        mat = np.zeros([self.N, self.I])
        for n in range(self.N):
            for i in range(self.I):
                mat[n, i] = distance(self.N_coords[n], self.I_coords[i])
        return mat

    def im_dist_mat(self):
        """
        Produces a array of the distance between the points on the reflector face (i) and
        the points in the space between the transducer and the reflector (m).
        :return: I x M array of floats
        """
        mat = np.zeros([self.I, self.M])
        for i in range(self.I):
            for m in range(self.M):
                mat[i, m] = distance(self.I_coords[i], self.M_coords[m])
        return mat

    def n_coords(self):
        """
        Calculates coordinates of centre point of each cell on the transducer face from the chosen parameters.
        :return: Dictionary:
                    'trans_x':     array of x-coordinates of edges of cells
                    'A_n':         array of areas of each cell on the transducer face
                    'N_coords':    array of coordinates of centres of cells
        """
        trans_x = np.arange(-self.trans_dia / 2, self.trans_dia / 2 + 0.002, 0.002)
        a_n = [(trans_x[n + 1] - trans_x[n]) / 2 for n in range(self.N)]
        cx_n = [trans_x[n] + (trans_x[n + 1] - trans_x[n]) / 2 for n in range(self.N)]
        coords = [(x, 0) for x in cx_n]
        d = {'trans_x': trans_x, 'A_n': a_n, 'N_coords': coords}
        return d

    def i_coords(self):
        """
        Calculates coordinates of centre point of each cell on the reflector face from the chosen parameters.
        :return: Dictionary:
                    'ref_x':       array of x-coordinates of edges of cells
                    'A_i':         array of areas of each cell on the reflector face
                    'I_coords':    array of coordinates of centres of cells
                    'cx_i':        array of x-coordinates of centres of cells
        """
        ref_x = np.arange(-self.ref_w / 2, self.ref_w / 2 + 0.002, 0.002)

        if self.ref_shape == 'c':  # Curved reflector
            dist_coords1 = [(ref_x[i], pos_on_semicircle(ref_x[i], self.R, self.c_xy)) for i in range(self.I)]
            dist_coords2 = [(ref_x[i + 1], pos_on_semicircle(ref_x[i + 1], self.R, self.c_xy)) for i in range(self.I)]
            a_i = [distance(dist_coords1[i], dist_coords2[i]) for i in range(self.I)]

            cx_i = [ref_x[i] + (ref_x[i + 1] - ref_x[i]) / 2 for i in range(self.I)]
            cy_i = [pos_on_semicircle(x, self.R, self.c_xy) for x in cx_i]
            i_coords = list(zip(cx_i, cy_i))
        else:  # Flat reflector
            a_i = [(ref_x[i + 1] - ref_x[i]) / 2 for i in range(self.I)]
            cx_i = [ref_x[i] + (ref_x[i + 1] - ref_x[i]) / 2 for i in range(self.I)]
            i_coords = [(x, self.h) for x in cx_i]
        d = {'ref_x': ref_x, 'A_i': a_i, 'I_coords': i_coords, 'cx_i': cx_i}

        return d

    def m_coords(self, obs_dims=2):
        """
        Calculates coordinates of points in space to calculate sound pressure at.
        :return: Array of coordinates in space between transducer nd reflector
                    (If obs_dims == 1 then coordinates form a straight line (z-axis))
        """
        if obs_dims == 2:
            d = self.i_coords()
            cx_i = d['cx_i']
            print('\n>>>  Discretizing space')
            coords = []
            heights = np.arange(self.ss, self.h, self.ss)
            for i in range(self.I):
                for height in heights:
                    if height >= pos_on_semicircle(i, self.R, self.c_xy):
                        break
                    else:
                        coords.append((cx_i[i], height))
            print('\n>>> Space has been discretized'
                  '\n >> {m} points created'.format(m=len(coords)))

        else:
            print('\n>>> Discretizing z-axis')
            coords = []
            heights = np.arange(self.ss, self.h, self.ss)
            for height in heights:
                if height >= pos_on_semicircle(0, self.R, self.c_xy):
                    break
                else:
                    coords.append((0, round(height, 5)))
            print('\n>>> Z-axis has been discretized'
                  '\n >> {m} points have been created along z-axis'.format(m=len(coords)))

        return np.array(coords)

    def pd(self):
        print('\n>>> Calculating Components of Direct Pressure Wave '
              '\n >> {x} Components to calculate\n'.format(x=(self.N * self.M)))
        total = self.N * self.M
        bar = progressbar.ProgressBar(max_value=total)
        cnt = 0
        press = np.ones([self.N, self.M], dtype=complex)
        for n in range(self.N):
            for m in range(self.M):
                dist = self.NM_dist_mat[n, m]
                press[n, m] = self.A_n[n] * np.exp(-1j * self.k * dist) / np.array(dist)
                bar.update(cnt)
                cnt += 1
        press *= (self.dens * self.c_air / self.wl)
        press = np.sum(press, 0)
        if self.obs_dims == 2:
            # noinspection PyTypeChecker
            press = griddata(self.M_coords, press, (self.xi, self.yi), method='cubic')

        print('\n\n>>> Completed Calculation of Direct Pressure Wave')

        if self.save:
            save_loc = 'saved_arrays\\{d}-D\\pressure\\{rs}\\{h}m height\\{ss} step\\PD @ M-coords Complex.npy'.format(
                rs=self.ref_shape,
                h=round(self.h, 3),
                ss=self.ss,
                d=self.obs_dims
            )

            if not os.path.isdir(os.path.split(save_loc)[0]):
                os.makedirs(os.path.split(save_loc)[0])

            np.save(save_loc, press)

        return press

    def pr1(self):
        print('\n>>> Calculating Components of First Reflected Pressure Wave'
              '\n >> {x} Components to calculate\n'.format(x=(self.N * self.I * self.M)))
        total = self.N * self.I * self.M
        bar = progressbar.ProgressBar(max_value=total)
        cnt = 0
        press = np.ones([self.N, self.I, self.M], dtype=complex)
        for n in range(self.N):
            for i in range(self.I):
                for m in range(self.M):
                    dist = self.NI_dist_mat[n, i] + self.IM_dist_mat[i, m]
                    press[n, i, m] = self.A_n[n] * self.A_i[i] * np.exp(-1j * self.k * dist) / np.array(dist)
                    bar.update(cnt)
                    cnt += 1
        press *= ((self.dens * self.c_air / self.wl) * (1j / self.wl))
        press = np.sum(press, 0)
        press = np.sum(press, 0)
        if self.obs_dims == 2:
            # noinspection PyTypeChecker
            press = griddata(self.M_coords, press, (self.xi, self.yi), method='cubic')

        print('\n\n>>> Completed Calculation of First Reflected Pressure Wave')

        if self.save:
            save_loc = 'saved_arrays\\{d}-D\\pressure\\{rs}\\{h}m height\\{ss} step\\PR1 @ M-coords Complex.npy'.format(
                rs=self.ref_shape,
                h=round(self.h, 3),
                ss=self.ss,
                d=self.obs_dims
            )

            if not os.path.isdir(os.path.split(save_loc)[0]):
                os.makedirs(os.path.split(save_loc)[0])

            np.save(save_loc, press)

        return press

    def pr2(self):
        print('\n>>> Calculating Components of Second Reflected Pressure Wave'
              '\n >> {x} Components to calculate\n'.format(x=(self.N * self.I * self.N * self.M)))
        total = self.N * self.I * self.N * self.M
        bar = progressbar.ProgressBar(max_value=total)
        cnt = 0
        press = np.ones([self.N, self.I, self.N, self.M], dtype=complex)
        for n in range(self.N):
            for i in range(self.I):
                for n2 in range(self.N):
                    for m in range(self.M):
                        dist = self.NI_dist_mat[n, i] + self.NI_dist_mat[n2, i] + self.NM_dist_mat[n2, m]
                        press[n, i, n2, m] = self.A_n[n] * self.A_i[i] * self.A_n[n2] * np.exp(
                            -1j * self.k * dist) / np.array(dist)
                        bar.update(cnt)
                        cnt += 1
        press *= ((self.dens * self.c_air / self.wl) * (1j / self.wl) ** 2)
        press = np.sum(press, 0)
        press = np.sum(press, 0)
        press = np.sum(press, 0)

        if self.obs_dims == 2:
            # noinspection PyTypeChecker
            press = griddata(self.M_coords, press, (self.xi, self.yi), method='cubic')

        print('\n\n>>> Completed Calculation of Second Reflected Pressure Wave')

        if self.save:
            save_loc = 'saved_arrays\\{d}-D\\pressure\\{rs}\\{h}m height\\{ss} step\\PR2 @ M-coords Complex.npy'.format(
                rs=self.ref_shape,
                h=round(self.h, 3),
                ss=self.ss,
                d=self.obs_dims
            )

            if not os.path.isdir(os.path.split(save_loc)[0]):
                os.makedirs(os.path.split(save_loc)[0])

            np.save(save_loc, press)

        return press

    def pr3(self):
        print('\n>>> Calculating Components of Third Reflected Pressure Wave'
              '\n >> {x} Components to calculate\n'.format(x=((self.N ** 2) * (self.I ** 2) * self.M)))
        total = (self.N ** 2) * (self.I ** 2) * self.M
        bar = progressbar.ProgressBar(max_value=total)
        cnt = 0
        press = np.ones([self.N, self.I, self.N, self.I, self.M], dtype=complex)
        for n in range(self.N):
            for i in range(self.I):
                for n2 in range(self.N):
                    for i2 in range(self.I):
                        for m in range(self.M):
                            dist = self.NI_dist_mat[n, i] + self.NI_dist_mat[n2, i] \
                                   + self.NI_dist_mat[n2, i2] + self.IM_dist_mat[i2, m]
                            press[n, i, n2, i2, m] = self.A_n[n] * self.A_i[i] * self.A_n[n2] * self.A_i[i2] \
                                                     * np.exp(-1j * self.k * dist) / np.array(dist)
                            bar.update(cnt)
                            cnt += 1
        press *= ((self.dens * self.c_air / self.wl) * (1j / self.wl) ** 3)
        press = np.sum(press, 0)
        press = np.sum(press, 0)
        press = np.sum(press, 0)
        press = np.sum(press, 0)

        if self.obs_dims == 2:
            # noinspection PyTypeChecker
            press = griddata(self.M_coords, press, (self.xi, self.yi), method='cubic')

        print('\n\n>>> Completed Calculation of Third Reflected Pressure Wave')

        if self.save:
            save_loc = 'saved_arrays\\{d}-D\\pressure\\{rs}\\{h}m height\\{ss} step\\PR3 @ M-coords Complex.npy'.format(
                rs=self.ref_shape,
                h=round(self.h, 3),
                ss=self.ss,
                d=self.obs_dims
            )

            if not os.path.isdir(os.path.split(save_loc)[0]):
                os.makedirs(os.path.split(save_loc)[0])

            np.save(save_loc, press)

        return press

    def pressure_to_velocity_potential(self, p_field):
        vp = -p_field / (1j * self.ang_freq * self.dens)
        return vp

    @staticmethod
    def velocity_potential_to_velocity(velocity_potential):
        grad = np.array(np.gradient(np.array(velocity_potential)))
        vel_x, vel_y = grad[1], grad[0]
        return vel_x, vel_y

    @staticmethod
    def velocity_magnitude(vx, vy):
        return np.sqrt(vx**2 + vy**2)

    def levitation_potential(self, input_p_field):
        """

        :param input_p_field:
        :return:
        """
        n_samples = 4 ** 3
        times = np.linspace(0, self.period, n_samples)

        def p_at_tm(p_field, tm):
            a = 1j * np.array(np.exp(1j * self.ang_freq * tm))
            complex_p = p_field * a
            return complex_p

        pressures = v_potentials = list(np.zeros(n_samples).tolist())

        for idx, time in enumerate(times):
            pressures[idx] = p_at_tm(input_p_field, time)

            v_potentials[idx] = self.pressure_to_velocity_potential(pressures[idx])

        pressures = np.array(pressures)
        v_potentials = np.array(v_potentials)

        rms_pressure = np.sqrt(pressures.sum(0) / n_samples)

        rms_v_potential = np.sqrt(v_potentials.sum(0) / n_samples)
        rms_vx, rms_vy = self.velocity_potential_to_velocity(rms_v_potential)

        rms_v_mag = self.velocity_magnitude(rms_vx, rms_vy)

        acoustic_potential = (rms_v_mag / (3 * self.dens * self.c_air**2)) - ((self.dens * rms_pressure) / 2)

        return acoustic_potential

    @staticmethod
    def acoustophoretic_force(acoustic_potential):
        """

        :param acoustic_potential:
        :return:
        """

        af = np.real(np.array(np.gradient(acoustic_potential)))

        return af
