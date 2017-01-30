import tkinter as tk
import numpy as np


class VariableSelectionBox(object):
    """
    var = self.variables
    var['ref_shape'] = self.shape.get()
    var['freq'] = self.freq_val.get()
    var['ang_freq'] = 2 * self.freq_val.get() * np.pi
    var['c_air'] = self.c_air
    var['dens'] = self.dens
    var['wl'] = self.c_air / self.freq_val.get()
    var['k'] = 2 * np.pi / var['wl']
    var['h'] = self.height.get() / 1000
    var['R'] = self.radius.get() / 1000
    var['c_xy'] = (0, round(var['h'] - var['R'], 6))
    var['trans_dia'] = self.trans_dia.get() / 1000
    var['ref_w'] = self.ref_w.get() / 1000
    var['ss'] = self.ss.get() / 1000
    """
    def __init__(self):
        self.variables = {}

        self.root = tk.Tk()
        self.root.wm_title("Values for Pressure Field Calculation")

        tk.Label(self.root,
                 text="What Reflector Shape is being used?",
                 justify='center').grid(row=0, sticky='s', pady=2, columnspan=4)

        self.shape = tk.StringVar()
        self.shape.set('c')

        curved_radio = tk.Radiobutton(text='Curved',
                                      variable=self.shape,
                                      value='c',
                                      indicatoron=0,
                                      width=25)
        curved_radio.grid(row=1, column=0, columnspan=2, sticky='e')

        flat_radio = tk.Radiobutton(text='Flat',
                                    variable=self.shape,
                                    value='f',
                                    indicatoron=0,
                                    width=25)
        flat_radio.grid(row=1, column=2, columnspan=2, stick='w')

        tk.Label(self.root).grid(row=2)

        tk.Label(self.root,
                 text='Operating Frequency (f):').grid(row=3, column=1, sticky='e')

        self.freq_val = tk.IntVar()
        self.freq_val.set(28000)

        freq_entry = tk.Entry(self.root,
                              textvariable=self.freq_val,
                              width=8)
        freq_entry.bind('<Return>', self.change_val)
        freq_entry.grid(row=3, column=2, sticky='w')
        tk.Label(self.root, text='Hz', justify='left').grid(row=3, column=2, sticky='s')

        tk.Label(self.root,
                 text='Angular Frequency (w):').grid(row=4, column=0, columnspan=2, sticky='e')

        self.ang_freq_str = tk.StringVar()
        self.ang_freq_str.set(str(round(self.freq_val.get() * 2 * np.pi, 3)) + ' rad/s')
        tk.Label(self.root,
                 textvariable=self.ang_freq_str).grid(row=4, column=2, stick='w')

        tk.Label(self.root, text='').grid(row=5)

        self.dens = 1.225
        tk.Label(self.root, text='Density of Air (dens):').grid(row=6, column=0, columnspan=2, sticky='e')
        tk.Label(self.root, text='1.225 kg/m^3').grid(row=6, column=2, columnspan=2, sticky='w')

        self.c_air = 343
        tk.Label(self.root, text='Speed of Sound in Air (c):').grid(row=7, column=0, columnspan=2, sticky='e')
        tk.Label(self.root, text='343 m/s').grid(row=7, column=2, columnspan=2, sticky='w')

        self.wl_str = tk.StringVar()
        self.wl_str.set(str(round((self.c_air / self.freq_val.get()) * 1000, 3)) + ' mm')
        tk.Label(self.root, text='Wavelength (wl):').grid(row=8, column=0, columnspan=2, sticky='e')
        tk.Label(self.root, textvariable=self.wl_str).grid(row=8, column=2, columnspan=2, sticky='w')

        self.k_str = tk.StringVar()
        self.k_str.set(str(round((2 * np.pi / (self.c_air / self.freq_val.get())) / 1000, 3)) + ' rad/mm')
        tk.Label(self.root, text='Wavenumber (k):').grid(row=9, column=0, columnspan=2, sticky='e')
        tk.Label(self.root, textvariable=self.k_str).grid(row=9, column=2, columnspan=2, sticky='w')

        tk.Label(self.root, text='').grid(row=10)

        tk.Label(self.root, text='Reflector-Transducer Dist (h):').grid(row=11, column=0, columnspan=2, sticky='e')
        self.height = tk.DoubleVar()
        self.height.set(round((self.c_air / self.freq_val.get()) * 1000 * 2.5, 3))
        height_entry = tk.Entry(self.root,
                                textvariable=str(self.height),
                                width=8)
        height_entry.grid(row=11, column=2, stick='w')
        height_entry.bind('<Return>', self.change_val)
        tk.Label(self.root, text='mm').grid(row=11, column=2, sticky='s')

        tk.Label(self.root, text='Radius of Curvature (R):').grid(row=12, column=0, columnspan=2, sticky='e')
        self.radius = tk.DoubleVar()
        self.radius.set(27.44)
        radius_entry = tk.Entry(self.root,
                                textvariable=str(self.radius),
                                width=8)
        radius_entry.grid(row=12, column=2, sticky='w')
        radius_entry.bind('<Return>', self.change_val)
        tk.Label(self.root, text='mm').grid(row=12, column=2, sticky='s')

        tk.Label(self.root, text='').grid(row=13)

        tk.Label(self.root, text='Transducer Diameter (td):').grid(row=14, column=0, columnspan=2, sticky='e')
        self.trans_dia = tk.DoubleVar()
        self.trans_dia.set(26)
        trans_dia_entry = tk.Entry(self.root,
                                   textvariable=str(self.trans_dia),
                                   width=8)
        trans_dia_entry.grid(row=14, column=2, sticky='w')
        trans_dia_entry.bind('<Return>', self.change_val)
        tk.Label(self.root, text='mm').grid(row=14, column=2, sticky='s')

        tk.Label(self.root, text='Reflector Width (rw):').grid(row=15, column=0, columnspan=2, sticky='e')
        self.ref_w = tk.DoubleVar()
        self.ref_w.set(45)
        ref_w_entry = tk.Entry(self.root,
                               textvariable=str(self.ref_w),
                               width=8)
        ref_w_entry.grid(row=15, column=2, stick='w')
        ref_w_entry.bind('<Return>', self.change_val)
        tk.Label(self.root, text='mm').grid(row=15, column=2, sticky='s')

        tk.Label(self.root, text='Step Size (ss):').grid(row=16, column=0, columnspan=2, sticky='e')
        self.ss = tk.DoubleVar()
        ss_spin = tk.Spinbox(self.root,
                             values=(0.5, 1, 2, 2.5),
                             textvariable=str(self.ss),
                             width=5)
        self.ss.set(1)
        ss_spin.grid(row=16, column=2, stick='w')
        ss_spin.bind('<Return>', self.change_val)
        tk.Label(self.root, text='mm').grid(row=16, column=2, sticky='s')

        tk.Label(self.root, text='').grid(row=17)

        submit_button = tk.Button(self.root,
                                  text='Submit Values',
                                  command=self.submit_click,
                                  width=15,
                                  pady=5)
        submit_button.grid(row=20, column=0, columnspan=2, sticky='s')

        cancel_button = tk.Button(self.root,
                                  text='Cancel',
                                  command=self.cancel,
                                  width=15,
                                  pady=5)
        cancel_button.grid(row=20, column=2, columnspan=2, sticky='s')

        self.root.bind('<Return>', self.submit_click)

        self.root.mainloop()

    def submit_click(self, *args):
        var = self.variables
        var['ref_shape'] = self.shape.get()
        var['freq'] = self.freq_val.get()
        var['ang_freq'] = 2 * self.freq_val.get() * np.pi
        var['c_air'] = self.c_air
        var['dens'] = self.dens
        var['wl'] = self.c_air / self.freq_val.get()
        var['k'] = 2 * np.pi / var['wl']
        var['h'] = self.height.get() / 1000
        var['R'] = self.radius.get() / 1000
        var['c_xy'] = (0, round(var['h'] - var['R'], 6))
        var['trans_dia'] = self.trans_dia.get() / 1000
        var['ref_w'] = self.ref_w.get() / 1000
        var['ss'] = self.ss.get() / 1000

        for k, v in var.items():
            if isinstance(v, str) or isinstance(v, tuple):
                continue
            else:
                var[k] = round(v, 6)

        self.root.destroy()
        return var

    def cancel(self):
        var = None
        self.root.destroy()
        return var

    def change_val(self, event):
        self.ang_freq_str.set(str(round(self.freq_val.get() * 2 * np.pi, 3)) + ' rad/s')
        self.wl_str.set(str(round((self.c_air / self.freq_val.get()) * 1000, 3)) + ' mm')
        self.k_str.set(str(round((2 * np.pi / (self.c_air / self.freq_val.get())) / 1000, 3)) + ' rad/mm')
