#!/usr/bin/env python

from reader.gsr import *
from include.option_parser import *
from include.plot_data import *

if __name__ == '__main__':
    op = OptionParser().get_options()

    Udist = 1.23427e+20  # In cm
    Uvel = 1.5646e+07  # In cm/s
    Umasa = 9.9455e+42  # In grams
    Utime = Udist/Uvel

    GRAVITY = 6.672e-8
    BOLTZMANN = 1.3806e-16
    PROTONMASS = 1.6726e-24
    Msun = 1.989e33
    parsec = 3.08567758e18
    year = 31556926
    THOMPSON = 6.65245e-25
    C = 2.9979e10

    G = GRAVITY/(Udist**3) * Umasa * (Utime**2)
    Xh = 0.76  # Mass fraction of hydrogen
    mu = 4.0/(3 * Xh + 1) * PROTONMASS  # Neutral gas
    gamma = 5.0/3.0
    Uenergy = Umasa * (Udist**2) / (Utime**2)

    number_of_files = len(op.fname)
    times = np.zeros(number_of_files)
    bh_mass1 = np.zeros(number_of_files)
    bh_mass2 = np.zeros(number_of_files)

    for index, f in enumerate(op.fname):
        print("-> Processing snapshot", f,"...")
        snap = Snapshot(f)
        Header = snap.get_header()
        times[index] = Header['Time']

        bh_data = snap.get_data_by_type(5)
        bh_mass1[index], bh_mass2[index] = bh_data[1]  # Mass

    # Now we have the time of the snapshot in `times`
    # and the masses of the BH in `bh_mass`


    Medd_constant = 4 * np.pi * GRAVITY * C * PROTONMASS / (0.1 * C * C * THOMPSON)
    print("Medd_constant: ", Medd_constant)

    bh_mass1 *= Umasa/Msun
    bh_mass2 *= Umasa/Msun
    times *= Utime/year

    times_frac = np.zeros(number_of_files -1)
    mdot1_medd = np.zeros(number_of_files -1)
    mdot2_medd = np.zeros(number_of_files -1)

    for i in range(number_of_files-1):

        dm1 = bh_mass1[i+1] - bh_mass1[i]
        dm2 = bh_mass2[i+1] - bh_mass2[i]
        dt = times[i+1] - times[i]

        Medd1 = (bh_mass1[i]*Msun)*Medd_constant
        Medd1 *= year/Msun

        Medd2 = (bh_mass2[i]*Msun)*Medd_constant
        Medd2 *= year/Msun

        Mdot1 = dm1/dt
        Mdot2 = dm2/dt

        times_frac[i] = (times[i+1] + times[i])/2.0
        mdot1_medd[i] = Mdot1/Medd1
        mdot2_medd[i] = Mdot2/Medd2

    print(times_frac)
    print(mdot1_medd)
    print(mdot2_medd)
