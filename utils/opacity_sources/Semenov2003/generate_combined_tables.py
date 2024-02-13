"""
generate_combined_tables.py

Authors: Alessandro A. Trani, Daria Gangardt

This file uses the Semenov2003 opacity Python interface and the Badnell2005 table
to generate a combined table of opacities (Semenov2003 at low T, Badnell2005 at high T),
following the table structure of Badnell2005

Check pysemenov_module.sh for instructions how to compile the Semenov2003 opacity interface,
and check test_pysemenov.py for instructions on how to use the interface
"""
import numpy as np
import pysemenov


def get_badnel2005_tables(fname="../AGS05.OP17"):
    """ This function parses the Badnel2005 tables and returns the table following table:
    TABLE # 74     200904130004       X=0.7000 Y=0.2700 Z=0.0300 dX1=0.0000 dX2=0.0000
    x-axis is log10(T), temperature in kelvin,
    y-axis is log10(R), where R=density in g/cm**3 / T6**3 and T6 = 1.e-6 * T in kelvin

    Parameters
    ----------
    fname: str, optional (default: "../AGS05.OP17")
        name and path to the Badnell2005 table file

    Returns
    ----------
    TR_matrix: 2d ndarray
        Opacity table, in log10(kappa[cm**2/g])
    logR_col: 1d ndarray
        Array of the x-axis in log10(T[kelvin]), 70 values from 3.75 to 8.70
    logT_row: 1d ndarray
        Array of the y-axis in log10(R), R=rho[g/cm**3]/T6**3, where T6=1.e-6*T[kelvin], 19 values from -8.0 to +1.0
    """
    f = open(fname, "r")
    istart = 5792
    istop = 5859
    TR_matrix = np.zeros(shape=(66,18), dtype=np.float64)
    logT_row = []
    k = 0
    for i, l in enumerate(f):
        if i == istart-1:
            coll = l.strip().split()
            #print(coll)
            logR_col = np.array(coll[2:], dtype=np.float64)

        if i>istart and i<istop:
            coll = l.strip().split()
            #print(coll)
            logT_row.append(coll[0])
            TR_matrix[k,:] = np.array(coll[2:], dtype=np.float64)
            k+=1

    logT_row = np.array(logT_row, dtype=np.float64)
    return TR_matrix, logR_col, logT_row


def generate_combined_table():
    """ This function combines the Badnell2005 tables (at high temperature) with the Semenov2003 tables (at low temperature)
    Everything in still in cgs units
    Saves the tables, the x-axis, and y-axis as logRT_combined_table.dat, logRT_combined_Taxis.dat and logRT_combined_Raxis.dat
    The table follows the Badnell2005 convention, so the table is in log10(kappa[cm**2/g]), the x-axis is in log10(T[kelvin])
    and the y-axis is in log10(R), R=rho[g/cm**3]/T6**3, where T6=1.e-6*T[kelvin],
    """
    TRmat, logR, logT = get_badnel2005_tables()

    print("Badnel logT range", logT)
    print("Badnel logR range", logR)
    print("Badnel min T:", logT.min(), 10**logT.min(), "log10[kelvin], [kelvin]")

    # Semenov has limits in 5[K]<T<~10,000[K] and for gas
    # rho = 2*10^-18 [g/cm^3] and ~2*10^-7[g/cm^3].

    # Let's use Semenov for T<1e4
    # So we cut the first 5 columns of the Badnell2005 table
    Tcut = np.log10(7e3)   # This seems the best value that ensures the smoothest transition
    iTcut = np.argwhere(logT>=Tcut)[0][0]
    print("logTcut[kelvin]", Tcut)
    print("fist temperature column index included:", iTcut)
    firstThi = logT[iTcut]
    print("first temperature column value included:", firstThi, 10**firstThi, "log10[kelvin], [kelvin]")

    # We keep only the high temperature from Badnell2005
    logT_hiT = logT[iTcut:]
    TRmat_hiT = TRmat[iTcut:,:]
    print("highT x-axis:", logT_hiT)

    # Continue the T, R table at low T
    # Then we decide the temperature range for T<1e4
    Tmin = 5  # kelvin
    logT_loT = np.linspace(np.log10(Tmin), firstThi*0.95, 100)
    print("lowT x-axis:", logT_loT)

    # Options for the Semenov2003 opacities
    ross = True  # False for Planck opacities
    model = 'nrm'  # Normal iron abundances
    top = 'c'  # composite grains
    shap = 'a'  # aggregate shape

    TRmat_loT = np.zeros(shape=(len(logT_loT), len(logR)), dtype=np.float64)
    for iT, lT in enumerate(logT_loT):
        logrho_loT = logR + lT * 3 - 18  # Badnell2005 relation between R, rho and T in logarithmic form
        for irho, lrho in enumerate(logrho_loT):
            rho = 10**lrho
            T = 10**lT
            kappaS = pysemenov.compute_kappa(top=top, ross=ross, model=model, shap=shap, rho=rho, T=T)
            logK = 9.999 if kappaS == 0.0 else np.log10(kappaS)
            TRmat_loT[iT,irho] = logK

    print(TRmat_loT.shape)
    print(TRmat_hiT.shape)

    # Here we glue the tables together
    TRmat_full = np.vstack((TRmat_loT,TRmat_hiT))
    logT_full = np.hstack((logT_loT, logT_hiT))
    logR_full = logR

    np.savetxt("../logRT_combined_table.dat", TRmat_full)
    np.savetxt("../logRT_combined_Raxis.dat", logR_full)
    np.savetxt("../logRT_combined_Taxis.dat", logT_full)


if __name__ == "__main__":
    generate_combined_table()




