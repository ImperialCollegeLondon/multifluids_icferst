#!/usr/bin/env python3
# Well doublet, compressible brine. uniform_*: injecting the aquifer T and C must leave
# them uniform. doublet_*: injecting hot salty brine must stay bounded, not reach the
# producer, and match reference_<case>.txt (written on the first run, commit it).
# doublet_beta0_incomp_adapt also carries a PassiveTracer, checked like C.
import csv
import glob
import os
import re

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

T0, C0, TI, CI, T_END = 324.25, 0.15, 423.15, 0.20, 345600.0
BIN = os.getcwd()[: os.getcwd().index("ICFERST")] + "bin/icferst"
# name: (tolerance on T, tolerance on C)
UNIFORM = {"uniform_beta1_fixed": (1e-7, 1e-10), "uniform_beta1_adapt": (1e-6, 1e-10)}
DOUBLET = {
    "doublet_beta1_fixed": (1e-3, 1e-7),
    "doublet_beta1_adapt": (2.0, 1e-6),
    "doublet_beta0_fixed": (1e-3, 1e-7),
    "doublet_beta0_incomp_adapt": (2.0, 1e-6),
}


def run(name):
    for f in glob.glob(name + "_*"):
        os.remove(f)
    os.system(f"{BIN} {name}.mpml > log_{name}.txt 2>&1")
    rows = []
    if os.path.exists(name + "_outfluxes.csv"):
        rows = [r for r in csv.reader(open(name + "_outfluxes.csv")) if r]
    vtus = glob.glob(name + "_[0-9]*.vtu")
    vtus.sort(key=lambda f: int(re.findall(r"\d+", f)[-1]))
    return len(rows) > 1 and float(rows[-1][0]) >= T_END - 1.0, vtus


def read(f):
    r = vtk.vtkXMLUnstructuredGridReader()
    r.SetFileName(f)
    r.Update()
    g, pd = r.GetOutput(), r.GetOutput().GetPointData()
    d = {
        pd.GetArrayName(i): vtk_to_numpy(pd.GetArray(i)).ravel()
        for i in range(pd.GetNumberOfArrays())
    }
    x = vtk_to_numpy(g.GetPoints().GetData())
    cells = g.GetCells()
    if hasattr(cells, "GetConnectivityArray"):
        conn = vtk_to_numpy(cells.GetConnectivityArray()).reshape(-1, 4)
    else:
        conn = vtk_to_numpy(cells.GetData()).reshape(-1, 5)[:, 1:]
    p = g.GetCellData().GetArray("phase1::Porosity")
    phi = vtk_to_numpy(p).ravel() if p else np.ones(len(conn))
    pipe = (np.abs(np.abs(x[:, 0]) - 100.0) < 1e-3) & (np.abs(x[:, 1]) < 1e-3)
    pipe &= x[:, 2] > -200.001
    T = np.append(d["phase1::Temperature"], d["phase2::Temperature"][pipe])
    C = np.append(d["phase1::Concentration"], d["phase2::Concentration"][pipe])
    if "phase1::PassiveTracer" in d:
        P = np.append(d["phase1::PassiveTracer"], d["phase2::PassiveTracer"][pipe])
        return x, conn, phi, T, C, P
    return x, conn, phi, T, C, None


def excess(x, conn, phi, v, ref):
    a, b, c, e = (x[conn[:, k]] for k in range(4))
    vol = np.abs(np.einsum("ij,ij->i", b - a, np.cross(c - a, e - a))) / 6.0
    return float(np.sum(phi * vol * np.mean(v[conn] - ref, axis=1)))


def check_uniform(vtus, tolT, tolC):
    dT = max(np.max(np.abs(read(f)[3] - T0)) for f in vtus)
    dC = max(np.max(np.abs(read(f)[4] - C0)) for f in vtus)
    print(f"   max|T-T0| {dT:.2e} K   max|C-C0| {dC:.2e}")
    return dT <= tolT and dC <= tolC


def check_doublet(name, vtus, tolT, tolC):
    ok = True
    for f in vtus:
        x, conn, phi, T, C, P = read(f)
        n = len(x)
        print(f"   {f}: T {T.min():.6f} {T.max():.6f} C {C.min():.10f} {C.max():.10f}")
        ok = ok and T.min() >= T0 - tolT and T.max() <= TI + tolT
        ok = ok and C.min() >= C0 - tolC and C.max() <= CI + tolC
        prod = np.hypot(x[:, 0] - 100.0, x[:, 1]) < 30.0
        ok = ok and np.max(np.abs(C[:n][prod] - C0)) <= 1e-9
        ok = ok and np.max(np.abs(T[:n][prod] - T0)) <= 1e-5
        if P is not None:
            print(f"   {f}: PassiveTracer {P.min():.10f} {P.max():.10f}")
            ok = ok and P.min() >= C0 - tolC and P.max() <= CI + tolC
            ok = ok and np.max(np.abs(P[:n][prod] - C0)) <= 1e-9
    ok = ok and C[:n].max() > 0.5 * (C0 + CI)
    vals = {
        "heat": excess(x, conn, phi, T[:n], T0),
        "salt": excess(x, conn, phi, C[:n], C0),
    }
    if P is not None:
        vals["tracer"] = excess(x, conn, phi, P[:n], C0)
    ref = "reference_" + name + ".txt"
    if not os.path.exists(ref):
        open(ref, "w").writelines(f"{k} {v:.10e}\n" for k, v in vals.items())
    for k, v in (line.split() for line in open(ref)):
        print(f"   {k} excess {vals[k]:.6e}, reference {float(v):.6e}")
        ok = ok and abs(vals[k] - float(v)) <= 2e-2 * abs(float(v))
    return ok


Passed = True
for name, (tolT, tolC) in {**UNIFORM, **DOUBLET}.items():
    print("== " + name)
    finished, vtus = run(name)
    if not finished:
        print("".join(open("log_" + name + ".txt").readlines()[-15:]))
        ok = False
    elif name in UNIFORM:
        ok = check_uniform(vtus, tolT, tolC)
    else:
        ok = check_doublet(name, vtus, tolT, tolC)
    print("   " + ("ok" if ok else "FAILED"))
    Passed = Passed and ok
print("Well doublet consistency: " + ("OK" if Passed else "FAILED"))
