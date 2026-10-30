#!/usr/bin/env python3
# ATES doublet. Seasons: well 1 injects 17 C, rest, well 2 injects 9 C, rest, well 1
# injects again (fixed mesh: five 6-day seasons; adaptive: three 3-day seasons, events_*.csv).
# Checks: temperatures within 9-17 C, both wells move the same volume, and the mean
# production temperatures (warm well in season 3, cold well in season 5) match
# reference_<case>.txt (written on the first run, commit it).
import csv
import glob
import os
import re

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

BIN = os.getcwd()[: os.getcwd().index("ICFERST")] + "bin/icferst"
# name: (tolerance on the means [C], season length [s], number of seasons)
CASES = {"ates_fixed": (1e-3, 518400.0, 5), "ates_adapt": (1e-2, 259200.0, 3)}


def run(name):
    for f in glob.glob(name + "_*"):
        os.remove(f)
    os.system(f"{BIN} {name}.mpml > log_{name}.txt 2>&1")
    rows = []
    if os.path.exists(name + "_outfluxes.csv"):
        rd = csv.reader(open(name + "_outfluxes.csv"))
        head = [h.strip() for h in next(rd)]
        rows = [dict(zip(head, map(float, r))) for r in rd if len(r) == len(head)]
    vtus = glob.glob(name + "_[0-9]*.vtu")
    vtus.sort(key=lambda f: int(re.findall(r"\d+", f)[-1]))
    return rows, vtus


def temperatures(f):
    r = vtk.vtkXMLUnstructuredGridReader()
    r.SetFileName(f)
    r.Update()
    pd = r.GetOutput().GetPointData()
    x = vtk_to_numpy(r.GetOutput().GetPoints().GetData())
    d = {
        pd.GetArrayName(i): vtk_to_numpy(pd.GetArray(i)).ravel()
        for i in range(pd.GetNumberOfArrays())
    }
    pipe = (np.abs(np.abs(x[:, 0]) - 30.0) < 1e-3) & (np.abs(x[:, 1]) < 1e-3)
    pipe &= x[:, 2] > -160.001
    return np.append(d["aquifer::Temperature"], d["well::Temperature"][pipe])


def production_mean(rows, col, season, P):
    t = np.array([r["Time[s]"] for r in rows])
    v = np.array([r[col] for r in rows])
    dt = np.diff(np.concatenate(([0.0], t)))
    sel = (t > (season - 1) * P) & (t <= season * P + 1.0)
    return float(np.sum(v[sel] * dt[sel]) / np.sum(dt[sel]))


def check(name, rows, vtus, tol, P, nseason):
    T = np.concatenate([temperatures(f) for f in vtus])
    print(f"   T in [{T.min():.4f}, {T.max():.4f}] C")
    ok = T.min() >= 9.0 - 0.05 and T.max() <= 17.0 + 0.05
    v1 = np.array([r["well[S22]CumuVolume"] for r in rows])
    v2 = np.array([r["well[S21]CumuVolume"] for r in rows])
    imbalance = np.max(np.abs(v1 + v2)) / np.max(np.abs(v1))
    print(f"   well volume imbalance {imbalance:.1e}")
    ok = ok and imbalance <= 1e-3
    vals = {
        "warm": production_mean(rows, "well[S22]AVG(Temperature)", 3, P),
    }
    if nseason >= 5:
        vals["cold"] = production_mean(rows, "well[S21]AVG(Temperature)", 5, P)
    ok = ok and 13.0 < vals["warm"] < 17.0 and 9.0 < vals.get("cold", 11.0) < 13.0
    ref = "reference_" + name + ".txt"
    if not os.path.exists(ref):
        open(ref, "w").writelines(f"{k} {v:.10e}\n" for k, v in vals.items())
    for k, v in (line.split() for line in open(ref)):
        print(f"   {k} well production {vals[k]:.6f} C, reference {float(v):.6f} C")
        ok = ok and abs(vals[k] - float(v)) <= tol
    return ok


Passed = True
for name, (tol, P, nseason) in CASES.items():
    print("== " + name)
    rows, vtus = run(name)
    if not rows or rows[-1]["Time[s]"] < nseason * P - 1.0:
        print("".join(open("log_" + name + ".txt").readlines()[-15:]))
        ok = False
    else:
        ok = check(name, rows, vtus, tol, P, nseason)
    print("   " + ("ok" if ok else "FAILED"))
    Passed = Passed and ok
print("ATES doublet: " + ("OK" if Passed else "FAILED"))
