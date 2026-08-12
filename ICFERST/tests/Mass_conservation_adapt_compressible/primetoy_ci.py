#!/usr/bin/env python3
"""Extraction for the test case, for all three mpmls.

For each dump family it computes, per tracked field, the worst per-adapt relative jumps of
the exact P1 integral (dF) and of the lumped rho*F product (lump), using exact P1 quadrature
and a fresh vtk reader per file.
"""
import glob, re
import numpy as np
if not hasattr(np, "bool"):
    np.bool = bool
import vtk
from vtk.util.numpy_support import vtk_to_numpy

FIELDS = ['Concentration', 'Temperature']


def _tet_connectivity(g):
    ncells = g.GetNumberOfCells()
    ca = g.GetCells()
    if hasattr(ca, 'GetConnectivityArray'):
        return vtk_to_numpy(ca.GetConnectivityArray()).reshape(ncells, 4)
    conn = np.empty((ncells, 4), dtype=np.int64)
    for i in range(ncells):
        ids = g.GetCell(i).GetPointIds()
        conn[i] = [ids.GetId(0), ids.GetId(1), ids.GetId(2), ids.GetId(3)]
    return conn


def _read(path):
    r = vtk.vtkXMLUnstructuredGridReader()
    r.SetFileName(path)
    r.Update()
    g = r.GetOutput()
    pd = g.GetPointData()
    names = [pd.GetArrayName(i) for i in range(pd.GetNumberOfArrays())]
    pts = vtk_to_numpy(g.GetPoints().GetData()).astype(float)
    conn = _tet_connectivity(g)
    a = pts[conn[:, 0]]
    vol = np.abs(np.einsum('ij,ij->i',
                           np.cross(pts[conn[:, 1]] - a, pts[conn[:, 2]] - a),
                           pts[conn[:, 3]] - a)) / 6.0
    dname = [n for n in names if n.endswith('::Density') or n == 'Density'][0]
    Rn = vtk_to_numpy(pd.GetArray(dname)).astype(float).ravel()[conn]
    out = {}
    for f in FIELDS:
        arr = [n for n in names if n.endswith('::' + f) or n == f]
        if not arr:
            continue
        Fn = vtk_to_numpy(pd.GetArray(arr[0])).astype(float).ravel()[conn]
        out[f] = (np.sum(vol / 4.0 * Fn.sum(axis=1)),
                  np.sum(vol / 4.0 * (Rn * Fn).sum(axis=1)))
    return len(pts), out


def family(prefix):
    vtus = sorted(glob.glob(prefix + '_[0-9]*.vtu'),
                  key=lambda f: int(re.search(r'_(\d+)\.vtu$', f).group(1)))
    worst = {}
    prev = None
    for fp in vtus:
        n, vals = _read(fp)
        if prev is not None and n != prev[0]:
            for f in vals:
                pF, pL = prev[1][f]
                iF, iL = vals[f]
                w = worst.setdefault(f, {'dF': 0.0, 'lump': 0.0})
                w['dF'] = max(w['dF'], abs((iF - pF) / pF))
                w['lump'] = max(w['lump'], abs((iL - pL) / pL))
        prev = (n, vals)
    return worst


def all_ledgers():
    return {'rhoC': family('z_primetoy_box_compressible_rhoC'),
            'rhoCT': family('z_primetoy_box_compressible_rhoCT'),
            'rhoC_disabled': family('z_primetoy_box_compressible_rhoC_disabled')}


if __name__ == '__main__':
    for k, v in all_ledgers().items():
        print(k, v)
