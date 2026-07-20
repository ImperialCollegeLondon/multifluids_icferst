#!/usr/bin/env python3

# REGRESSION TEST: tracer mass conservation across mesh adapt.
#
# Physical background: the Galerkin mesh-to-mesh projection used during
# adaptation conserves the volume integral (int C dV), not physical tracer
# mass (int phi*C dV). With a porosity contrast (here 0.5 / 0.06), an
# unweighted projection moves physical mass across the contrast at every
# adapt. 
#
#   Config: Impose_min_max OFF, bounded::Diffuse ON in galerkin_projection
#
# The model is a closed box, ZERO FLOW, with a sharp tracer
# blob through the porosity contrast.
# With zero flow, transport does nothing: ANY change in
# sum(phi * V * C) between dumps is caused by the adapt chain alone.

import vtk
import sys
import os
import math
import numpy as np

# TOLERANCE OF THE CHECKING
# Measured with the fixes: ~1e-13 relative drift per adapt event.
# Fail threshold set 4 orders above measured, still ~4 orders below
# the pre-fix defect (which lost O(1e-4..1e-2) relative per adapt).
Tolerance_per_event = 1.0e-9   # max |d(phi*V*C)| / M0 between consecutive dumps
Tolerance_total     = 1.0e-9   # |M_end - M_0| / M0
Min_dumps           = 5        # the run must produce at least this many dumps
Min_mesh_changes    = 3        # at least this many dumps must show a changed
                               # node count (i.e. adapt actually happened)

print('Running the model')

# Get path and run IC-FERST
path = os.getcwd()
binpath = path[:path.rindex('ICFERST')] + 'bin/icferst'
os.system('rm -f ' + path + '/*.vtu')
os.system(binpath + ' ' + path + '/tracer_mass_conservation_adapt.mpml')

# ---------------- locate the dump series ----------------
AutoNumber = -1
AutoFile = ''
for files in os.listdir(path):
    if files.endswith(".vtu") and 'checkpoint' not in files.lower():
        pos = files.rfind('_')
        pos2 = files.rfind('.')
        tail = files[pos + 1:pos2]
        if not tail.isdigit():
            continue
        AutoFile = files[:pos]
        AutoNumber = max(AutoNumber, int(tail))

Passed = True
if AutoNumber < Min_dumps:
    print('Only %d dumps found (need >= %d): the run did not complete'
          % (AutoNumber + 1, Min_dumps))
    Passed = False


def read_vtu(fname):
    r = vtk.vtkXMLUnstructuredGridReader()
    r.SetFileName(fname)
    r.Update()
    return r.GetOutput()


def vtk_np(arr):
    from vtk.util import numpy_support
    return None if arr is None else numpy_support.vtk_to_numpy(arr)


def get_array(container, suffix):
    """Return the first array whose name ends with `suffix`
    (handles any phase prefix, e.g. 'aquifer::Concentration')."""
    for i in range(container.GetNumberOfArrays()):
        n = container.GetArrayName(i)
        if n and (n == suffix or n.endswith('::' + suffix)):
            return container.GetArray(i)
    return None


def tet_cells(data):
    """Connectivity (ntets, 4) of tetrahedral cells + their cell ids."""
    try:
        types = vtk_np(data.GetCellTypesArray())
    except AttributeError:  # very old/new VTK API differences
        types = np.array([data.GetCellType(c)
                          for c in range(data.GetNumberOfCells())])
    ids = np.nonzero(types == vtk.VTK_TETRA)[0]
    try:  # VTK >= 9 fast path
        conn = vtk_np(data.GetCells().GetConnectivityArray())
        offs = vtk_np(data.GetCells().GetOffsetsArray())
        starts = offs[ids]
        tconn = np.stack([conn[starts + k] for k in range(4)], axis=1)
    except AttributeError:  # older VTK: per-cell fallback
        tconn = np.empty((len(ids), 4), dtype=np.int64)
        idl = vtk.vtkIdList()
        for j, c in enumerate(ids):
            data.GetCellPoints(int(c), idl)
            tconn[j] = [idl.GetId(k) for k in range(4)]
    return tconn, ids


def phi_v_c(fname):
    """Return (sum phi*V*C, n_nodes) for one dump.
    Cell C = average of the 4 nodal values (PointDataToCellData convention);
    porosity is the P0 cell array; V from coordinates (exact tet volume)."""
    data = read_vtu(fname)
    pts = vtk_np(data.GetPoints().GetData()).astype(np.float64)
    conc = get_array(data.GetPointData(), 'Concentration')
    if conc is None:
        raise RuntimeError('Concentration point array not found in ' + fname)
    conc = vtk_np(conc).astype(np.float64)
    tconn, tids = tet_cells(data)
    por_arr = get_array(data.GetCellData(), 'Porosity')
    por = (vtk_np(por_arr)[tids].astype(np.float64)
           if por_arr is not None else np.ones(len(tids)))
    a = pts[tconn[:, 0]]
    vol = np.abs(np.einsum('ij,ij->i',
                           pts[tconn[:, 1]] - a,
                           np.cross(pts[tconn[:, 2]] - a,
                                    pts[tconn[:, 3]] - a))) / 6.0
    c_cell = conc[tconn].mean(axis=1)
    mass = math.fsum(np.asarray(vol * por * c_cell, dtype=np.float64))
    return mass, data.GetNumberOfPoints()


# ---------------- compute the series ----------------
masses = []
nnodes = []
if Passed:
    for i in range(AutoNumber + 1):
        m, n = phi_v_c(os.path.join(path, '%s_%d.vtu' % (AutoFile, i)))
        masses.append(m)
        nnodes.append(n)
        print('  dump %3d: nodes = %7d, phi*V*C = %.10e' % (i, n, m))

    M0 = masses[0]
    if abs(M0) < 1e-30:
        print('Initial mass is zero: broken initial condition')
        Passed = False

if Passed:
    worst = 0.0
    for i in range(1, len(masses)):
        rel = abs(masses[i] - masses[i - 1]) / M0
        worst = max(worst, rel)
    total = abs(masses[-1] - M0) / M0
    changes = sum(1 for i in range(1, len(nnodes)) if nnodes[i] != nnodes[i - 1])

    print('worst per-event relative drift : %.3e (tolerance %.1e)'
          % (worst, Tolerance_per_event))
    print('total relative drift           : %.3e (tolerance %.1e)'
          % (total, Tolerance_total))
    print('dumps with changed node count  : %d (require >= %d)'
          % (changes, Min_mesh_changes))

    if worst > Tolerance_per_event:
        Passed = False
    if total > Tolerance_total:
        Passed = False
    # Guard against the test silently passing because adapt never ran:
    # with no adapts, conservation is trivial and nothing is being tested.
    if changes < Min_mesh_changes:
        print('Mesh barely changed: adaptation did not exercise the projection')
        Passed = False

if Passed:
    print('Tracer mass conservation across adaptivity works OK')
else:
    print('Tracer mass conservation across adaptivity does NOT work')
