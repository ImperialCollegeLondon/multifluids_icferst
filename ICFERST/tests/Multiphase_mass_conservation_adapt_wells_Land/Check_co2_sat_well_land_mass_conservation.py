#!/usr/bin/env python3

# REGRESSION TEST: per-phase SATURATION mass conservation across mesh adapt
#
# Physical background: the Galerkin projection conserves int(S)dV, not the pore
# volume int(phi*S)dV. With a porosity contrast an unweighted
# projection moves phase volume across the contrast at every adapt. The chain
# under test: phi-scaling of PhaseVolumeFraction around the projection +
# bounded::Diffuse in the projection + BoundedSolutionCorrections +
# Set_Saturation_to_sum_one + Restore_Phase_Mass_SumOne.
#
# The model is a closed box, ZERO FLOW, gravity off, with an S=0.5 CO2 blob
# (centre (100,100,80), r=35) crossing both porosity interfaces. With zero flow
# ANY change of int(phi*V*S) per phase between dumps is the adapt chain alone.
#
# 15 fixed 1 s timesteps, adapt every step: ~15 adapt events, target < 1 min.

import vtk
import sys
import os
import math
import numpy as np

# TOLERANCES
# Measured with the fix chain: <= ~1e-11 relative per adapt event (toys).
# Fail thresholds 2-3 orders above measured, still 5+ orders below the pre-fix
# defect (O(1e-3..1e-2) relative per adapt on this contrast).
Tolerance_per_event = 1.0e-6    # per phase: |d int(phi V S)| / M0_phase per dump
Tolerance_total     = 1.0e-6    # per phase: |M_end - M_0| / M0_phase
Tolerance_sum_one   = 1.0e-8    # max |1 - sum_p S_p| over nodes, every dump
Bounds_slack        = 1.0e-12   # S must stay in [0-eps, 1+eps]
Far_radius          = 75.0      # [m] from the blob centre. Blob r=35 + ~2
                                # elements of the COARSE regression mesh
                                # (~20-25 m): one element of legitimate
                                # interface smearing must stay inside, or the
                                # guard fires
Far_S_max           = 1.0e-2    # coarse placement guard
Blob_centre         = np.array([100.0, 100.0, 80.0])
Min_dumps           = 5
Min_mesh_changes    = 3
Phases              = ["water", "CO2"]   # reservoir phases to check
Land_phase          = "CO2"    # reservoir phase carrying the Land immobile frac
Imm_drop_max        = 0.40     # max relative drop of max CV_ImmobileFraction per
                               # dump AND first->last. With the fix the immobile
                               # fraction holds; without it, it collapses by ~the
                               # porosity factor (~0.8) at the first adapt.
Imm_floor           = 1.0e-3   # a dump whose max immobile is below this counts
                               # as collapsed regardless of ratio

print('Running the model')

path = os.getcwd()
binpath = path[:path.rindex('ICFERST')] + 'bin/icferst'
os.system('rm -f ' + path + '/*.vtu')
os.system(binpath + ' ' + path + '/co2_sat_well_land_mass_conservation_adapt.mpml')

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


def get_array(container, name):
    for i in range(container.GetNumberOfArrays()):
        n = container.GetArrayName(i)
        if n and (n == name or n.endswith('::' + name) or n == name.replace('::', '')):
            return container.GetArray(i)
    return None


def get_phase_array(container, phase, suffix):
    want = phase + '::' + suffix
    for i in range(container.GetNumberOfArrays()):
        n = container.GetArrayName(i)
        if n == want:
            return container.GetArray(i)
    return None


def tet_cells(data):
    try:
        types = vtk_np(data.GetCellTypesArray())
    except AttributeError:
        types = np.array([data.GetCellType(c)
                          for c in range(data.GetNumberOfCells())])
    ids = np.nonzero(types == vtk.VTK_TETRA)[0]
    try:
        conn = vtk_np(data.GetCells().GetConnectivityArray())
        offs = vtk_np(data.GetCells().GetOffsetsArray())
        starts = offs[ids]
        tconn = np.stack([conn[starts + k] for k in range(4)], axis=1)
    except AttributeError:
        tconn = np.empty((len(ids), 4), dtype=np.int64)
        idl = vtk.vtkIdList()
        for j, c in enumerate(ids):
            data.GetCellPoints(int(c), idl)
            tconn[j] = [idl.GetId(k) for k in range(4)]
    return tconn, ids


def measure(fname):
    """Per-phase sum(phi*V*<S>_e), max|1-sum S|, S range, far-field max S,
    node count -- one dump."""
    data = read_vtu(fname)
    pts = vtk_np(data.GetPoints().GetData()).astype(np.float64)
    tconn, tids = tet_cells(data)
    por_arr = get_array(data.GetCellData(), 'Porosity')
    por = (vtk_np(por_arr)[tids].astype(np.float64)
           if por_arr is not None else np.ones(len(tids)))
    a = pts[tconn[:, 0]]
    vol = np.abs(np.einsum('ij,ij->i',
                           pts[tconn[:, 1]] - a,
                           np.cross(pts[tconn[:, 2]] - a,
                                    pts[tconn[:, 3]] - a))) / 6.0
    Ssum = None
    out = {}
    far = 0.0
    smin, smax = 1e30, -1e30
    d2 = ((pts - Blob_centre[None, :]) ** 2).sum(axis=1)
    for p in Phases:
        arr = get_phase_array(data.GetPointData(), p, 'PhaseVolumeFraction')
        if arr is None:
            raise RuntimeError('%s::PhaseVolumeFraction not found in %s'
                               % (p, fname))
        S = vtk_np(arr).astype(np.float64).ravel()
        out[p] = math.fsum(np.asarray(vol * por * S[tconn].mean(axis=1),
                                      dtype=np.float64))
        Ssum = S if Ssum is None else Ssum + S
        smin = min(smin, float(S.min()))
        smax = max(smax, float(S.max()))
        if p == 'CO2':
            m = d2 > Far_radius ** 2
            if m.any():
                far = float(np.abs(S[m]).max())
    sum_res = float(np.abs(Ssum - 1.0).max())
    imm_arr = get_phase_array(data.GetPointData(), Land_phase, 'CV_ImmobileFraction')
    imm_max = (float(np.abs(vtk_np(imm_arr)).max()) if imm_arr is not None else None)
    return out, sum_res, smin, smax, far, imm_max, data.GetNumberOfPoints()


masses = {p: [] for p in Phases}
nnodes = []
imm_hist = []
worst_sum = 0.0
worst_far = 0.0
smin_all, smax_all = 1e30, -1e30
if Passed:
    for i in range(AutoNumber + 1):
        out, sres, smin, smax, far, imm_max, n = measure(
            os.path.join(path, '%s_%d.vtu' % (AutoFile, i)))
        for p in Phases:
            masses[p].append(out[p])
        imm_hist.append(imm_max)
        nnodes.append(n)
        worst_sum = max(worst_sum, sres)
        worst_far = max(worst_far, far)
        smin_all = min(smin_all, smin)
        smax_all = max(smax_all, smax)
        print('  dump %3d: nodes = %7d, ' % (i, n) +
              ', '.join('phi*V*S[%s] = %.10e' % (p, out[p]) for p in Phases) +
              ', max|1-sum| = %.2e' % sres)

if Passed:
    for p in Phases:
        M0 = masses[p][0]
        if abs(M0) < 1e-30:
            print('Initial %s volume is zero: broken initial condition' % p)
            Passed = False
            continue
        worst = max((abs(masses[p][i] - masses[p][i - 1]) / M0
                     for i in range(1, len(masses[p]))), default=0.0)
        total = abs(masses[p][-1] - M0) / M0
        print('%s: worst per-event drift %.3e (tol %.1e), total %.3e (tol %.1e)'
              % (p, worst, Tolerance_per_event, total, Tolerance_total))
        if worst > Tolerance_per_event or total > Tolerance_total:
            Passed = False

    changes = sum(1 for i in range(1, len(nnodes)) if nnodes[i] != nnodes[i - 1])
    print('max |1 - sum_p S_p|          : %.3e (tolerance %.1e)'
          % (worst_sum, Tolerance_sum_one))
    print('saturation range             : [%.3e, %.3e]' % (smin_all, smax_all))
    print('max far-field CO2 S (>%.0f m) : %.3e (tolerance %.1e)'
          % (Far_radius, worst_far, Far_S_max))
    print('dumps with changed node count: %d (require >= %d)'
          % (changes, Min_mesh_changes))
    if worst_sum > Tolerance_sum_one:
        Passed = False
    if smin_all < -Bounds_slack or smax_all > 1.0 + Bounds_slack:
        Passed = False
    if worst_far > Far_S_max:
        Passed = False
    if changes < Min_mesh_changes:
        print('Mesh barely changed: adaptation did not exercise the projection')
        Passed = False

# NOTE: the well-block saturation (Well_H2O / Well_CO2) is NOT checked. In this
# toy those prognostic phases collapse to zero on the first ordinary timestep,
# before any adapt, identically in stock IC-FERST and every campaign build : it is
# pre-existing well-solver behaviour, unrelated to the mass fix. This test guards
# reservoir mass conservation and the Land immobile fraction (below).

# ---- Land verdict: the CO2 immobile fraction must not collapse across adapts.
# With fix F1 (deferred flip/cap after the unscale) it holds; without it the cap
# runs on the scaled phi*S inside the adapt window and it collapses.
# Run this independently of the mass verdict so it always prints (the mass and
# Land verdicts are separate concerns; a marginal mass reading must not hide Land).
if imm_hist:
    if any(v is None for v in imm_hist):
        print('CV_ImmobileFraction not found: add the diagnostic to the CO2 '
              'phase, or the binary does not populate it')
        Passed = False
    else:
        imm0 = imm_hist[0]
        if imm0 < Imm_floor:
            print('Initial CO2 immobile fraction %.3e below floor %.1e: the '
                  'Land_coefficient is not producing trapping' % (imm0, Imm_floor))
            Passed = False
        else:
            worst_event = max(((imm_hist[i - 1] - imm_hist[i]) / imm_hist[i - 1]
                               for i in range(1, len(imm_hist))
                               if imm_hist[i - 1] > Imm_floor), default=0.0)
            total_drop = (imm0 - imm_hist[-1]) / imm0
            collapsed = any(v < Imm_floor for v in imm_hist)
            print('CO2 immobile fraction: first %.6e last %.6e' % (imm0, imm_hist[-1]))
            print('  worst per-event drop %.3e, total drop %.3e (tol %.2f)'
                  % (worst_event, total_drop, Imm_drop_max))
            if worst_event > Imm_drop_max or total_drop > Imm_drop_max or collapsed:
                Passed = False

if Passed:
    print('Well saturation mass and Land trapping across adaptivity works OK')
else:
    print('Well saturation mass and Land trapping across adaptivity does NOT work')
