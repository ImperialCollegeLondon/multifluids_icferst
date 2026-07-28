#!/usr/bin/env python3

# REGRESSION TEST: metal mass conservation across mesh adapt, WITH a live
# dissolution reaction.
#
# The Galerkin mesh-to-mesh projection used during adaptation conserves the
# plain volume integral (int C dV), not physical mass. For a metal exchange
# reaction the conserved invariant is the TOTAL metal mass
#
#     I = int( phi * S * rho_f * C_fluid  +  (1 - phi) * rho_s * C_solid ) dV
#
#   Config: Impose_min_max OFF, bounded::Diffuse ON, zero flow, porosity
#   contrast through a metal blob, LIVE reaction, correction rescale disabled.

import vtk
import sys
import os
import math
import numpy as np

if not hasattr(np, 'bool'):
    np.bool = bool

# TOLERANCE OF THE CHECKING
# Fail threshold set ~4 orders above the
# measured post-fix drift, still far below the pre-fix defect (which moves
# O(1e-4..1e-2) relative of the solid part per adapt across a 0.5/0.06 contrast).
Tolerance_per_event = 1.0e-9   # max |d I| / I0 between consecutive dumps
Tolerance_total     = 1.0e-9   # |I_end - I_0| / I0
Min_dumps           = 5        # the run must produce at least this many dumps
Min_mesh_changes    = 3        # at least this many dumps must show a changed
                               # node count (i.e. adapt actually happened)

# Constant physical parameters from the mpml (single-phase, incompressible)
RHO_FLUID = 1000.0
RHO_SOLID = 2650.0

print('Running the model')

path = os.getcwd()
binpath = path[:path.rindex('ICFERST')] + 'bin/icferst'
os.system('rm -f ' + path + '/*.vtu')
os.system(binpath + ' ' + path + '/metal_mass_conservation_adapt.mpml')

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
    (handles any phase prefix, e.g. 'aquifer::Cu_solid')."""
    for i in range(container.GetNumberOfArrays()):
        n = container.GetArrayName(i)
        if n and (n == suffix or n.endswith('::' + suffix)):
            return container.GetArray(i)
    return None


def tet_cells(data):
    """Connectivity (ntets, 4) of tetrahedral cells + their cell ids."""
    try:
        types = vtk_np(data.GetCellTypesArray())
    except AttributeError:
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


def metal_invariant(fname):
    """Return (I, n_nodes) for one dump, where
       I = sum over tets of V * [ phi*rho_f*C_fluid + (1-phi)*rho_s*C_solid ],
    each cell field being the average of its 4 nodal values (the
    PointDataToCellData convention used by the tracer test)."""
    data = read_vtu(fname)
    pts = vtk_np(data.GetPoints().GetData()).astype(np.float64)

    c_fluid = get_array(data.GetPointData(), 'PassiveTracer_Cu')
    c_solid = get_array(data.GetPointData(), 'Cu_solid')
    if c_fluid is None:
        raise RuntimeError('PassiveTracer_Cu point array not found in ' + fname)
    if c_solid is None:
        raise RuntimeError('Cu_solid point array not found in ' + fname)
    c_fluid = vtk_np(c_fluid).astype(np.float64)
    c_solid = vtk_np(c_solid).astype(np.float64)

    tconn, tids = tet_cells(data)
    por_arr = get_array(data.GetCellData(), 'Porosity')
    por = (vtk_np(por_arr)[tids].astype(np.float64)
           if por_arr is not None else np.ones(len(tids)))

    a = pts[tconn[:, 0]]
    vol = np.abs(np.einsum('ij,ij->i',
                           pts[tconn[:, 1]] - a,
                           np.cross(pts[tconn[:, 2]] - a,
                                    pts[tconn[:, 3]] - a))) / 6.0

    cf_cell = c_fluid[tconn].mean(axis=1)
    cs_cell = c_solid[tconn].mean(axis=1)

    # single-phase incompressible: S = 1, rho_f constant
    integrand = vol * (por * RHO_FLUID * cf_cell
                       + (1.0 - por) * RHO_SOLID * cs_cell)
    I = math.fsum(np.asarray(integrand, dtype=np.float64))
    return I, data.GetNumberOfPoints()


def split_parts(fname):
    """Diagnostic: return (fluid_part, solid_part) separately, to show which
    side drifts if the test fails."""
    data = read_vtu(fname)
    pts = vtk_np(data.GetPoints().GetData()).astype(np.float64)
    c_fluid = vtk_np(get_array(data.GetPointData(), 'PassiveTracer_Cu')).astype(np.float64)
    c_solid = vtk_np(get_array(data.GetPointData(), 'Cu_solid')).astype(np.float64)
    tconn, tids = tet_cells(data)
    por_arr = get_array(data.GetCellData(), 'Porosity')
    por = (vtk_np(por_arr)[tids].astype(np.float64)
           if por_arr is not None else np.ones(len(tids)))
    a = pts[tconn[:, 0]]
    vol = np.abs(np.einsum('ij,ij->i',
                           pts[tconn[:, 1]] - a,
                           np.cross(pts[tconn[:, 2]] - a,
                                    pts[tconn[:, 3]] - a))) / 6.0
    fluid = math.fsum(np.asarray(vol * por * RHO_FLUID * c_fluid[tconn].mean(axis=1), dtype=np.float64))
    solid = math.fsum(np.asarray(vol * (1.0 - por) * RHO_SOLID * c_solid[tconn].mean(axis=1), dtype=np.float64))
    return fluid, solid


# Porosity contrast in the mpml: region A ~ 0.5, region B ~ 0.06. Classify each
# tet by its OWN porosity value (closer to 0.5 or closer to 0.06), rather than
# relying on a dumped RegionIds array, so this works regardless of what cell
# arrays IC-FERST happens to write out.
PHI_MID = 0.28  # midpoint between the 0.06 and 0.5 porosities in this mpml


def solid_mass_by_region(fname):
    """Return solid metal mass integrated separately over the high-porosity
    region and the low-porosity region."""
    data = read_vtu(fname)
    pts = vtk_np(data.GetPoints().GetData()).astype(np.float64)
    c_solid = vtk_np(get_array(data.GetPointData(), 'Cu_solid')).astype(np.float64)
    tconn, tids = tet_cells(data)
    por_arr = get_array(data.GetCellData(), 'Porosity')
    por = (vtk_np(por_arr)[tids].astype(np.float64)
           if por_arr is not None else np.ones(len(tids)))
    a = pts[tconn[:, 0]]
    vol = np.abs(np.einsum('ij,ij->i',
                           pts[tconn[:, 1]] - a,
                           np.cross(pts[tconn[:, 2]] - a,
                                    pts[tconn[:, 3]] - a))) / 6.0
    cs_cell = c_solid[tconn].mean(axis=1)
    mass_cell = vol * (1.0 - por) * RHO_SOLID * cs_cell
    high_phi = por >= PHI_MID
    m_high = math.fsum(np.asarray(mass_cell[high_phi], dtype=np.float64))
    m_low = math.fsum(np.asarray(mass_cell[~high_phi], dtype=np.float64))
    return m_high, m_low


# ---------------- compute the series ----------------
masses = []
nnodes = []
fluid_masses = []
solid_masses = []
region_high = []
region_low = []
if Passed:
    for i in range(AutoNumber + 1):
        m, n = metal_invariant(os.path.join(path, '%s_%d.vtu' % (AutoFile, i)))
        fl, so = split_parts(os.path.join(path, '%s_%d.vtu' % (AutoFile, i)))
        m_high, m_low = solid_mass_by_region(os.path.join(path, '%s_%d.vtu' % (AutoFile, i)))
        masses.append(m)
        nnodes.append(n)
        fluid_masses.append(fl)
        solid_masses.append(so)
        region_high.append(m_high)
        region_low.append(m_low)
        print('  dump %3d: nodes = %7d, I = %.10e  (fluid %.6e | solid %.6e)  '
              '[solid by region: high-phi %.6e | low-phi %.6e]'
              % (i, n, m, fl, so, m_high, m_low))

    I0 = masses[0]
    if abs(I0) < 1e-30:
        print('Initial invariant is zero: broken initial condition')
        Passed = False


def worst_and_total(series):
    ref = series[0]
    worst = 0.0
    for i in range(1, len(series)):
        worst = max(worst, abs(series[i] - series[i - 1]) / abs(ref))
    return worst, abs(series[-1] - ref) / abs(ref)


if Passed:
    # ---- THE pass/fail criterion: conservation of the total invariant I ----
    worst_I, total_I = worst_and_total(masses)
    changes = sum(1 for i in range(1, len(nnodes)) if nnodes[i] != nnodes[i - 1])

    print('worst per-event relative drift, total I : %.3e (tolerance %.1e)'
          % (worst_I, Tolerance_per_event))
    print('total relative drift, total I           : %.3e (tolerance %.1e)'
          % (total_I, Tolerance_total))
    print('dumps with changed node count            : %d (require >= %d)'
          % (changes, Min_mesh_changes))

    if worst_I > Tolerance_per_event or total_I > Tolerance_total:
        print('Total metal mass I is not conserved across adapt: the metal '
              'tracers are not being correctly weighted around the mesh-to-mesh '
              'projection (check the (1-phi)*rho_s weight for the solid tracer '
              'and the phi*rho*S weight for the fluid tracer).')
        Passed = False
    if changes < Min_mesh_changes:
        print('Mesh barely changed: adaptation did not exercise the projection')
        Passed = False

    # ---- DIAGNOSTICS ONLY (not pass/fail) ----
    worst_fluid, total_fluid = worst_and_total(fluid_masses)
    worst_solid, total_solid = worst_and_total(solid_masses)
    worst_high, total_high = worst_and_total(region_high)
    worst_low, total_low = worst_and_total(region_low)
    print('(diagnostic) M_fluid drift : worst %.3e, total %.3e (reaction-driven, not a gate)'
          % (worst_fluid, total_fluid))
    print('(diagnostic) M_solid drift : worst %.3e, total %.3e (reaction-driven, not a gate)'
          % (worst_solid, total_solid))
    print('(diagnostic) solid mass, high-phi region : worst %.3e, total %.3e'
          % (worst_high, total_high))
    print('(diagnostic) solid mass, low-phi region  : worst %.3e, total %.3e'
          % (worst_low, total_low))

if Passed:
    print('Metal mass conservation across adaptivity works OK')
else:
    print('Metal mass conservation across adaptivity does NOT work')


