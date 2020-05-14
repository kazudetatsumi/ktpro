#!/usr/bin/env python

import argparse
import numpy as np
import h5py
from phonopy.structure.symmetry import Symmetry
from phonopy.interface import read_crystal_structure
from phonopy.structure.cells import get_primitive
from phonopy.harmonic.force_constants import similarity_transformation
from phono3py.phonon3.triplets import (get_ir_grid_points,
                                       get_grid_points_by_rotations)
from phono3py.file_IO import write_kappa_to_hdf5

def fracval(frac):
    if frac.find('/') == -1:
        return float(frac)
    else:
        x = frac.split('/')
        return float(x[0]) / float(x[1])

def get_grid_symmetry(data):
    symmetry = data['symmetry']
    mesh = data['mesh']
    weights = data['weight']
    qpoints = data['qpoint']
    rotations = symmetry.get_pointgroup_operations()
    (ir_grid_points,
     weights_for_check,
     grid_address,
     grid_mapping_table) = get_ir_grid_points(mesh, rotations)

    np.testing.assert_array_equal(weights, weights_for_check)
    qpoints_for_check = grid_address[ir_grid_points] / mesh.astype('double')
    diff_q = qpoints - qpoints_for_check
    np.testing.assert_almost_equal(diff_q, np.rint(diff_q))

    return ir_grid_points, grid_address, grid_mapping_table

def expand(data):
    gv = data['group_velocity']
    qpoint = data['qpoint']
    frequency = data['frequency']
    cv = data['heat_capacity']
    if 'gamma_N' in data:
        g_N = data['gamma_N']
        g_U = data['gamma_U']
    symmetry = data['symmetry']
    primitive = data['cell']
    mesh = data['mesh']
    ir_grid_points = data['ir_grid_points']
    grid_address = data['grid_address']

    point_operations = symmetry.get_reciprocal_operations()
    rec_lat = np.linalg.inv(primitive.get_cell())
    rotations_cartesian = np.array(
        [similarity_transformation(rec_lat, r)
         for r in point_operations], dtype='double')

    gv_bz = np.zeros((len(grid_address),) + gv.shape[1:],
                     dtype='double', order='C')
    qpt_bz = np.zeros((len(grid_address), 3), dtype='double', order='C')
    freq_bz = np.zeros((len(grid_address), frequency.shape[1]),
                       dtype='double', order='C')
    cv_bz = np.zeros((cv.shape[0], len(grid_address), cv.shape[2]),
                     dtype='double', order='C')
    if 'gamma_N' in data:
        g_N_bz = np.zeros_like(cv_bz)
        g_U_bz = np.zeros_like(cv_bz)
    else:
        g_N_bz = None
        g_U_bz = None

    num_band = gv.shape[1]
    for i, gp in enumerate(ir_grid_points):
        rotation_map = get_grid_points_by_rotations(
            grid_address[gp],
            point_operations,
            mesh)
        multi = len(rotation_map) // len(np.unique(rotation_map))
        assert len(np.unique(rotation_map)) * multi == len(rotation_map)
        
        for rgp, r_c, r in zip(rotation_map,
                               rotations_cartesian,
                               point_operations):
            gv_bz[rgp] += np.dot(gv[i], r_c.T) / multi
            qpt_bz[rgp] = np.dot(r, qpoint[i])
            freq_bz[rgp] = frequency[i]
            cv_bz[:, rgp, :] = cv[:, i, :]
            if 'gamma_N' in data:
                g_N_bz[:, rgp, :] = g_N[:, i, :]
                g_U_bz[:, rgp, :] = g_U[:, i, :]

    return gv_bz, qpt_bz, freq_bz, cv_bz, g_N_bz, g_U_bz

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate cummulative physical properties")
    parser.add_argument("--pa", dest="primitive_matrix",
                        default="1 0 0 0 1 0 0 0 1", help="Primitive matrix")
    parser.add_argument("--yaml",
                        dest="write_yaml", action="store_true",
                        help="Write data to text in yaml format")
    parser.add_argument("--pwscf", action="store_true",
                        help="Pwscf mode")
    parser.add_argument('filenames', nargs='*')
    args = parser.parse_args()
    return args

def get_data(args, interface_mode=None):
    args = parse_args()
    cell, _ = read_crystal_structure(args.filenames[0],
                                     interface_mode=interface_mode)
    f = h5py.File(args.filenames[1])
    primitive_matrix = np.reshape(
        [fracval(x) for x in args.primitive_matrix.split()], (3, 3))
    primitive = get_primitive(cell, primitive_matrix)
    symmetry = Symmetry(primitive)

    data = {}
    data['cell'] = primitive
    data['symmetry'] = symmetry
    data['mesh'] = np.array(f['mesh'][:], dtype='intc') # (3)
    data['weight'] = f['weight'][:] # (gp)
    data['group_velocity'] = f['group_velocity'][:] # (gp, band, 3)
    data['qpoint'] = f['qpoint'][:] # (gp, 3)
    data['frequency'] = f['frequency'][:] # (gp, band)
    if 'gamma_N' in f:
        data['gamma_N'] = f['gamma_N'][:] # (temps, gp, band)
        data['gamma_U'] = f['gamma_U'][:] # (temps, gp, band)
    data['heat_capacity'] = f['heat_capacity'][:] # (temps, gp, band)
    data['temperature'] = np.array(f['temperature'][:], dtype='double') # (temps)

    ir_grid_points, grid_address, _ = get_grid_symmetry(data)
    data['ir_grid_points'] = ir_grid_points
    data['grid_address'] = grid_address

    return data

def write_yaml(data):
    gv_bz, qpt_bz, freq_bz, cv_bz, g_N_bz, g_U_bz = expand(data)
    print("temperature:")
    for t in data['temperature']:
        print("- %-10.1f" % t)
    print("phonon:")
    for i, (gv_q, f_q, q) in enumerate(zip(gv_bz, freq_bz, qpt_bz)):
        print("- q-position: [ %12.7f, %12.7f, %12.7f ] # %d" %
              (tuple(q) + (i + 1,)))
        print("  band:")
        for j, (gv_band, f_band) in enumerate(zip(gv_q, f_q)):
            print("  - frequency: %-15.10f" % f_band)
            print("    group_velocity: [ %13.7f, %13.7f, %13.7f ]" % tuple(gv_band))
            print("    heat_capacity:")
            for cv_temp in cv_bz[:, i, j]:
                print("    - %-15.10f" % cv_temp)
            print("    gamma_N_and_U:")
            for g_N_temp, g_U_temp in zip(g_N_bz[:, i, j], g_U_bz[:, i, j]):
                print("    - [ %12.7f, %12.7f ]" % (g_N_temp, g_U_temp))

def main():
    args = parse_args()
    if args.pwscf:
        data = get_data(args, interface_mode='pwscf')
    else:
        data = get_data(args)
    if args.write_yaml:
        write_yaml(data)
    else:
        gv_bz, qpt_bz, freq_bz, cv_bz, g_N_bz, g_U_bz = expand(data)
        filename = write_kappa_to_hdf5(data['temperature'],
                                       data['mesh'],
                                       frequency=freq_bz,
                                       group_velocity=gv_bz,
                                       heat_capacity=cv_bz,
                                       gamma_N=g_N_bz,
                                       gamma_U=g_U_bz,
                                       filename="bz",
                                       verbose=False)
        print("The data are written to %s." % filename)
                           
if __name__ == '__main__':
    main()
