#!/usr/bin/env python
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Tuple, Optional, Union
import numpy as np
import spglib
from scipy.sparse.linalg import LinearOperator, eigsh

# Physical constants (names kept close to your original code)
A0: float = 0.5291772   # Bohr [Å]
EH: float = 27.211      # 1 Hartree [eV]
MH_H: float = 1836.0    # proton mass [me]


# ----------------------------- Data container --------------------------------
@dataclass(frozen=True)
class Structure:
    """Lightweight structure container.

    Attributes
    ----------
    axis : ndarray, shape (3, 3)
        Lattice vectors in Å.
    positions : ndarray, shape (N, 3)
        Fractional coordinates.
    numbers : ndarray, shape (N,)
        Species indices or atomic numbers.
    """
    axis: np.ndarray
    positions: np.ndarray
    numbers: np.ndarray


# ------------------------------- POSCAR I/O ----------------------------------
class PoscarReader:
    """Read POSCAR/CONTCAR and return a normalized Structure (fractional)."""

    def read(self, path: str) -> Tuple[Structure, dict]:
        with open(path, 'r') as f:
            lines = f.readlines()
        #lines = open(path).read().splitlines()
        scale = float(lines[1].split()[0])
        lat = np.array(
            [list(map(float, lines[2 + i].split()[:3])) for i in range(3)]
        ) * scale

        nspc = np.asarray(lines[6].split(), dtype=int)
        numbers = np.array([i + 1 for i, c in enumerate(nspc)
                            for _ in range(c)])

        if "Di" in lines[7] or "di" in lines[7]:
            iniline = 8
        elif "Sel" in lines[7] or "sel" in lines[7]:
            iniline = 9

        pos = np.array(
            [list(map(float, lines[iniline + i].split()[:3]))
                for i in range(numbers.size)]
        )

        s = Structure(axis=lat, positions=pos, numbers=numbers)
        meta = {"lines": lines, "iniline": iniline, }
        return s, meta


# -------------------------------- Symmetry -----------------------------------
class SymmetryAnalyzer:
    """Thin spglib wrapper (primitive/refined cell and operations)."""

    def __init__(self, prim: bool = True, symprec: float = 1e-5):
        self.prim = prim
        self.symprec = symprec

    def refine(self, s: Structure) -> Tuple[Structure, dict]:
        """Return primitive/refined cell with numbers-sorted positions."""
        cell = (s.axis, s.positions, s.numbers)
        if self.prim:
            lat, pos, num = spglib.find_primitive(cell)
        else:
            lat, pos, num = spglib.refine_cell(cell)
        idx = np.argsort(num)
        s2 = Structure(axis=lat, positions=pos[idx], numbers=np.sort(num))
        return s2, {"sorted_index": idx}

    def operations(self, s: Structure) -> Tuple[np.ndarray, np.ndarray]:
        """Return symmetry rotations and translations."""
        cell = (s.axis, s.positions, s.numbers)
        sym = spglib.get_symmetry(cell, symprec=self.symprec)
        return sym["rotations"], sym["translations"]


def strip_target_H(s: Structure, center_frac: Iterable[float],
                   tol: float = 1e-5) -> Structure:
    """Remove the H at center_frac (fractional basis) from Structure."""
    c = np.asarray(center_frac, float)
    keep = np.sum(np.abs((s.positions - c) % 1.0), axis=1) >= tol
    return Structure(axis=s.axis, positions=s.positions[keep],
                     numbers=s.numbers[keep])


# --------------------------------- H-grid ------------------------------------
def make_hgrid_cube(center_frac: Iterable[float],
                    edgelength: float | Iterable[float],
                    nx: int, *, old_style: bool = False) -> np.ndarray:
    """Return a cubic grid around H in fractional coordinates.

    Parameters
    ----------
    center_frac : iterable of float
        Fractional center (e.g., (0.5, 0.5, 0.5)).
    edgelength : float or iterable of float
        Edge length in fractional units (scalar or per-axis).
    nx : int
        Grid count per axis.
    old_style : bool, default False
        If True, spacing = L/(nx - 1); else spacing = L/nx.

    Returns
    -------
    ndarray, shape (nx^3, 3)
        Fractional grid points.
    """
    c = np.asarray(center_frac, float)
    L = (np.asarray(edgelength, float)
         if np.ndim(edgelength) else np.array([edgelength] * 3, float))
    dx = L / (nx - 1) if old_style else L / nx
    pts = []
    for ix in range(nx):
        for iy in range(nx):
            for iz in range(nx):
                p = c - L / 2 + np.array([ix, iy, iz]) * dx
                pts.append(p % 1.0)
    return np.asarray(pts)


# -------------- Irreducible reps (from original GetIrreducibleShift) ---------
def get_irreducible_points(hpos: np.ndarray,
                           rotations: np.ndarray,
                           translations: np.ndarray,
                           tol: float = 1e-5) -> np.ndarray:
    """Remove symmetry-equivalent points and return irreducible reps.

    Notes
    -----
    Logic matches the original class method:
      - apply each non-identity operation
      - delete points that become identical to others by that operation
      - keep self-point at each step
    """
    irr_hpos = np.array(hpos)
    print('chkchk, irr_hpos', irr_hpos.shape)
    isym = 0
    for rot, trans in zip(rotations[1:], translations[1:]):
        isym += 1
        sym_hpos = ((np.matmul(rot, irr_hpos.T)).T + trans) % 1.0
        sym_hpos[np.abs(sym_hpos - 1.0) <= 1e-14] = 0.
        cond = np.ones((irr_hpos.shape[0]), dtype=bool)
        for i, _irr_hpos in enumerate(irr_hpos):
            if cond[i]:
                cond *= np.sum(np.abs(_irr_hpos - sym_hpos) % 1.0,
                               axis=1) >= 1e-5
                cond[i] = True
        irr_hpos = irr_hpos[cond]
    print('chk, irr_hpos', irr_hpos.shape)
    # print(irr_hpos)
    return irr_hpos


# ------- Mapping all-points -> rep index (from original GetDataOverAllHpos4) -
def get_mapping(hpos: np.ndarray,
                irr_hpos: np.ndarray,
                rotations: np.ndarray,
                translations: np.ndarray,
                tol: float = 1e-5) -> np.ndarray:
    """Return mapping (irr_idx): for each all-point, the rep index it belongs to.

    Notes
    -----
    Vectorized, equivalent to your GetDataOverAllHpos4:
      - build all symmetry images of rep points
      - compare with all points under periodicity
      - if any op makes them equal (< tol), assign that rep
    """
    # rot(na,3,3) @ irr(3,Nir) -> (na,3,Nir) -> (Nir,na,3)
    test = (rotations @ irr_hpos.T).transpose(2, 0, 1)
    test = (test + translations[None, :, :]) % 1.0
    test[np.abs(test - 1.0) <= 1e-14] = 0.0

    # Broadcast to (Nir, na, Nall, 3)
    diff = (hpos[None, None, :, :] - test[:, :, None, :]) % 1.0
    cond = np.sum(np.abs(diff), axis=3) >= tol       # (Nir, na, Nall)
    cond = np.prod(cond, axis=1, dtype=bool)         # (Nir, Nall)

    pairs = np.where(~cond)                          # indices where equal
    irr_idx = pairs[0][np.argsort(pairs[1])]         # (Nall,)
    return irr_idx

    ##original GetDataOverAllHpos4 moidifed as a function.
    #test = (np.matmul(rotations,
    #                  irr_hpos.T).transpose(2, 0, 1)
    #        + translations) % 1.0
    #test[np.abs(test - 1.) <= 1e-15] = 0.
    #test2 = np.repeat(np.expand_dims(test, 2), hpos.shape[0], axis=2)
    #cond = np.prod(np.sum(np.abs(hpos - test2) % 1.0,
    #                      axis=3) >= 1e-5, axis=1, dtype=bool)
    #pairs = np.where(np.invert(cond))
    #irr_idx = pairs[0][np.argsort(pairs[1])]
    #print('chk, irr_idx', irr_idx.shape, type(irr_idx))
    #return irr_idx


# ----- Write POSCAR files with target H at each irreducible position -----
def generate_shifted_poscar(structure: Structure,
                            irr_hpos: np.ndarray,
                            meta: dict,
                            center_frac) -> None:
    """
    Write POSCAR_### by replacing the target H (at center_frac, fractional)
    with each irreducible position in `irr_hpos`.
    Only this function is added/modified to minimize diffs.
    """
    lines = list(meta['lines'])
    head = int(meta['iniline'])
    # Update lattice vectors
    lat = np.asarray(structure.axis, float)
    for il in range(3):
        lines[il+2] = "     {:.16f}    {:.16f}    {:.16f} \n"\
                      .format(lat[il, 0], lat[il, 1], lat[il, 2])
        # Overwrite positions and flags (fractional)
    pos = np.asarray(structure.positions, float)
    for il in range(pos.shape[0]):
        lines[head+il] = "  {:.16f}  {:.16f}  {:.16f} \n"\
                         .format(pos[il, 0], pos[il, 1], pos[il, 2])
        # locate target H index
    c = np.asarray(center_frac, float)
    diffs = np.sum(np.abs((pos - c) % 1.0), axis=1)
    hits = np.where(diffs < 1e-5)[0]
    h_idx = int(hits[0]) if hits.size > 0 else int(np.argmin(diffs))
    # filename width
    nd = int(irr_hpos.shape[0])
    width = 5 if nd > 9999 else (4 if nd > 999 else 3)
    for iridx, ir in enumerate(irr_hpos):
        out_lines = lines.copy()
        out_lines[head+h_idx] = "  {:.16f}  {:.16f}  {:.16f} \n"\
                                .format(ir[0], ir[1], ir[2])
        outfile = 'POSCAR_{}'.format(str(iridx+1).zfill(width))
        with open(outfile, 'w') as f:
            f.writelines(out_lines)
    return None


# ---------------------------- Potential reconstruction -----------------------
def build_potential_3d(irr_E: np.ndarray, mapping: np.ndarray, nx: int
                       ) -> np.ndarray:
    """Fill a 3D cube with irreducible energies via mapping."""
    return irr_E[mapping].reshape(nx, nx, nx)


# ----------------------------------- G-space ---------------------------------
def gen_gset(nx: int, rg_idx: int | None = None) -> np.ndarray:
    """Generate index-space sphere-cut G set."""
    if rg_idx is None:
        rg_idx = (nx - 1) // 2 - 1
    Gs = []
    for i in range(-nx // 2, nx // 2 + 1):
        for j in range(-nx // 2, nx // 2 + 1):
            for k in range(-nx // 2, nx // 2 + 1):
                if i * i + j * j + k * k < rg_idx * rg_idx:
                    Gs.append([i, j, k])
    return np.asarray(Gs, dtype=int)


def kinetic_diag(Gs: np.ndarray, lattice: np.ndarray, mh: float,
                 *, a0: float = A0, kvec: Iterable[float] | None = None
                 ) -> np.ndarray:
    """Return diagonal kinetic energy in Hartree."""
    A = lattice / a0
    BinvT = np.linalg.inv(A).T
    Gcart = (BinvT @ Gs.T).T
    if kvec is not None:
        kcart = BinvT @ np.asarray(kvec, float)
        Gcart = Gcart + kcart
    return (2 * np.pi) ** 2 * np.sum(Gcart ** 2, axis=1) / (2.0 * mh)


# ------------------------------ V operator (FFT) -----------------------------
def fft_coeff(V3d: np.ndarray, *, Eh: float = EH,
              subtract_mean: bool = True) -> np.ndarray:
    """Return FFT coefficients of V/Eh (norm='forward')."""
    arr = (V3d - np.mean(V3d)) / Eh if subtract_mean else (V3d / Eh)
    return np.fft.fftn(arr, norm='forward')


def make_matvec(Gs: np.ndarray, z: np.ndarray, Kdiag: np.ndarray, nx: int):
    """Return Hψ matvec in the G basis."""
    gmod = (Gs % nx).astype(int)
    gidx = (gmod[:, 0], gmod[:, 1], gmod[:, 2])

    def mv(x: np.ndarray) -> np.ndarray:
        X = np.zeros((nx, nx, nx), dtype=np.complex128)
        X[gidx] = x
        phi_r = np.fft.ifftn(X)
        Y_full = np.fft.fftn(phi_r) * z
        y = Y_full[gidx] + Kdiag * x
        return y.astype(np.complex128)

    return mv


# ---------------------------------- Solver -----------------------------------
def solve_krylov(matvec, Ng: int, *, k: int = 30, which: str = "SA",
                 tol: float = 1e-8, Eh: float = EH,
                 maxiter: int | None = None) -> Tuple[np.ndarray, np.ndarray]:
    """Compute lowest eigenpairs via scipy.sparse.linalg.eigsh."""
    A = LinearOperator((Ng, Ng), matvec=matvec, dtype=np.complex128)
    E_h, U = eigsh(A, k=k, which=which, tol=tol, maxiter=maxiter)
    return E_h * Eh * 1000.0, U  # meV


# ----------------------------------- main ------------------------------------
def main(poscar_path="CONTCAR", energies_path="ENERGIES",
         center_frac=(0.5, 0.5, 0.5), edgelength=0.4, nx=14,
         prim=True, symprec=1e-5, old_style=False,
         kvec=None, subtract_mean=True, mh: float = MH_H):
    """Minimal pipeline: POSCAR -> symmetry -> H-grid -> irreps -> eigenpairs.

    Steps
    -----
    1) Read POSCAR into Structure.
    2) Refine structure and get symmetry operations on the H-removed cell.
    3) Make a cubic fractional H-grid of size nx^3 around center_frac.
    4) Compute irreducible representatives and the all->rep mapping (irr_idx).
    5) Load irreducible energies and reconstruct a 3D potential.
    6) Build G set, kinetic diagonal, FFT coeffs, and Hψ matvec.
    7) Solve k lowest eigenpairs by Krylov method.

    Returns
    -------
    E_meV : ndarray
        Eigenvalues in meV.
    U : ndarray
        Eigenvectors (G-basis).
    extras : dict
        Contains irr_pts, irr_idx, V3d, Gs for further analysis.
    """
    # (1) POSCAR -> Structure
    s, _ = PoscarReader().read(poscar_path)

    # (2) Symmetry (remove one H at center_frac before getting operations)
    sym = SymmetryAnalyzer(prim=prim, symprec=symprec)
    s_ref, _ = sym.refine(s)
    s_wo_H = strip_target_H(s_ref, np.asarray(center_frac, float))
    R, T = sym.operations(s_wo_H)

    # (3) H-grid (fractional)
    hpoints = make_hgrid_cube(center_frac=center_frac,
                              edgelength=edgelength,
                              nx=nx, old_style=old_style)

    # (4) Irreps and mapping (irr_idx)
    irr_pts = get_irreducible_points(hpoints, R, T, tol=1e-5)
    irr_idx = get_mapping(hpoints, irr_pts, R, T, tol=1e-5)

    # (5) Potential reconstruction from irreducible energies
    irr_E = np.loadtxt(energies_path)  # one-column assumed
    V3d = build_potential_3d(irr_E, irr_idx, nx)

    # (6) G set, kinetic, and V operator
    Gs = gen_gset(nx=nx)
    Kdiag = kinetic_diag(Gs, lattice=s_ref.axis, mh=mh, kvec=kvec)
    z = fft_coeff(V3d, Eh=EH, subtract_mean=subtract_mean)
    matvec = make_matvec(Gs=Gs, z=z, Kdiag=Kdiag, nx=nx)

    # (7) Krylov eigen-solve
    E_meV, U = solve_krylov(matvec, Ng=Gs.shape[0], k=30,
                            which="SA", tol=1e-8)
    print("min(E) [meV]:", float(np.min(E_meV)))
    extras = dict(irr_pts=irr_pts, irr_idx=irr_idx, V3d=V3d, Gs=Gs)
    return E_meV, U, extras


if __name__ == "__main__":
    main()


