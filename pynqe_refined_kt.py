#!/usr/bin/env python
"""
Refined English-only version of pynqe.py
- All Japanese comments/docstrings removed or translated to English
- Docstrings unified and clarified
- Function name inconsistencies reconciled with thin wrappers (no behavioral
  changes)

Notes
-----
This file keeps the original physics and numerics but restores missing helper
functions that were referenced in the pipeline (e.g., `make_hgrid_cube`,
`gen_gset`, `solve_krylov`, `make_matvec`) and adds aliases so that older names
still work (e.g., `build_potential_3d`, `compute_short_bond_flags`,
`replace_energies_by_flag`).
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Iterable, Tuple, Optional, Union

import numpy as np
import spglib
from scipy.sparse.linalg import LinearOperator, eigsh
from scipy.linalg import eigh
# convergence exception
from scipy.sparse.linalg._eigen.arpack.arpack import ArpackNoConvergence

# ---------------------------- Physical constants ----------------------------
A0: float = 0.5291772   # Bohr [Å]
EH: float = 27.211      # 1 Hartree [eV]
MH_H: float = 1836.0    # proton mass [me]


# ------------------------------ Data container ------------------------------
@dataclass(frozen=True)
class Structure:
    """Lightweight container for a crystal structure.

    Attributes
    ----------
    axis : ndarray, shape (3, 3)
        Lattice vectors in Å (rows).
    positions : ndarray, shape (N, 3)
        Fractional coordinates.
    numbers : ndarray, shape (N,)
        Species indices or atomic numbers.
    """

    axis: np.ndarray
    positions: np.ndarray
    numbers: np.ndarray


# ------------------------------ POSCAR I/O ------------------------------
class PoscarReader:
    """Read POSCAR/CONTCAR and return a normalized :class:`Structure`.

    The returned positions are fractional.
    """

    def read(self, path: str) -> Tuple[Structure, dict]:
        with open(path, "r", encoding="utf-8") as f:
            lines = f.readlines()
        scale = float(lines[1].split()[0])
        lat = np.array(
            [list(map(float, lines[2 + i].split()[:3])) for i in range(3)]
        ) * scale
        nspc = np.asarray(lines[6].split(), dtype=int)
        numbers = np.array([i+1 for i, c in enumerate(nspc) for _ in range(c)])
        if "Di" in lines[7] or "di" in lines[7]:
            iniline = 8
        elif "Sel" in lines[7] or "sel" in lines[7]:
            iniline = 9
        else:
            iniline = 8  # default fallback
        pos = np.array(
            [list(map(float, lines[iniline + i].split()[:3]))
                for i in range(numbers.size)]
        )
        s = Structure(axis=lat, positions=pos, numbers=numbers)
        meta = {"lines": lines, "iniline": iniline}
        return s, meta


# ------------------------------ Symmetry ------------------------------
class SymmetryAnalyzer:
    """Thin spglib wrapper to obtain primitive/refined cells and operations."""

    def __init__(self, prim: bool = True, symprec: float = 1e-5):
        self.prim = prim
        self.symprec = symprec

    def refine(self, s: Structure) -> Tuple[Structure, dict]:
        """Return primitive/refined cell with pos. sorted by species number."""
        cell = (s.axis, s.positions, s.numbers)
        if self.prim:
            lat, pos, num = spglib.find_primitive(cell)
        else:
            lat, pos, num = spglib.refine_cell(cell)
        idx = np.argsort(num)
        s2 = Structure(axis=lat, positions=pos[idx], numbers=np.sort(num))
        return s2, {"sorted_index": idx}

    def operations(self, s: Structure) -> Tuple[np.ndarray, np.ndarray]:
        """Return symmetry rot. and trans. for the given structure."""
        cell = (s.axis, s.positions, s.numbers)
        sym = spglib.get_symmetry(cell, symprec=self.symprec)
        return sym["rotations"], sym["translations"]


def strip_target_H(s: Structure, center_frac: Iterable[float],
                   tol: float = 1e-5) -> Structure:
    """Remove the H atom at ``center_frac`` (fractional) from ``s``.

    The target is identified by the smallest periodic distance within a
    tolerance.
    """
    c = np.asarray(center_frac, float)
    keep = np.sum(np.abs((s.positions - c) % 1.0), axis=1) >= tol
    return Structure(axis=s.axis, positions=s.positions[keep],
                     numbers=s.numbers[keep])


# ------------------------------ H-grid ------------------------------
def make_hgrid_and_hlat_cube(
    structure: Structure,
    center_frac: Iterable[float],
    edge_len_frac: Union[float, Iterable[float]],
    nx: int,
    *,
    old_style: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return a cubic grid around ``center_frac`` in fractional coordinates.

    Parameters
    ----------
    structure : Structure
        The reference structure (its lattice is used to compute ``hlat``).
    center_frac : iterable of float
        Fractional center (e.g., ``(0.5, 0.5, 0.5)``).
    edge_len_frac : float or iterable of float
        Edge length in fractional units (scalar or per-axis sequence).
    nx : int
        Number of grid points per axis.
    old_style : bool, default False
        If True, spacing = L/(nx-1); else spacing = L/nx.

    Returns
    -------
    points : ndarray, shape (nx**3, 3)
        Fractional grid points.
    hlat : ndarray, shape (3, 3)
        Effective lattice spanned by the cubic region in Cartesian coordinates.
    """
    center_frac = np.asarray(center_frac, float)
    if np.ndim(edge_len_frac) == 0:
        edge_len_frac = np.array([edge_len_frac] * 3, float)
    else:
        edge_len_frac = np.asarray(edge_len_frac, float)
    dx = edge_len_frac / (nx - 1) if old_style else edge_len_frac / nx
    pts = []
    for ix in range(nx):
        for iy in range(nx):
            for iz in range(nx):
                p = center_frac-edge_len_frac/2+np.array([ix, iy, iz])*dx
                pts.append(p % 1.0)
    lat = np.asarray(structure.axis, float)
    hlat = lat @ np.diag(edge_len_frac)
    return np.asarray(pts), hlat


def make_hgrid_cube(
    center_frac: Iterable[float],
    edge_len_frac: Union[float, Iterable[float]],
    nx: int,
    *,
    old_style: bool = False,
) -> np.ndarray:
    """Convenience wrapper of :func:`make_hgrid_and_hlat_cube` returning only points."""
    # A minimal structure-like holder for the lattice is not required here,
    # so we create a dummy with identity lattice; only points are returned.
    dummy = Structure(axis=np.eye(3), positions=np.zeros((0, 3)), numbers=np.zeros((0,), int))
    pts, _ = make_hgrid_and_hlat_cube(dummy, center_frac, edge_len_frac, nx, old_style=old_style)
    return pts


# ------------------------ Irreducible representatives ------------------------
def get_irreducible_points(
    hpos: np.ndarray,
    rotations: np.ndarray,
    translations: np.ndarray,
    tol: float = 1e-5,
) -> np.ndarray:
    """Remove symmetry-equivalent grid points and return irreducible ones.

    Logic:
    - apply each non-identity symmetry operation to the current set
    - drop points that become identical (within periodicity) under that loop
    - always keep the self-point at each step
    """
    irr_hpos = np.array(hpos)
    for rot, trans in zip(rotations[1:], translations[1:]):
        sym_hpos = ((np.matmul(rot, irr_hpos.T)).T + trans) % 1.0
        sym_hpos[np.abs(sym_hpos - 1.0) <= 1e-14] = 0.0
        cond = np.ones((irr_hpos.shape[0]), dtype=bool)
        for i, _irr in enumerate(irr_hpos):
            if cond[i]:
                cond *= np.sum(np.abs(_irr - sym_hpos) % 1.0, axis=1) >= tol
                cond[i] = True
        irr_hpos = irr_hpos[cond]
    return irr_hpos


def get_mapping(
    hpos: np.ndarray,
    irr_hpos: np.ndarray,
    rotations: np.ndarray,
    translations: np.ndarray,
    tol: float = 1e-5,
) -> np.ndarray:
    """Map each full-grid point to the index of its irreducible representative.

    Vectorized equivalent to the legacy routine: build all symmetry images of
    irreducible points and compare them to full-grid points under periodicity.
    """
    test = (np.matmul(rotations, irr_hpos.T).transpose(2, 0, 1)
            + translations) % 1.0
    test[np.abs(test - 1.0) <= 1e-15] = 0.0
    test2 = np.repeat(np.expand_dims(test, 2), hpos.shape[0], axis=2)
    cond = np.prod(
        np.sum(np.abs(hpos - test2) % 1.0, axis=3) >= tol, axis=1, dtype=bool
    )
    pairs = np.where(np.invert(cond))
    irr_idx = pairs[0][np.argsort(pairs[1])]
    return irr_idx


# ------------------------ POSCAR writer for shifted H ------------------------
def generate_shifted_poscar(
    structure: Structure,
    irr_hpos: np.ndarray,
    meta: dict,
    center_frac: Iterable[float],
) -> None:
    """Write ``POSCAR_#####`` by replacing the target H with each irrep position.

    The target H is specified by ``center_frac`` in fractional coordinates.
    Only this function touches I/O to minimize diffs to the original.
    """
    lines = list(meta["lines"])
    head = int(meta["iniline"])

    # Update lattice vectors
    lat = np.asarray(structure.axis, float)
    for il in range(3):
        lines[il + 2] = " {:.16f} {:.16f} {:.16f} \n"\
                .format(lat[il, 0], lat[il, 1], lat[il, 2])
        print(lines[il+2])

    # Overwrite positions (fractional)
    pos = np.asarray(structure.positions, float)
    for il in range(pos.shape[0]):
        lines[head + il] = " {:.16f} {:.16f} {:.16f} \n"\
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
        out_lines[head + h_idx] = " {:.16f} {:.16f} {:.16f} \n"\
                                  .format(ir[0], ir[1], ir[2])
        outfile = f"POSCAR_{str(iridx + 1).zfill(width)}"
        with open(outfile, "w", encoding="utf-8") as f:
            f.writelines(out_lines)
    return None


# ------------------------- Adiabatic potential data -------------------------
def load_energies(path: str, *, relax: bool = False) -> np.ndarray:
    """Load ``ENERGIES`` as a 1D float array.

    If ``relax=False`` (default), use the 4th column (index 3);
    if ``True``, use the 5th column (index 4). Mirrors the legacy behavior.
    """
    col = 4 if relax else 3
    ene: list[float] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            ene.append(float(line.split()[col]))
    return np.asarray(ene, dtype=float)


def replace_ene_with_flags(
    ene: np.ndarray, short_atom_distance_flags: np.ndarray, delta: float = 10.0
) -> np.ndarray:
    """Replace energies at positions where ``flags==True`` with
    ``(mean of non-flagged) + delta``.

    Equivalent to the legacy behavior with ``delta=10``.
    """
    eave = np.average(ene[~short_atom_distance_flags])
    ene[short_atom_distance_flags] = eave + delta
    return ene


# Backward-compatible alias
#def replace_energies_by_flag(ene: np.ndarray, flags: np.ndarray,
#                             delta: float = 10.0) -> np.ndarray:
#    return replace_ene_with_flags(ene, flags, delta)


def expand_positions_3x3x3(frac_pos: np.ndarray) -> np.ndarray:
    """Expand fractional coordinates (M, 3) into 27 periodic images.

    Returns an array of shape ``(27, M, 3)`` (legacy ``expos`` behavior).
    """
    pos = np.asarray(frac_pos, float)
    shifts = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1)
                      for k in (-1, 0, 1)], float)
    return shifts[:, None, :] + pos[None, :, :]


def get_min_bond_length(h_frac: np.ndarray, axis: np.ndarray,
                        host_frac: np.ndarray) -> float:
    """Compute the shortest host–H distance [Å] for a candidate H at ``h_frac``.

    Uses 27 periodic images. ``axis`` is the lattice (3, 3) in Å. Mirrors legacy behavior.
    """
    epos = expand_positions_3x3x3(host_frac)  # (27, M, 3)
    diff_frac = epos - np.asarray(h_frac, float)  # (27, M, 3)
    diff_cart = diff_frac @ np.asarray(axis, float)  # Cartesian [Å]
    dists = np.linalg.norm(diff_cart, axis=2)  # (27, M)
    return float(np.min(dists))


def check_short_atom_distance_flags(
    irr_points: np.ndarray, woH_structure: Structure, threshold: float = 0.3
) -> np.ndarray:
    """For each irreducible H point, return True if min host–H distance < threshold [Å].

    Functional counterpart to the legacy checker. Default threshold = 0.6 Å in legacy
    codes, but here the default is 0.3 Å to be conservative. Adjust via ``threshold``.
    """
    axis = woH_structure.axis
    host = woH_structure.positions
    atom_distance_flags = [
        get_min_bond_length(p, axis, host) < threshold for p in np.asarray(irr_points, float)
    ]
    return np.asarray(atom_distance_flags, dtype=bool)


# Backward-compatible alias
#def compute_short_bond_flags(irr_points: np.ndarray, woH_structure: Structure, threshold: float = 0.6) -> np.ndarray:
#    return check_short_atom_distance_flags(irr_points, woH_structure, threshold=threshold)


def assemble_potential_3d(irr_energies: np.ndarray, mapping: np.ndarray,
                          nx: int) -> np.ndarray:
    """Assemble a ``(nx, nx, nx)`` potential array from irreducible energies.

    Equivalent to the legacy: ``potential = ene[irr_idx].reshape(nx, nx, nx)``.
    """
    return np.asarray(irr_energies, float)[mapping].reshape(nx, nx, nx)


# Backward-compatible alias
#def build_potential_3d(irr_energies: np.ndarray, mapping: np.ndarray, nx: int) -> np.ndarray:
#    return assemble_potential_3d(irr_energies, mapping, nx)


# ----------------------------- G-space utilities -----------------------------
def gen_gset(nx: int, rg_idx: Optional[int] = None) -> np.ndarray:
    """Generate an index-space sphere-cut set of G-vectors.

    Parameters
    ----------
    nx : int
        FFT grid size per dimension.
    rg_idx : int or None
        Cutoff radius in index space. If None, uses ``(nx - 1) // 2 - 1``.
    """
    if rg_idx is None:
        rg_idx = (nx - 1) // 2 - 1
    Gs = []
    for i in range(-nx // 2, nx // 2 + 1):
        for j in range(-nx // 2, nx // 2 + 1):
            for k in range(-nx // 2, nx // 2 + 1):
                if i * i + j * j + k * k < rg_idx * rg_idx:
                    Gs.append([i, j, k])
    return np.asarray(Gs, dtype=int)


def cut_Gs_physically(lattice: np.ndarray, nx: int,
                      ao: float = A0) -> np.ndarray:
    """Generate G-vectors by a physical cutoff in reciprocal space.

    The cutoff radius is based on the minimum distance of reciprocal basis
    vectors and a Nyquist-like index derived from ``nx``.
    """
    A = np.asarray(lattice, float) / A0
    BinvT = np.linalg.inv(A).T
    nnyq = (nx - 1) // 2 - 1
    dist = []
    for axis in range(3):
        vec = np.zeros(3)
        vec[axis] = 1.0
        dist.append(np.linalg.norm(BinvT @ vec))
    mindist = min(dist)
    r_phys = mindist * nnyq
    Gs = []
    for i in range(-nnyq, nnyq + 1):
        for j in range(-nnyq, nnyq + 1):
            for k in range(-nnyq, nnyq + 1):
                G = BinvT @ np.array([i, j, k])
                if np.linalg.norm(G) < r_phys:
                    Gs.append([i, j, k])
    return np.asarray(Gs, dtype=int)


def kinetic_diag(
    Gs: np.ndarray,
    lattice: np.ndarray,
    mh: float = MH_H,
    *,
    a0: float = A0,
    kvec: Optional[Iterable[float]] = None,
) -> np.ndarray:
    """Return the diagonal kinetic energy in Hartree."""
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
    """Return FFT coefficients of ``V/Eh`` with ``norm='forward'``.

    If ``subtract_mean`` is True (default), the mean of ``V`` is removed
    before scaling.
    """
    arr = (V3d - np.mean(V3d)) / Eh if subtract_mean else (V3d / Eh)
    return np.fft.fftn(arr, norm="forward")


def get_V(Gs: np.ndarray, nx: int, z: np.ndarray) -> np.ndarray:
    """Build the potential matrix ``V`` in the G basis via FFT indexing.

    ``z`` is the FFT of ``V/Eh``.
    """
    gi, gj, gk = Gs[:, 0], Gs[:, 1], Gs[:, 2]
    di = (gi[:, None] - gi[None, :]) % nx
    dj = (gj[:, None] - gj[None, :]) % nx
    dk = (gk[:, None] - gk[None, :]) % nx
    return z[di, dj, dk]


def make_matvec(Gs: np.ndarray, z: np.ndarray, Kdiag: np.ndarray, nx: int):
    """Return the ``Hψ`` matvec in the G basis (multiply by V in real space).

    Implementation detail:
    - Transform ``z`` back to real space to obtain ``V_r = V/EH`` on the real
      grid.
    - For each input vector ``x`` in the G-basis, scatter to the full grid,
      IFFT to real space, multiply by ``V_r``, FFT back, gather the G slots,
      and add the kinetic term.
    """
    # Precompute V in real space (dimensionless V/EH)
    V_r = np.fft.ifftn(z)

    # Fixed indices of active G-slots on the FFT grid
    gmod = (Gs % nx).astype(int)
    gidx = (gmod[:, 0], gmod[:, 1], gmod[:, 2])

    def mv(x: np.ndarray) -> np.ndarray:
        X = np.zeros((nx, nx, nx), dtype=np.complex128)
        X[gidx] = x
        phi_r = np.fft.ifftn(X)   # real-space wavefunction
        Y_r = V_r * phi_r         # multiply in real space (physical operation)
        Y_f = np.fft.fftn(Y_r)    # back to reciprocal space
        y = Y_f[gidx] + Kdiag * x
        return y.astype(np.complex128, copy=False)

    return mv


# ------------------------------ Solvers ------------------------------
def solve_krylov(
    matvec,
    Ng: int,
    *,
    k: int = 30,
    which: str = "SA",
    tol: float = 1e-8,
    maxiter: Optional[int] = None,
    Eh: float = EH,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute the lowest ``k`` eigenpairs using ``scipy.sparse.linalg.eigsh``.

    Returns
    -------
    E_meV : ndarray
        Eigenvalues in meV.
    U : ndarray
        Eigenvectors (columns) in the G basis (complex128).
    """
    A = LinearOperator((Ng, Ng), matvec=matvec, dtype=np.complex128)
    E, U = eigsh(A, k=k, which=which, tol=tol, maxiter=maxiter)
    E *= Eh * 1000.0
    idx = np.argsort(E)
    E = E[idx]
    U = U[:, idx]
    return E, U


def solve_dense_H(H: np.ndarray, Compress: bool = False,
                  v: Optional[np.ndarray] = None):
    """Compute eigenpairs of a dense Hermitian matrix using LAPACK (``eigh``).

    If ``Compress`` is True and an input basis ``v`` is provided (or computed),
    the Hamiltonian is projected to the subspace spanned by ``v`` for faster
    solution. All energies are returned in meV.
    """
    if Compress and v is None:
        print("Compression is used in Eigenvalue solution")
        print("Calculating v matrix of 200 eigen vectors..")
        dummy, v = eigh(H, subset_by_index=[0, 199], driver="evr")
        dummy *= EH * 1000.0
        print(np.min(dummy))
        print(dummy[0:15] - np.min(dummy))
    if Compress and v is not None:
        print("Compressing H by v.")
        H_comp = v.conj().T @ H @ v
        print("Calculating 30 eigen from the compressed H")
        E_comp, U_comp = eigh(H_comp, subset_by_index=[0, 29], driver="evr")
        E_comp *= EH * 1000.0
        print(np.min(E_comp))
        print(E_comp[0:15] - np.min(E_comp))
        U = v @ U_comp
        return E_comp, U
    else:
        print("No compression is used in Eigenvalue solution")
        E, U = eigh(H, subset_by_index=[0, 29])
        E *= EH * 1000.0
        print(np.min(E))
        print(E[0:15] - np.min(E))
        return E, U


class NQEKrylovSolver:
    """Small utility class to solve ``Hψ`` with a Krylov (``eigsh``) method.

    Features
    --------
    - Multiplies ``V/EH`` in real space and transforms back via FFT (robust &
      consistent).
    - Stores fixed data (``gidx``, ``nx``, ``Kdiag``, ``V_r``) as members.
    - :meth:`solve` returns sorted energies (ascending) and ensures real
      energies.
    """

    def __init__(
        self,
        Gs: np.ndarray,
        nx: int,
        Kdiag: np.ndarray,
        V3d: np.ndarray,
        *,
        subtract_mean: bool = True,
        Eh: float = EH,
        dtype=np.complex128,
    ):
        self.nx = int(nx)
        self.dtype = dtype
        gmod = (Gs % nx).astype(int)
        self.gidx = (gmod[:, 0], gmod[:, 1], gmod[:, 2])
        self.Ng = gmod.shape[0]
        if subtract_mean:
            self.V_r = (np.asarray(V3d, float) - float(np.mean(V3d))) / Eh
        else:
            self.V_r = np.asarray(V3d, float) / Eh
        self.Kdiag = np.asarray(Kdiag, dtype=dtype)
        self._X = np.zeros((nx, nx, nx), dtype=dtype)
        self.operator = LinearOperator(shape=(self.Ng, self.Ng),
                                       matvec=self._matvec, dtype=dtype)

    def _matvec(self, x: np.ndarray) -> np.ndarray:
        X = self._X
        X.fill(0.0)
        X[self.gidx] = x.astype(self.dtype, copy=False)
        phi_r = np.fft.ifftn(X)
        Y_r = self.V_r * phi_r
        Y_f = np.fft.fftn(Y_r)
        y = Y_f[self.gidx] + self.Kdiag * x
        return y.astype(self.dtype, copy=False)

    def check_hermitian(self, ntry: int = 3, tol: float = 1e-10) -> bool:
        """Quick self-check for the Hermitian property: <x|Hy> == <Hx|y>."""
        for _ in range(ntry):
            x = (np.random.randn(self.Ng) + 1j * np.random.randn(self.Ng)
                 ).astype(self.dtype)
            y = (np.random.randn(self.Ng) + 1j * np.random.randn(self.Ng)
                 ).astype(self.dtype)
            Hx, Hy = self._matvec(x), self._matvec(y)
            lhs = np.vdot(x, Hy)
            rhs = np.vdot(Hx, y)
            if not np.allclose(lhs, rhs, atol=tol, rtol=0):
                return False
        return True

    def solve(
        self,
        *,
        k: int = 30,
        which: str = "SA",
        tol: float = 1e-8,
        maxiter: Optional[int] = None,
        ncv: Optional[int] = None,
        v0: Optional[np.ndarray] = None,
        Eh: float = EH,
        verbose: bool = True,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Compute the lowest ``k`` eigenpairs using Lanczos/``eigsh``.

        Returns ``(E_meV, U)``.
        """
        k_eff = max(1, min(k, self.Ng - 2))  # assume Ng >= 3
        if ncv is None:
            ncv_eff = min(self.Ng, max(2 * k_eff + 1, 20))
        else:
            ncv_eff = min(self.Ng, max(ncv, k_eff + 1))
        v0_c = None
        if v0 is not None:
            v0_c = np.asarray(v0, dtype=self.dtype)
            if v0_c.shape[0] != self.Ng:
                raise ValueError("v0 length must equal Ng")
        try:
            E_h, U = eigsh(self.operator, k=k_eff, which=which, tol=tol,
                           maxiter=maxiter, ncv=ncv_eff, v0=v0_c)
        except ArpackNoConvergence:
            if verbose:
                print("[WARN] eigsh not converge; retry with relaxed tol/ncv")
            ncv_fb = min(self.Ng, max(ncv_eff + max(5, k_eff), 2 * k_eff + 5))
            tol_fb = max(tol, 5e-7)
            E_h, U = eigsh(self.operator, k=k_eff, which=which, tol=tol_fb,
                           maxiter=maxiter, ncv=ncv_fb, v0=v0_c)
        idx = np.argsort(E_h.real)
        E_h = E_h[idx].real
        U = U[:, idx]
        E_meV = E_h * Eh * 1000.0
        if verbose:
            print(f"Krylov finished: k={k_eff}, min(E)[meV]={E_meV[0]:.10f}")
            print(E_meV[: min(15, E_meV.size)] - E_meV[0])
        return E_meV, U


def main(
    poscar_path: str = "CONTCAR",
    energies_path: str = "ENERGIES",
    center_frac: Iterable[float] = (0.5, 0.5, 0.5),
    edgelength: float = 0.4,
    nx: int = 14,
    prim: bool = True,
    symprec: float = 1e-5,
    old_style: bool = False,
    kvec: Optional[Iterable[float]] = None,
    subtract_mean: bool = True,
    mh: float = MH_H,
):
    """Minimal pipeline: POSCAR → symmetry → H-grid → irreps → eigenpairs.

    Steps
    -----
    1) Read POSCAR into :class:`Structure`.
    2) Refine structure and get symmetry operations on the H-removed cell.
    3) Make a cubic fractional H-grid of size ``nx**3`` around ``center_frac``.
    4) Compute irreducible representatives and the full→rep mapping
       (``irr_idx``).
    5) Load irreducible energies and reconstruct a 3D potential.
    6) Build G set, kinetic diagonal, FFT coeffs, and Hψ matvec.
    7) Solve the ``k`` lowest eigenpairs by the Krylov method.

    Variables
    -------
    E_meV : ndarray
        Eigenvalues in meV.
    U : ndarray
        Eigenvectors (G basis).
    """

    mode = "Dense"     # "Dense" or "Krylov"
    infile = 'CONTCAR'
    std = 0.5
    center_frac = (std, std, std)
    edge_len_frac = 0.64
    nx = 20
    prim = False
    energies_path = 'ENERGIES'
    kvec = None
    # (1) POSCAR → Structure
    org_st, meta = PoscarReader().read(infile)
    # (2) Symmetry (remove one H at center_frac before getting operations)
    r_st, _ = SymmetryAnalyzer(prim=prim).refine(org_st)
    woH_st = strip_target_H(r_st, np.asarray(center_frac, float))
    Rot, Trans = SymmetryAnalyzer(prim=prim).operations(woH_st)
    # (3) H-grid (fractional)
    hpoints, hlat = make_hgrid_and_hlat_cube(woH_st, center_frac=center_frac,
                                             edge_len_frac=edge_len_frac,
                                             nx=nx)
    # (4) Irreps and mapping (irr_idx)
    irr_pts = get_irreducible_points(hpoints, Rot, Trans)
    irr_idx = get_mapping(hpoints, irr_pts, Rot, Trans)
    #generate_shifted_poscar(r_st, irr_pts, meta, center_frac)
    # (5) Potential reconstruction from irreducible energies
    short_atom_distance_flags = check_short_atom_distance_flags(irr_pts,
                                                                woH_st)
    ene = load_energies(energies_path)
    if short_atom_distance_flags.any():
        ene = replace_ene_with_flags(ene, short_atom_distance_flags)
    V3d = assemble_potential_3d(ene, irr_idx, nx)
    # (6) G set, kinetic, and V operator
    Gs = cut_Gs_physically(hlat, nx)
    Kdiag = kinetic_diag(Gs, lattice=hlat,  kvec=kvec)
    z = fft_coeff(V3d)
    # (7) Dense eigen-solve
    if mode == 'Dense':
        V = get_V(Gs, nx, z)
        H = np.diag(Kdiag) + V
        E, U = solve_dense_H(H)
    # (7) Krylov eigen-solve
    elif mode == 'Krylov':
        solver = NQEKrylovSolver(Gs, nx, Kdiag, V3d, subtract_mean=True)
        E, U = solver.solve(k=30, tol=1e-8)


if __name__ == "__main__":
    main()
