#!/usr/bin/env python
"""
Calculate Eigenvalues and Eigenvectors of Hydrogen nucleus quantum states
in a crystal.
- Main parts in the original change_h_pos_in_poscar_class.py are rearranged
  to be like the parts in pypolymlp codes.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Tuple, Optional, Union

from pathlib import Path 
import numpy as np
import spglib
from scipy.sparse.linalg import LinearOperator, eigsh
from scipy.linalg import eigh
# convergence exception
from scipy.sparse.linalg._eigen.arpack.arpack import (
    ArpackNoConvergence,
)

# ------------------------------ Physical constants ------------------------------
A0: float = 0.5291772  # Bohr [Å]
EH: float = 27.211  # 1 Hartree [eV]
MH_H: float = 1836.0  # proton mass [me]


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
        lat = (
            np.array(
                [list(map(float, lines[2 + i].split()[:3])) for i in range(3)]
            )
            * scale
        )
        nspc = np.asarray(lines[6].split(), dtype=int)
        numbers = np.array(
            [i + 1 for i, c in enumerate(nspc) for _ in range(c)]
        )
        if "Di" in lines[7] or "di" in lines[7]:
            iniline = 8
        elif "Sel" in lines[7] or "sel" in lines[7]:
            iniline = 9
        else:
            iniline = 8  # default fallback
        pos = np.array(
            [
                list(map(float, lines[iniline + i].split()[:3]))
                for i in range(numbers.size)
            ]
        )
        s = Structure(axis=lat, positions=pos, numbers=numbers)
        meta = {"lines": lines, "iniline": iniline}
        return s, meta


# ------------------------------ Symmetry ------------------------------
class SymmetryAnalyzer:
    """Thin spglib wrapper to obtain primitive/refined cells and
    operations.
    """

    def __init__(self, prim: bool = True, symprec: float = 1e-5):
        self.prim = prim
        self.symprec = symprec

    def refine(self, s: Structure) -> Tuple[Structure, dict]:
        """Return primitive/refined cell with positions sorted by
        species.
        """
        cell = (s.axis, s.positions, s.numbers)
        if self.prim:
            lat, pos, num = spglib.find_primitive(cell)
        else:
            lat, pos, num = spglib.refine_cell(cell)
        idx = np.argsort(num)
        s2 = Structure(
            axis=lat, positions=pos[idx], numbers=np.sort(num)
        )
        return s2, {"sorted_index": idx}

    def operations(self, s: Structure) -> Tuple[np.ndarray, np.ndarray]:
        """Return symmetry rotations and translations for the
        structure.
        """
        cell = (s.axis, s.positions, s.numbers)
        sym = spglib.get_symmetry(cell, symprec=self.symprec)
        return sym["rotations"], sym["translations"]


def strip_target_H(
    s: Structure, center_frac: Iterable[float], tol: float = 1e-5
) -> Structure:
    """Remove the H atom nearest to ``center_frac`` (fractional) from
    ``s``. Selection uses the smallest periodic distance within a
    tolerance.
    """
    c = np.asarray(center_frac, float)
    keep = np.sum(np.abs((s.positions - c) % 1.0), axis=1) >= tol
    return Structure(
        axis=s.axis, positions=s.positions[keep], numbers=s.numbers[keep]
    )


# ------------------------------ H-grid ------------------------------
def make_hgrid_and_hlat_cube(
    structure: Structure,
    center_frac: Iterable[float],
    edge_len_frac: Union[float, Iterable[float]],
    nx: int,
    *,
    old_style: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return a cubic grid around ``center_frac`` in fractional coords.

    Returns points (nx**3, 3) and an effective lattice (3, 3) for the
    cube.
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
                base = center_frac - edge_len_frac / 2
                p = base + np.array([ix, iy, iz]) * dx
                pts.append(p % 1.0)
    lat = np.asarray(structure.axis, float)
    hlat = lat @ np.diag(edge_len_frac)
    return np.asarray(pts), hlat


# ------------------------------ Irreducible representatives ---------------------
def get_irreducible_points(
    hpos: np.ndarray,
    rotations: np.ndarray,
    translations: np.ndarray,
    tol: float = 1e-5,
) -> np.ndarray:
    """Remove symmetry-equivalent grid points and return irreducible
    ones.

    Algorithm
    ---------
    - Apply each non-identity symmetry operation to the current set.
    - Drop points that become identical under that operation (with PBC).
    - Always keep the self-point at each step.
    """
    irr_hpos = np.array(hpos)
    for rot, trans in zip(rotations[1:], translations[1:]):
        sym_hpos = ((np.matmul(rot, irr_hpos.T)).T + trans) % 1.0
        sym_hpos[np.abs(sym_hpos - 1.0) <= 1.0e-14] = 0.0
        cond = np.ones((irr_hpos.shape[0]), dtype=bool)
        for i, _irr in enumerate(irr_hpos):
            if cond[i]:
                diff = np.sum(np.abs(_irr - sym_hpos) % 1.0, axis=1)
                cond *= diff >= tol
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
    """Map each full-grid point to the index of its irreducible repr."""
    test = (
        np.matmul(rotations, irr_hpos.T).transpose(2, 0, 1) + translations
    ) % 1.0
    test[np.abs(test - 1.0) <= 1.0e-15] = 0.0
    test2 = np.repeat(np.expand_dims(test, 2), hpos.shape[0], axis=2)
    cond = np.prod(
        np.sum(np.abs(hpos - test2) % 1.0, axis=3) >= tol,
        axis=1,
        dtype=bool,
    )
    pairs = np.where(np.invert(cond))
    irr_idx = pairs[0][np.argsort(pairs[1])]
    return irr_idx


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


# ------------------------------ Adiabatic potential data -----------------------
def load_energies(path: str, *, relax: bool = False) -> np.ndarray:
    """Load ``ENERGIES`` as a 1D float array.

    If ``relax=False`` use col 4 (idx 3); if True use col 5 (idx 4).
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
    """Replace energies where flag True with mean(non-flagged) + delta."""
    eave = np.average(ene[~short_atom_distance_flags])
    ene[short_atom_distance_flags] = eave + delta
    return ene


def expand_positions_3x3x3(frac_pos: np.ndarray) -> np.ndarray:
    """Expand fractional coords (M, 3) into 27 periodic images."""
    pos = np.asarray(frac_pos, float)
    shifts = np.array(
        [(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)],
        float,
    )
    return shifts[:, None, :] + pos[None, :, :]


def get_min_bond_length(
    h_frac: np.ndarray, axis: np.ndarray, host_frac: np.ndarray
) -> float:
    """Shortest host–H distance [Å] for a candidate H at ``h_frac``."""
    epos = expand_positions_3x3x3(host_frac)
    diff_frac = epos - np.asarray(h_frac, float)
    diff_cart = diff_frac @ np.asarray(axis, float)
    dists = np.linalg.norm(diff_cart, axis=2)
    return float(np.min(dists))


def check_short_atom_distance_flags(
    irr_points: np.ndarray, woH_structure: Structure, threshold: float = 0.3
) -> np.ndarray:
    """Return True where min host–H distance is below ``threshold`` [Å]."""
    axis = woH_structure.axis
    host = woH_structure.positions
    flags = [
        get_min_bond_length(p, axis, host) < threshold
        for p in np.asarray(irr_points, float)
    ]
    return np.asarray(flags, dtype=bool)


def assemble_potential_3d(
    irr_energies: np.ndarray, mapping: np.ndarray, nx: int
) -> np.ndarray:
    """Assemble (nx, nx, nx) potential array from irreducible energies."""
    return np.asarray(irr_energies, float)[mapping].reshape(nx, nx, nx)


# ==== minimal additions for three-way energies backend ====

def _reconstruct_woH_from_meta(meta: dict, center_frac) -> tuple[Structure, Structure]:
    """
    meta['lines'] と meta['iniline'] から refined 構造を復元し、
    その上で target H を strip した構造（woH）を返す。
    """
    # avoid E741: ‘l’ (lowercase L)
    lat = np.array([list(map(float, meta["lines"][2 + i].split()[:3])) for i in range(3)], float)
    pos_lines = meta["lines"][meta["iniline"]:]
    positions = np.array([list(map(float, ln.split()[:3])) for ln in pos_lines], float)
    numbers = np.arange(1, positions.shape[0] + 1, dtype=int)
    struct = Structure(axis=lat, positions=positions, numbers=numbers)
    woH = strip_target_H(struct, np.asarray(center_frac, float))
    return struct, woH


# ---- pypolymlp hooks (safe stubs; replace later) ----
def polymlp_eval_batch(polymlp_opts: dict, woH_structure: Structure, irr_pts: np.ndarray) -> np.ndarray:
    """
    TODO: 実装時に pypolymlp を呼ぶ。
    ひとまず ‘未実装’ エラーで明確に落とす（VASP ルートはこの関数を通らない）。
    """
    raise NotImplementedError("polymlp_eval_batch is not implemented yet. Provide model & eval logic.")


def polymlp_relax(polymlp_opts: dict, woH_structure: Structure, irr_pts: np.ndarray, nx: int) -> np.ndarray:
    """
    TODO: 実装時に pypolymlp + forces + L-BFGS-B を呼ぶ。
    """
    raise NotImplementedError("polymlp_relax is not implemented yet. Provide relax logic.")


def get_energies(
    backend: str,
    *,
    irr_pts: np.ndarray,
    structure: Structure,          # refined 構造（r_st）
    woH_st: Structure,             # H を抜いた構造
    center_frac: tuple[float,float,float],
    nx: int,
    energies_path: str = "ENERGIES",
    polymlp_opts: dict | None = None,
    meta: dict | None = None,      # VASP で POSCAR 出力が必要なときだけ渡す
) -> np.ndarray | None:
    """
    Return:
      np.ndarray(shape=(N_irrep,))  … energy array
      None … VASP initial run (POSCARs are generated and end temporaly）
    """
    if backend == "vasp":
        p = Path(energies_path)
        if p.exists():
            return load_energies(str(p))
        if meta is None:
            raise ValueError("VASP ルートで POSCAR を出すには meta が必要です。")
        generate_shifted_poscar(structure, irr_pts, meta, center_frac)
        print("[INFO] POSCAR_### を生成。VASP 実行後に同じコマンドで再実行してください。")
        return None

    if backend == "polymlp-batch":
        # TODO: pypolymlp 実装を後で差し込み
        return polymlp_eval_batch(polymlp_opts or {}, woH_st, irr_pts)

    if backend == "polymlp-relax":
        # TODO: pypolymlp+forces の緩和実装を後で差し込み
        return polymlp_relax(polymlp_opts or {}, woH_st, irr_pts, nx)

    raise ValueError("backend must be 'polymlp-batch' | 'polymlp-relax' | 'vasp'")


# ----------------------------- G-space utilities -----------------------------
def cut_Gs_physically(lattice: np.ndarray, nx: int, ao: float = A0
                      ) -> np.ndarray:
    """Generate G-vectors by a physical cutoff in reciprocal space."""
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


# ------------------------------ V operator (FFT) ------------------------------
def fft_coeff(
    V3d: np.ndarray, *, Eh: float = EH, subtract_mean: bool = True
) -> np.ndarray:
    """Return FFT coefficients of ``V/Eh`` with ``norm='forward'``."""
    arr = (V3d - np.mean(V3d)) / Eh if subtract_mean else (V3d / Eh)
    return np.fft.fftn(arr, norm="forward")


def get_V(Gs: np.ndarray, nx: int, z: np.ndarray) -> np.ndarray:
    """Build the potential matrix ``V`` in the G basis via FFT indexing."""
    gi, gj, gk = Gs[:, 0], Gs[:, 1], Gs[:, 2]
    di = (gi[:, None] - gi[None, :]) % nx
    dj = (gj[:, None] - gj[None, :]) % nx
    dk = (gk[:, None] - gk[None, :]) % nx
    return z[di, dj, dk]


# ------------------------------ Krylov solver ------------------------------
class NQEKrylovSolver:
    """Small utility class to solve ``Hψ`` with a Krylov (``eigsh``)
    method.
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
    ) -> None:
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
        self.operator = LinearOperator(
            shape=(self.Ng, self.Ng), matvec=self._matvec, dtype=dtype
        )

    def _matvec(self, x: np.ndarray) -> np.ndarray:
        X = self._X
        X.fill(0.0)
        X[self.gidx] = x.astype(self.dtype, copy=False)
        phi_r = np.fft.ifftn(X)
        Y_r = self.V_r * phi_r
        Y_f = np.fft.fftn(Y_r)
        y = Y_f[self.gidx] + self.Kdiag * x
        return y.astype(self.dtype, copy=False)

    def check_hermitian(self, ntry: int = 3, tol: float = 1.0e-10) -> bool:
        """Quick self-check for Hermitian: <x|Hy> == <Hx|y>."""
        for _ in range(ntry):
            rnd = np.random.randn
            x = (rnd(self.Ng) + 1j * rnd(self.Ng)).astype(self.dtype)
            y = (rnd(self.Ng) + 1j * rnd(self.Ng)).astype(self.dtype)
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
        tol: float = 1.0e-8,
        maxiter: Optional[int] = None,
        ncv: Optional[int] = None,
        v0: Optional[np.ndarray] = None,
        Eh: float = EH,
        verbose: bool = True,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Compute the lowest ``k`` eigenpairs using ``eigsh``."""
        k_eff = max(1, min(k, self.Ng - 2))
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
            E_h, U = eigsh(
                self.operator,
                k=k_eff,
                which=which,
                tol=tol,
                maxiter=maxiter,
                ncv=ncv_eff,
                v0=v0_c,
            )
        except ArpackNoConvergence:
            if verbose:
                print("[WARN] eigsh not converge; retry with relaxed tol/ncv")
            ncv_fb = min(self.Ng, max(ncv_eff + max(5, k_eff), 2 * k_eff + 5))
            tol_fb = max(tol, 5e-7)
            E_h, U = eigsh(
                self.operator,
                k=k_eff,
                which=which,
                tol=tol_fb,
                maxiter=maxiter,
                ncv=ncv_fb,
                v0=v0_c,
            )
        idx = np.argsort(E_h.real)
        E_h = E_h[idx].real
        U = U[:, idx]
        E_meV = E_h * Eh * 1000.0
        if verbose:
            msg = "Krylov finished: k={k}, min(E)[meV]={e:.10f}"
            print(msg.format(k=k_eff, e=E_meV[0]))
            print(E_meV[: min(15, E_meV.size)] - E_meV[0])
        return E_meV, U


# ------------------------------ Solvers (Dense path) ---------------------------
def solve_dense_H(
    H: np.ndarray, Compress: bool = False, v: Optional[np.ndarray] = None,
    verbose: bool = True
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute eigenpairs of a dense Hermitian matrix using ``eigh``."""
    if Compress and v is None:
        if verbose:
            print("Compression is used in Eigenvalue solution")
            print("Calculating v matrix of 200 eigen vectors..")
        dummy, v = eigh(H, subset_by_index=[0, 199], driver="evr")
        dummy *= EH * 1000.0
        if verbose:
            print(np.min(dummy))
            print(dummy[0:15] - np.min(dummy))
    if Compress and v is not None:
        if verbose:
            print("Compressing H by v.")
        H_comp = v.conj().T @ H @ v
        if verbose:
            print("Calculating 30 eigen from the compressed H")
            E_comp, U_comp = eigh(H_comp, subset_by_index=[0, 29], driver="evr")
            E_comp *= EH * 1000.0
        if verbose:
            print(np.min(E_comp))
            print(E_comp[0:15] - np.min(E_comp))
        U = v @ U_comp
        return E_comp, U
    if not Compress:
        E, U = eigh(H, subset_by_index=[0, 29])
        E *= EH * 1000.0
        if verbose:
            print("No compression is used in Eigenvalue solution")
            print(np.min(E))
            print(E[0:15] - np.min(E))
        return E, U


# ---------------------------- CLI entry (original flow preserved) ------------
def main(
    poscar_path: str = "CONTCAR",
    energies_path: str = "ENERGIES",
    center_frac: Iterable[float] = (0.5, 0.5, 0.5),
    edgelength: float = 0.4,
    nx: int = 14,
    prim: bool = True,
    symprec: float = 1.0e-5,
    old_style: bool = False,
    kvec: Optional[Iterable[float]] = None,
    subtract_mean: bool = True,
    mh: float = MH_H,
) -> Tuple[np.ndarray, np.ndarray]:
    """Minimal pipeline: POSCAR → symmetry → H-grid → irreps → eigenpairs.
       Note that thsi assumes total energies of vasp are stored in
       energies_path.
       If no file of the total energies, execute (1), (2), (3), (4), and (4)',
       prepare the file, execute (1), (2), (3), (4), (5), (6) and (7).
    """
    # --------- IMPORTANT: keep exactly as in your original file ----------
    mode = "Dense"  # "Dense" or "Krylov"
    infile = "CONTCAR"
    std = 0.5
    center_frac = (std, std, std)
    edge_len_frac = 0.64
    nx = 20
    prim = False
    energies_path = "ENERGIES"
    kvec = None
    # --------------------------------------------------------------------

    # (1) POSCAR → Structure
    org_st, meta = PoscarReader().read(infile)

    # (2) Symmetry (remove one H at center_frac before operations)
    r_st, _ = SymmetryAnalyzer(prim=prim).refine(org_st)
    woH_st = strip_target_H(r_st, np.asarray(center_frac, float))
    Rot, Trans = SymmetryAnalyzer(prim=prim).operations(woH_st)

    # (3) H-grid (fractional)
    hpoints, hlat = make_hgrid_and_hlat_cube(
        woH_st,
        center_frac=center_frac,
        edge_len_frac=edge_len_frac,
        nx=nx,
    )

    # (4) Irreps and mapping (irr_idx)
    irr_pts = get_irreducible_points(hpoints, Rot, Trans)
    irr_idx = get_mapping(hpoints, irr_pts, Rot, Trans)

    # (4)' POSCARS_xxxx if you do not have a file of total energies
    #      You should generate POSCARS_xxxx and do happy vasping first.
    #      generate_shifted_poscar(r_st, irr_pts, meta, center_frac)

    # (5) Potential reconstruction from irreducible energies
    # (4) Irreps は既にあると仮定: irr_pts, irr_idx
    # energies を取得

    ene = get_energies(
        backend="vasp",                   # "polymlp-batch" / "polymlp-relax" に切替可
        irr_pts=irr_pts,
        structure=r_st,
        woH_st=woH_st,
        center_frac=center_frac,
        nx=nx,
        energies_path=energies_path,
        polymlp_opts=None,                # 後で pypolymlp 用 dict を渡す
        meta=meta                         # VASP の POSCAR 出力にだけ使用
    )

    #ene = get_energies(
    #    backend="vasp",                      # ← "polymlp-batch" / "polymlp-relax" も可
    #    irr_pts=irr_pts,
    #    meta=meta,                           # PoscarReader.read の戻りの meta を渡す
    #    center_frac=center_frac,
    #    nx=nx,
    #    energies_path="ENERGIES",
    #    polymlp_opts={"model": "path/to/model.yaml"}  # 実装時に参照
    #)
    if ene is None:
        return  # ここで一旦終了（VASP を回してから同じコマンドを再実行）

    short_atom_distance_flags = check_short_atom_distance_flags(
        irr_pts, woH_st
    )
    #ene = load_energies(energies_path)
    if short_atom_distance_flags.any():
        ene = replace_ene_with_flags(ene, short_atom_distance_flags)
    V3d = assemble_potential_3d(ene, irr_idx, nx)

    # (6) G set, kinetic, and V operator
    Gs = cut_Gs_physically(hlat, nx)
    Kdiag = kinetic_diag(Gs, lattice=hlat, kvec=kvec)
    z = fft_coeff(V3d)

    # (7) Dense eigen-solve / Krylov eigen-solve (original branching)
    if mode == "Dense":
        V = get_V(Gs, nx, z)
        H = np.diag(Kdiag) + V
        E, U = solve_dense_H(H)
    elif mode == "Krylov":
        solver = NQEKrylovSolver(Gs, nx, Kdiag, V3d, subtract_mean=True)
        E, U = solver.solve(k=30, tol=1.0e-8)
    else:
        raise ValueError("mode must be 'Dense' or 'Krylov'")
    return E, U


if __name__ == "__main__":  # pragma: no cover
    e, _u = main()
    print("min(E) [meV]:", float(np.min(e)))
