# structure_compare.py
"""
Calculate eigenvalues and eigenvectors of hydrogen nucleus quantum states
in a crystal within a single particle picture under an adiabatic potential.
- Main parts in the original change_h_pos_in_poscar_class.py are rearranged
  to be like the parts in pypolymlp codes.
"""
from __future__ import annotations
from typing import Dict, Tuple, Optional
import numpy as np
from collections import Counter, defaultdict
from itertools import permutations

# =============================
# A-1/A-2/B/C: consistency checks
# =============================

try:
    # Optional speed-up for position matching
    from scipy.spatial import cKDTree as KDTree
    _HAS_KDTREE = True
except Exception:  # pragma: no cover
    _HAS_KDTREE = False


def _cell_volume(a_mat: np.ndarray) -> float:
    """Return absolute unit cell volume from a (3, 3) row-lattice matrix."""
    return float(abs(np.linalg.det(a_mat)))


def _orthogonal_procrustes_rows(
    a_rows: np.ndarray,
    b_rows: np.ndarray,
) -> Tuple[np.ndarray, float]:
    """
    Solve min || a_rows @ R - b_rows ||_F  s.t.  R^T R = I, det(R) = +1.
    Returns rotation R (3x3) and Frobenius error.
    """
    u, _, vt = np.linalg.svd(a_rows.T @ b_rows)
    r_mat = u @ vt
    if np.linalg.det(r_mat) < 0.0:
        # Enforce a proper rotation (no inversion)
        u[:, -1] *= -1.0
        r_mat = u @ vt
    err = float(np.linalg.norm(a_rows @ r_mat - b_rows, ord="fro"))
    return r_mat, err


def cells_equivalent_up_to_perm_rot(
    a0: np.ndarray,
    a1: np.ndarray,
    tol: float = 1.0e-8,
) -> Tuple[bool, Optional[Tuple[int, int, int]], Optional[np.ndarray], float]:
    """
    Check if two row-lattice matrices are equivalent up to
    (row) axis permutation and a proper rotation.

    Returns:
        ok: True if equivalent
        perm: best row permutation (tuple of 3) or None
        r_best: best rotation (3x3) or None
        err: Frobenius norm of the best match
    """
    a0 = np.asarray(a0, float)
    a1 = np.asarray(a1, float)

    best_ok = False
    best_perm: Optional[Tuple[int, int, int]] = None
    best_r: Optional[np.ndarray] = None
    best_err = np.inf

    for p in permutations((0, 1, 2)):
        a0p = a0[list(p), :]
        r_mat, err = _orthogonal_procrustes_rows(a0p, a1)
        if err < best_err:
            best_ok = bool(err <= tol)
            best_perm = (int(p[0]), int(p[1]), int(p[2]))
            best_r = r_mat
            best_err = err
        if best_ok:
            break

    return best_ok, best_perm, best_r, float(best_err)


def _pbc_match_cartesian(
    new_cart: np.ndarray,
    old_cart: np.ndarray,
    a_mat: np.ndarray,
    tol_pos: float,
    require_one_to_one: bool,
) -> Tuple[bool, Optional[np.ndarray], str]:
    """
    PBC nearest-neighbor matching in Cartesian space.

    Args:
        new_cart: (N, 3) new Cartesian coordinates.
        old_cart: (N, 3) old Cartesian coordinates.
        a_mat: (3, 3) row-lattice matrix.
        tol_pos: tolerance in Å for nearest distance.
        require_one_to_one: if True, enforce a 1:1 mapping.

    Returns:
        ok, mapping (new_idx -> old_idx), message
    """
    n_atoms = new_cart.shape[0]

    if _HAS_KDTREE:
        shifts = np.array(
            [[i, j, k] for i in (-1, 0, 1)
             for j in (-1, 0, 1) for k in (-1, 0, 1)],
            float,
        )
        images = (old_cart[None, :, :] + shifts[:, None, :] @ a_mat).reshape(
            -1, 3
        )
        tree = KDTree(images)
        dmin, idx = tree.query(new_cart, k=1)
        if np.any(dmin > tol_pos):
            msg = (
                f"position mismatch (max d={np.max(dmin):.3e} Å "
                f"> tol {tol_pos:.1e})"
            )
            return False, None, msg
        j_old = (idx % n_atoms).astype(int)
        if require_one_to_one:
            counts = np.bincount(j_old, minlength=n_atoms)
            if np.any(counts > 1):
                return False, None, "non-1:1 mapping (degeneracy)"
        return True, j_old, "ok"

    # Fallback: O(N^2 * 27)
    used = np.zeros(n_atoms, dtype=bool)
    j_old = np.empty(n_atoms, dtype=int)
    for i in range(n_atoms):
        best_j = -1
        best_d = 1.0e9
        for j in range(n_atoms):
            for si in (-1, 0, 1):
                for sj in (-1, 0, 1):
                    for sk in (-1, 0, 1):
                        diff = new_cart[i] - (
                            old_cart[j]
                            + si * a_mat[0]
                            + sj * a_mat[1]
                            + sk * a_mat[2]
                        )
                        dval = float(np.linalg.norm(diff))
                        if dval < best_d:
                            best_d = dval
                            best_j = j
        if best_d > tol_pos:
            msg = (
                f"position mismatch (d={best_d:.3e} Å "
                f"> tol {tol_pos:.1e})"
            )
            return False, None, msg
        if require_one_to_one and used[best_j]:
            return False, None, "non-1:1 mapping (degeneracy)"
        j_old[i] = best_j
        used[best_j] = True
    return True, j_old, "ok(slow)"


def _blockwise_match_by_elements(
    new_frac: np.ndarray,
    old_frac: np.ndarray,
    a_new: np.ndarray,
    elements_new,
    elements_old,
    tol_pos: float,
    require_one_to_one: bool,
) -> Tuple[bool, Optional[np.ndarray], str]:
    """
    Match by element buckets to avoid Pd<->H mis-assignments.
    """
    mapping = np.empty(len(new_frac), dtype=int)
    new_cart = np.asarray(new_frac, float) @ a_new
    old_cart = np.asarray(old_frac, float) @ a_new

    buckets_new: Dict[object, list] = defaultdict(list)
    buckets_old: Dict[object, list] = defaultdict(list)
    for i, elem in enumerate(elements_new):
        buckets_new[elem].append(i)
    for j, elem in enumerate(elements_old):
        buckets_old[elem].append(j)

    for elem in buckets_new.keys():
        idx_new = np.array(buckets_new[elem], dtype=int)
        idx_old = np.array(buckets_old[elem], dtype=int)
        ok, part, msg = _pbc_match_cartesian(
            new_cart[idx_new],
            old_cart[idx_old],
            a_new,
            tol_pos,
            require_one_to_one,
        )
        if not ok:
            return False, None, f"[{elem}] {msg}"
        mapping[idx_new] = idx_old[part]
    return True, mapping, "ok(elements-split)"


def check_structure_consistency_full(  # noqa: D401
    old_st,
    new_st,
    *,
    # ---- A: cell equivalence ----
    tol_vol: float = 1.0e-8,
    tol_permrot: float = 1.0e-8,
    # ---- B: counts / composition ----
    use_elements: bool = True,
    # ---- C: PBC matching ----
    tol_pos: float = 1.0e-4,
    by_elements: bool = True,
    require_one_to_one: bool = True,
) -> Tuple[bool, Dict[str, object]]:
    """
    Check A-1 (volume), A-2 (perm+rotation), B (counts), C (PBC mapping).

    Returns:
        ok: True if all checks passed
        rep: dict with diagnostics:
             'msg', 'vol_ratio', 'perm', 'R', 'permrot_err',
             'map_new_to_old', 'map_old_to_new', 'type_map'
    """
    rep: Dict[str, object] = {
        "vol_ratio": np.nan,
        "perm": None,
        "R": None,
        "permrot_err": np.nan,
        "map_new_to_old": None,
        "map_old_to_new": None,
        "type_map": None,
    }

    a0 = np.asarray(old_st.axis, float)
    a1 = np.asarray(new_st.axis, float)

    # ---- A-1: volume ratio ----
    v0 = _cell_volume(a0)
    v1 = _cell_volume(a1)
    rep["vol_ratio"] = v1 / v0 if v0 != 0.0 else np.inf
    if abs(rep["vol_ratio"] - 1.0) > tol_vol:
        return False, {**rep, "msg": "cell volume changed (primitive/basis?)"}

    # ---- A-2: axis permutation + rotation equivalence ----
    ok_eq, perm, r_best, err = cells_equivalent_up_to_perm_rot(
        a0, a1, tol=tol_permrot
    )
    rep["perm"] = perm
    rep["R"] = r_best
    rep["permrot_err"] = float(err)
    if not ok_eq:
        msg = (
            "cell differs (not equivalent up to axis-permutation + rotation)"
        )
        return False, {**rep, "msg": msg}

    # ---- B: counts / composition ----
    if len(old_st.positions) != len(new_st.positions):
        return False, {**rep, "msg": "number of atoms changed"}
    if use_elements and hasattr(old_st, "elements") and hasattr(
        new_st, "elements"
    ):
        if Counter(old_st.elements) != Counter(new_st.elements):
            return False, {**rep, "msg": "composition differs (elements)"}
    else:
        # Fallback: only total count via types (labels may change)
        if len(old_st.types) != len(new_st.types):
            return False, {**rep, "msg": "type counts mismatch"}

    # ---- C: PBC nearest-neighbor mapping ----
    # Map old fractional into the *new* lattice frame:
    #   f_old_in_new = f_old @ A0 @ inv(A1), then modulo 1
    t_mat = a0 @ np.linalg.inv(a1)
    old_frac_in_new = (np.asarray(old_st.positions, float) @ t_mat) % 1.0
    new_frac = np.asarray(new_st.positions, float)

    if by_elements and hasattr(old_st, "elements") and hasattr(
        new_st, "elements"
    ):
        ok, map_no, msg = _blockwise_match_by_elements(
            new_frac,
            old_frac_in_new,
            a1,
            new_st.elements,
            old_st.elements,
            tol_pos,
            require_one_to_one,
        )
    else:
        ok, map_no, msg = _pbc_match_cartesian(
            new_frac @ a1,
            old_frac_in_new @ a1,
            a1,
            tol_pos,
            require_one_to_one,
        )

    if not ok:
        return False, {**rep, "msg": msg}

    rep["map_new_to_old"] = map_no.astype(int)
    inv = np.empty_like(map_no)
    inv[map_no] = np.arange(len(map_no), dtype=int)
    rep["map_old_to_new"] = inv.astype(int)

    # Optional: majority-vote new_type -> old_type map
    new_types = np.asarray(new_st.types, int)
    old_types = np.asarray(old_st.types, int)
    type_map: Dict[int, int] = {}
    for ntype in np.unique(new_types):
        idx_n = np.where(new_types == ntype)[0]
        olds = old_types[rep["map_new_to_old"][idx_n]]
        vals, counts = np.unique(olds, return_counts=True)
        type_map[int(ntype)] = int(vals[np.argmax(counts)])
    rep["type_map"] = type_map

    return True, {**rep, "msg": "consistent"}


def reorder_like(
    old_st,
    new_st,
    map_old_to_new: np.ndarray,
    type_map: Optional[Dict[int, int]] = None,
):
    """
    Return a copy of new_st reordered to match old_st atom order.
    Optionally remap type labels with type_map (new_type -> old_type).
    """
    idx = np.asarray(map_old_to_new, int)
    pos2 = np.asarray(new_st.positions, float)[idx]
    types2 = np.asarray(new_st.types, int)[idx]
    if type_map is not None:
        # Vectorized dictionary lookup with fallback identity
        def _map_t(t: int) -> int:
            return int(type_map.get(int(t), int(t)))
        types2 = np.vectorize(_map_t)(types2)
    if hasattr(new_st, "elements"):
        elems2 = np.asarray(new_st.elements, object)[idx].tolist()
    else:
        elems2 = []

    return Structure(
        axis=new_st.axis,
        positions=pos2,
        types=types2.tolist(),
        n_atoms=new_st.n_atoms,
        elements=elems2,
        name=new_st.name,
        comment=new_st.comment,
        volume=new_st.volume,
    )
