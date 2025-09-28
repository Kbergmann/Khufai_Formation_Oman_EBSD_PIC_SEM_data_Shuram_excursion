#!/usr/bin/env python3
"""
XRD (.xrdml) peak extraction & dolomite ordering metrics

Features
- Robust XRDML parsing (PANalytical)
- Background correction: median (rolling) or ASLS
- 2θ correction: per-scan quartz(101) or manual constant offset
- Peak detection: tallest local max in window (fallback tallest point)
- FWHM estimation
- Metrics: I(015)/I(104), I(021)/I(104), I(101)/I(104), I(015)/I(110), FWHM(015)/FWHM(104|110)
- Single CSV output (one row per scan)
- Plots:
  * Grey histogram I(015)/I(104) with colored label points & legend
  * Scatter FWHM(015)/FWHM(104) vs I(015)/I(104) with labeled points
  * Box plot I(015)/I(104) & I(015)/I(110) with colored overlays
  * Grey histogram of applied 2θ offsets with labeled points
  * Cross-plot I(015)/I(104) vs I(015)/I(110)
  * Per-scan corrected-only spectrum PNG with seed bars (grey/orange/gold)
  * Per-scan 5-panel QC (Quartz101, Dol(104),(110),(015),(021)) with dashed seed & red used peak.

Label file
- One per line: "<filename_or_path>,<label>[,<hexcolor>]" (comma, no spaces)
- Colors like "#1f77b4" or "1f77b4" (hex optional). If omitted, auto colormap is used.
"""

import argparse
import glob
import os
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET

# ---------- Reference positions (Cu Kα) ----------
# Quartz from RUFF R040031
QUARTZ_REF: Dict[str, float] = {
    "(100)": 20.8943,
    "(101)": 26.6737,
    "(110)": 36.5833,
}

# Dolomite (R-3c). Ordering-sensitive: (015), (021), (dolomite 101) from R040030 RUFF
DOLOMITE_REF: Dict[str, float] = {
    "(101)": 22.08,  # superstructure (ordering; very weak in protodolomite; adjust if needed)
    "(012)": 24.1186,
    "(104)": 31.0073,  # fundamental
    "(110)": 37.4187,  # fundamental
    "(006)": 33.6015,
    "(015)": 35.3787,  # superstructure (ordering); tweak if your lab uses a slightly different value
    "(113)": 41.1963,  # fundamental
    "(021)": 43.8568,  # superstructure (ordering; weaker)
    "(202)": 44.9990,
    "(214)": 63.5073,
    "(300)": 67.4643,
    "(0012)": 70.58,
}

# Calcite (optional flagging only; not used for exclusion)
CALCITE_REF: Dict[str, float] = {
    "(012)": 23.04,
    "(104)": 29.42,  # strong
    "(006)": 35.97,
    "(110)": 39.42,
    "(113)": 43.14,
}

# ---------- Parsing ----------
def parse_xrdml_robust(path: str) -> Tuple[np.ndarray, np.ndarray]:
    tree = ET.parse(path)
    root = tree.getroot()
    ns = {"xrd": root.tag.split("}")[0].strip("{")}

    positions_nodes = root.findall(".//xrd:scan/xrd:dataPoints/xrd:positions", ns)
    intensities_node = root.find(".//xrd:scan/xrd:dataPoints/xrd:intensities", ns)
    if intensities_node is None or not intensities_node.text:
        raise ValueError(f"{path}: no <intensities> data")
    intensities = np.fromstring(intensities_node.text.strip(), sep=" ").astype(float)

    pos2t = None
    for n in positions_nodes:
        if n.attrib.get("axis", "").lower() in ("2theta", "twotheta"):
            pos2t = n
            break
    if pos2t is None:
        pos2t = positions_nodes[0] if positions_nodes else None
    if pos2t is None:
        raise ValueError(f"{path}: no <positions> node")

    start = pos2t.get("startPosition")
    end = pos2t.get("endPosition")
    if start is None:
        sp = pos2t.find("xrd:startPosition", ns)
        start = sp.text if sp is not None else None
    if end is None:
        ep = pos2t.find("xrd:endPosition", ns)
        end = ep.text if ep is not None else None

    if start is not None and end is not None:
        two_theta = np.linspace(float(start), float(end), intensities.size)
    else:
        # explicit positions?
        if pos2t.text and pos2t.text.strip():
            two_theta = np.fromstring(pos2t.text.strip(), sep=" ").astype(float)
            n = min(two_theta.size, intensities.size)
            two_theta = two_theta[:n]
            intensities = intensities[:n]
        else:
            raise ValueError(f"{path}: positions missing")

    return two_theta, intensities

# ---------- Peak utilities ----------
def simple_peak_find(x: np.ndarray,
                     y: np.ndarray,
                     min_rel_height: float = 0.02,
                     min_distance_pts: int = 5) -> List[int]:
    if y.size < 3:
        return []
    y_max = float(np.max(y))
    thresh = y_max * float(min_rel_height)
    cand = np.where((y[1:-1] > y[:-2]) & (y[1:-1] > y[2:]) & (y[1:-1] >= thresh))[0] + 1
    if cand.size == 0:
        return []
    selected: List[int] = []
    for idx in cand:
        if selected and (idx - selected[-1]) < int(min_distance_pts):
            if y[idx] > y[selected[-1]]:
                selected[-1] = idx
        else:
            selected.append(idx)
    return selected

def estimate_fwhm(x: np.ndarray, y: np.ndarray, peak_idx: int) -> Optional[float]:
    if peak_idx <= 0 or peak_idx >= len(x) - 1:
        return None
    half = y[peak_idx] / 2.0
    # left
    left = None
    for i in range(peak_idx - 1, 0, -1):
        if (y[i] - half) * (y[i + 1] - half) <= 0:
            x1, x2 = x[i], x[i + 1]
            y1, y2 = y[i], y[i + 1]
            left = x1 + (half - y1) * (x2 - x1) / (y2 - y1 + 1e-12)
            break
    # right
    right = None
    for i in range(peak_idx, len(x) - 1):
        if (y[i] - half) * (y[i + 1] - half) <= 0:
            x1, x2 = x[i], x[i + 1]
            y1, y2 = y[i], y[i + 1]
            right = x1 + (half - y1) * (x2 - x1) / (y2 - y1 + 1e-12)
            break
    if left is None or right is None:
        return None
    return float(right - left)

def match_targets(two_theta: np.ndarray,
                  intensities: np.ndarray,
                  peak_indices: List[int],
                  targets: Dict[str, float],
                  tolerance: float = 0.30) -> Dict[str, Dict[str, Optional[float]]]:
    """
    For each target hkl, pick the *tallest* detected local max within ±tolerance;
    if none, fallback to tallest raw point within the window.
    """
    results: Dict[str, Dict[str, Optional[float]]] = {}
    peaks_tt = two_theta[peak_indices] if len(peak_indices) else np.array([])
    peaks_ht = intensities[peak_indices] if len(peak_indices) else np.array([])

    for hkl, tt0 in targets.items():
        win = (two_theta >= tt0 - tolerance) & (two_theta <= tt0 + tolerance)
        if not np.any(win):
            results[hkl] = {"2theta": None, "height": None, "fwhm": None}
            continue

        if peaks_tt.size:
            in_win_idx = np.where((peaks_tt >= tt0 - tolerance) & (peaks_tt <= tt0 + tolerance))[0]
        else:
            in_win_idx = np.array([], dtype=int)

        if in_win_idx.size > 0:
            best = in_win_idx[np.argmax(peaks_ht[in_win_idx])]
            idx = peak_indices[best]
        else:
            idxs = np.where(win)[0]
            idx = idxs[np.argmax(intensities[idxs])]  # tallest point

        fwhm = estimate_fwhm(two_theta, intensities, idx)
        results[hkl] = {
            "2theta": float(two_theta[idx]),
            "height": float(intensities[idx]),
            "fwhm": float(fwhm) if fwhm is not None else None,
        }
    return results

def expand_matches_as_row(matches: Dict[str, Dict[str, Optional[float]]],
                          hkls: List[str]) -> Dict[str, Optional[float]]:
    row: Dict[str, Optional[float]] = {}
    for h in hkls:
        m = matches.get(h, {})
        tag = h.strip("()")
        row[f"tt_{tag}"] = m.get("2theta")
        row[f"I_{tag}"] = m.get("height")
        row[f"FWHM_{tag}"] = m.get("fwhm")
    return row

def compute_dolomite_metrics(matches: Dict[str, Dict[str, Optional[float]]]) -> Dict[str, Optional[float]]:
    def gI(h):
        v = matches.get(h, {}).get("height")
        return None if v is None else float(v)
    def gW(h):
        v = matches.get(h, {}).get("fwhm")
        return None if v is None else float(v)
    def gTT(h):
        v = matches.get(h, {}).get("2theta")
        return None if v is None else float(v)
    def safe(a, b):
        return (a / b) if (a is not None and b not in (None, 0)) else None

    # Intensities we care about
    I015 = gI("(015)")
    I021 = gI("(021)")
    I110 = gI("(110)")

    # Widths/positions we care about
    W015 = gW("(015)")
    tt104 = gTT("(104)")

    return {
        # NEW: ratios over (110)
        "I015_over_I110": safe(I015, I110),
        "I021_over_I110": safe(I021, I110),

        # NEW: half width of 015
        "HWHM_015": (W015 / 2.0) if W015 is not None else None,

        # NEW: 104 position offset relative to seed (stoichiometry proxy)
        "delta_tt_104": (tt104 - DOLOMITE_REF["(104)"]) if tt104 is not None else None,
    }

# ---------- Background correction ----------
def apply_rolling_median_baseline(two_theta: np.ndarray,
                                  intensities: np.ndarray,
                                  bg_window_deg: float) -> np.ndarray:
    # estimate step size robustly
    step = np.median(np.diff(two_theta)) if len(two_theta) > 1 else bg_window_deg
    win_pts = max(3, int(round(bg_window_deg / max(step, 1e-9))))
    if win_pts % 2 == 0:
        win_pts += 1
    s = pd.Series(intensities)
    baseline = s.rolling(window=win_pts, center=True, min_periods=1).median().to_numpy()
    y = intensities - baseline
    y[y < 0] = 0
    return y

def baseline_asls(y: np.ndarray, lam: float = 1e5, p: float = 0.01, niter: int = 10) -> Tuple[np.ndarray, np.ndarray]:
    """
    Asymmetric least squares baseline (Eilers & Boelens).
    Dense implementation (can be memory heavy for very long scans).
    """
    L = len(y)
    D = np.diff(np.eye(L), 2, axis=0)     # (L-2, L)
    w = np.ones(L)
    for _ in range(niter):
        W = np.diag(w)                    # (L, L)
        Z = W + lam * (D.T @ D)           # (L, L)
        z = np.linalg.solve(Z, w * y)     # baseline
        w = p * (y > z) + (1 - p) * (y < z)
    b = z
    yc = y - b
    yc[yc < 0] = 0
    return yc, b

# ---------- Quartz offset helpers ----------
def find_peak_near(two_theta: np.ndarray, intensity: np.ndarray,
                   target_tt: float, search_window: float = 0.8,
                   min_rel_height: float = 0.02, min_distance_pts: int = 5) -> Optional[Tuple[float, float]]:
    mask = (two_theta >= target_tt - search_window) & (two_theta <= target_tt + search_window)
    if not np.any(mask):
        return None
    x = two_theta[mask]; y = intensity[mask]
    if y.size < 3:
        return None
    # local peaks
    cand = np.where((y[1:-1] > y[:-2]) & (y[1:-1] > y[2:]))[0] + 1
    if cand.size == 0:
        return None
    tallest = cand[np.argmax(y[cand])]
    return float(x[tallest]), float(y[tallest])

def quartz_offset_for_scan(two_theta: np.ndarray, intensity: np.ndarray,
                           quartz_main_tt: float = 26.6737,
                           search_window: float = 0.8) -> Optional[float]:
    found = find_peak_near(two_theta, intensity, quartz_main_tt, search_window=search_window)
    if found is None:
        return None
    measured_tt, _ = found
    return float(measured_tt - quartz_main_tt)

# ---------- Per-scan plotters ----------
def save_corrected_spectrum_plot(path: str,
                                 two_theta: np.ndarray,
                                 intensity_corr: np.ndarray,
                                 spectra_outdir: str,
                                 quartz_main_tt: float,
                                 dolo_refs: Dict[str, float]):
    """Corrected-only spectrum with seed bars (grey/orange/gold) BEHIND the trace."""
    os.makedirs(spectra_outdir, exist_ok=True)
    base = os.path.splitext(os.path.basename(path))[0]

    fig, ax = plt.subplots(figsize=(16, 3))

    # --- reference bars (low zorder so they sit behind the spectrum) ---
    ordering = {"(015)", "(021)", "(101)"}
    for h, tt in sorted(dolo_refs.items(), key=lambda kv: kv[1]):
        color = "orange" if h in ordering else "gold"
        ax.axvline(tt, linestyle="--", color=color, linewidth=1.0, zorder=1,
                   label=("Dol (ordering)" if h == "(015)" else None))
    ax.axvline(quartz_main_tt, linestyle="--", color="grey", linewidth=1.0, zorder=1,
               label="Quartz (101)")

    # --- spectrum on top ---
    ax.plot(two_theta, intensity_corr, lw=1, zorder=5, label="baseline-corrected")

    # legend (dedupe)
    h, l = ax.get_legend_handles_labels()
    by = dict(zip(l, h))
    if by:
        ax.legend(by.values(), by.keys(), fontsize=8, loc="best")

    ax.set_xlabel("2θ ($^\circ$)")
    ax.set_ylabel("Intensity (corrected)")
    ax.set_title(f"Corrected spectrum: {base}")
    fig.tight_layout()
    out_png = os.path.join(spectra_outdir, f"{base}.png")
    fig.savefig(out_png, dpi=180)
    plt.close(fig)

def save_peak_check_panels(path: str,
                           two_theta: np.ndarray,
                           intensity_corr: np.ndarray,
                           window_half_deg: float,
                           dolo_refs: Dict[str, float],
                           quartz_refs: Dict[str, float],
                           outdir: str):
    """1x5 row: Quartz(101), Dol(104),(110),(015),(021) with dashed seed & red used peak."""
    os.makedirs(outdir, exist_ok=True)
    base = os.path.splitext(os.path.basename(path))[0]

    def _match_one(ref_tt: float) -> Tuple[Optional[float], Optional[float]]:
        mask = (two_theta >= ref_tt - window_half_deg) & (two_theta <= ref_tt + window_half_deg)
        if not np.any(mask):
            return None, None
        xw = two_theta[mask]; yw = intensity_corr[mask]
        if len(yw) >= 3:
            cand = np.where((yw[1:-1] > yw[:-2]) & (yw[1:-1] > yw[2:]))[0] + 1
        else:
            cand = np.array([], dtype=int)
        j = cand[np.argmax(yw[cand])] if cand.size else int(np.argmax(yw))
        return float(xw[j]), float(yw[j])

    panels = [
        ("Quartz (101)", "(101)", QUARTZ_REF, True),
        ("Dol (104)",    "(104)", DOLOMITE_REF, False),
        ("Dol (110)",    "(110)", DOLOMITE_REF, False),
        ("Dol (015)",    "(015)", DOLOMITE_REF, False),
        ("Dol (021)",    "(021)", DOLOMITE_REF, False),
    ]

    fig, axs = plt.subplots(1, 5, figsize=(16, 3), constrained_layout=True)
    seed_leg = used_leg = height_leg = False
    ordering = {"(015)", "(021)", "(101)"}

    for ax, (title, hkl, refdict, is_quartz) in zip(axs, panels):
        tt_seed = refdict.get(hkl, None)
        ax.set_title(title)
        ax.set_xlabel("2θ ($^\circ$)"); ax.set_ylabel("I (corr)")

        if tt_seed is None:
            ax.text(0.5, 0.5, f"No ref for {hkl}", ha="center", va="center", transform=ax.transAxes)
            continue

        # Window FIRST (so xmin/xmax exist for all later calls)
        xmin, xmax = tt_seed - window_half_deg, tt_seed + window_half_deg
        mask = (two_theta >= xmin) & (two_theta <= xmax)

        # Reference seed line behind
        sc = "grey" if is_quartz else ("orange" if hkl in ordering else "gold")
        ax.axvline(tt_seed, color=sc, linestyle="--", linewidth=1.2, zorder=1,
                   label="seed (ref)" if not seed_leg else None); seed_leg = True

        # Used peak (vertical) and height (horizontal)
        tt_pk, I_pk = _match_one(tt_seed)
        if tt_pk is not None and I_pk is not None:
            ax.axvline(tt_pk, color="tab:red", linestyle="-", linewidth=1.2, zorder=3,
                       label="used peak" if not used_leg else None); used_leg = True
            ax.hlines(I_pk, xmin, xmax, color="tab:red", linestyles="dotted", linewidth=1.0, zorder=3,
                      label="peak height" if not height_leg else None); height_leg = True
            ax.text(0.02, 0.92, f"{hkl}: {tt_pk:.2f}°, I={I_pk:.0f}",
                    transform=ax.transAxes, fontsize=8, ha="left", va="top",
                    bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))

        # Spectrum on top
        ax.plot(two_theta[mask], intensity_corr[mask], lw=1, zorder=5)
        ax.set_xlim(xmin, xmax)

    # dedupe legend
    handles, labels = fig.axes[-1].get_legend_handles_labels()
    by = dict(zip(labels, handles))
    if by:
        fig.legend(by.values(), by.keys(), fontsize=8, loc="upper right")

    fig.suptitle(
        f"Peak check: {base}\n"
        f"Method: tallest local max in ±{window_half_deg:.2f}° window (fallback tallest point). "
        f"Seed colors: grey=quartz, orange=ordering dolomite, gold=other dolomite.",
        fontsize=10, y=1.2
    )
    out_png = os.path.join(outdir, f"{base}_peaks.png")
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)

# ---------- Label file parsing (supports optional hex color) ----------
def build_label_maps(label_file: Optional[str]) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Returns:
      label_map: normalized file key (basename or stem, lowercased) -> display label
      preferred_color: display label -> hex color (if provided)
    """
    def _norm_keys(s: str):
        base = os.path.basename(str(s)).strip()
        base_lower = base.casefold()
        stem_lower = os.path.splitext(base_lower)[0]
        return base_lower, stem_lower

    def _norm_hex(h: str):
        if not h:
            return None
        h = h.strip()
        if h.startswith("#"):
            h = h[1:]
        if len(h) != 6:
            return None
        try:
            int(h, 16)
        except ValueError:
            return None
        return "#" + h.lower()

    label_map: Dict[str, str] = {}
    preferred_color: Dict[str, str] = {}
    if not label_file:
        return label_map, preferred_color

    try:
        with open(label_file, "r") as fh:
            for raw in fh:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split(",", 2)
                if len(parts) < 2:
                    continue
                left = parts[0].strip()
                disp = parts[1].strip()
                hexc = _norm_hex(parts[2]) if len(parts) == 3 else None
                base_lower, stem_lower = _norm_keys(left)
                label_map[base_lower] = disp
                label_map[stem_lower]  = disp
                if hexc:
                    preferred_color[disp] = hexc
    except Exception as e:
        print(f"[WARN] Could not read --label-files: {e}")
    return label_map, preferred_color

# ---------- Core per-file processing ----------
def process_file(path: str,
                 tol_deg: float,
                 min_rel_height: float,
                 min_distance_pts: int,
                 apply_offset: Optional[float],
                 quartz_main_tt: float,
                 cal_window: float,
                 baseline_mode: str,
                 bg_window_deg: float,
                 asls_lam: float,
                 asls_p: float) -> Dict[str, Any]:

    two_theta, intensity = parse_xrdml_robust(path)

    # --- Quartz(101) after median baseline, BEFORE any 2θ offset (for calibration QC) ---
    I_medcorr_pre = apply_rolling_median_baseline(two_theta, intensity, bg_window_deg=bg_window_deg)
    q_med = find_peak_near(two_theta, I_medcorr_pre, quartz_main_tt, search_window=cal_window)
    q101_tt_med_pre = float(q_med[0]) if q_med is not None else None
    q101_I_med_pre  = float(q_med[1]) if q_med is not None else None

    # record quartz measurement (pre-shift) for reporting
    q_found = find_peak_near(two_theta, intensity, quartz_main_tt, search_window=cal_window)
    measured_quartz = float(q_found[0]) if q_found is not None else None

    # apply 2θ offset (None => no shift)
    if apply_offset is not None:
        two_theta = two_theta - float(apply_offset)

    # background correction
    if baseline_mode == "median":
        intensity = apply_rolling_median_baseline(two_theta, intensity, bg_window_deg=bg_window_deg)
    elif baseline_mode == "asls":
        intensity, _ = baseline_asls(intensity, lam=asls_lam, p=asls_p, niter=10)
    # else: none

    # peaks
    peak_idxs = simple_peak_find(two_theta, intensity, min_rel_height=min_rel_height, min_distance_pts=min_distance_pts)

    # matches
    dol = match_targets(two_theta, intensity, peak_idxs, DOLOMITE_REF, tolerance=tol_deg)
    qtz = match_targets(two_theta, intensity, peak_idxs, QUARTZ_REF, tolerance=tol_deg)
    cal = match_targets(two_theta, intensity, peak_idxs, CALCITE_REF, tolerance=tol_deg)
    cal_heights = [m.get("height") for m in cal.values() if m.get("height") is not None]
    max_calcite_I = float(max(cal_heights)) if cal_heights else 0.0
    has_calcite = bool(max_calcite_I > 100.0)  # tweak as needed

    # build row
    hkls_dolo = ["(015)", "(021)", "(101)", "(104)", "(110)", "(113)"]
    row: Dict[str, Any] = {
        "file": path,
        "n_points": int(len(two_theta)),
        "n_peaks_detected": int(len(peak_idxs)),
        "Q101_tt_med_preOffset": q101_tt_med_pre,   # 2θ of quartz(101), median-corrected, pre-offset
        "Q101_I_med_preOffset":  q101_I_med_pre,    # intensity of that peak
        "applied_offset_deg": float(apply_offset) if apply_offset is not None else None,
        "measured_quartz_101_2theta": measured_quartz,
        "has_calcite": has_calcite,
        "max_calcite_intensity": max_calcite_I,
        "baseline_mode": baseline_mode,
        "bg_window_deg": bg_window_deg,
        "asls_lam": asls_lam,
        "asls_p": asls_p,
    }
    row.update(expand_matches_as_row(dol, hkls_dolo))
    row.update(compute_dolomite_metrics(dol))
    return row

# ---------- Main ----------
def main():
    p = argparse.ArgumentParser(description="Dolomite/quartz peak extraction & ordering metrics from .xrdml")
    p.add_argument("inputs", nargs="+", help="One or more .xrdml files (or globs)")

    # detection / matching
    p.add_argument("--tol", type=float, default=0.30, help="2θ tolerance for matching ($^\circ$)")
    p.add_argument("--min-rel-height", type=float, default=0.02, help="Peak threshold as fraction of max intensity")
    p.add_argument("--min-distance-pts", type=int, default=5, help="Minimum index separation between peaks")

    # baseline
    p.add_argument("--baseline-mode", choices=["none", "median", "asls"], default="median",
                   help="Background: none | median | asls")
    p.add_argument("--bg-window-deg", type=float, default=0.6, help="Median baseline window ($^\circ$)")
    p.add_argument("--asls-lam", type=float, default=1e5, help="ASLS smoothness (1e4..1e7)")
    p.add_argument("--asls-p", type=float, default=0.01, help="ASLS asymmetry (0.001..0.1)")

    # offsets
    p.add_argument("--per-scan-quartz", action="store_true",
                   help="Use per-scan quartz(101) offset; otherwise use --manual-offset if provided.")
    p.add_argument("--fallback-offset", type=float, default=0.0942,
                   help="Fallback offset ($^\circ$) if per-scan quartz peak not found")
    p.add_argument("--manual-offset", type=float, default=None, help="Apply a fixed global offset ($^\circ$) to all scans")

    # quartz reference for calibration display & seeds in spectra
    p.add_argument("--quartz-main", type=float, default=26.6737, help="Quartz (101) 2θ reference")
    p.add_argument("--cal-window", type=float, default=0.8, help="Half-window ($^\circ$) to search around quartz (101)")

    # outputs
    p.add_argument("--single-csv-out", type=str, default="out/summary_single.csv", help="Combined CSV path")
    p.add_argument("--plots-outdir", type=str, default="out/plots", help="Directory for summary plots")
    p.add_argument("--plot-spectra", action="store_true", help="Write corrected-only spectrum per scan")
    p.add_argument("--spectra-outdir", type=str, default="out/spectra_corrected", help="Folder for spectra PNGs")
    p.add_argument("--per-sample-peaks", action="store_true", help="Write 5-panel QC figure per scan")
    p.add_argument("--peaks-outdir", type=str, default="out/peaks_check", help="Folder for QC PNGs")
    p.add_argument("--peak-window-deg", type=float, default=0.7, help="Half-window for QC panels ($^\circ$)")

    # labeling
    p.add_argument("--label-files", type=str, default=None,
                   help="Text file with lines '<filename>,<label>[,<hexcolor>]'. Comma, no spaces.")

    args = p.parse_args()

    # expand globs
    files: List[str] = []
    for patt in args.inputs:
        expanded = glob.glob(patt)
        files.extend(expanded if expanded else [patt])

    # labels & colors
    label_map, preferred_color = build_label_maps(args.label_files)
    import matplotlib.cm as cm
    label_values = sorted(set(label_map.values()))
    cmap = cm.get_cmap("tab10", max(1, len(label_values)))
    auto_color = {lab: cmap(i) for i, lab in enumerate(label_values)}
    final_color_for_label = {**auto_color, **preferred_color}

    def get_label_and_color(fname: str) -> Tuple[Optional[str], Any]:
        base = os.path.basename(str(fname)).strip().casefold()
        stem = os.path.splitext(base)[0]
        disp = label_map.get(base) or label_map.get(stem)
        if disp:
            return disp, final_color_for_label.get(disp, "black")
        return None, "black"

    # extraction loop
    rows: List[Dict[str, Any]] = []
    os.makedirs(os.path.dirname(args.single_cvs_out) if hasattr(args, "single_cvs_out") else os.path.dirname(args.single_csv_out) or ".", exist_ok=True)
    for f in files:
        try:
            offset_to_apply: Optional[float] = None
            if args.per_scan_quartz:
                tt_tmp, ii_tmp = parse_xrdml_robust(f)
                per_scan = quartz_offset_for_scan(tt_tmp, ii_tmp, quartz_main_tt=args.quartz_main, search_window=args.cal_window)
                offset_to_apply = float(per_scan) if per_scan is not None else float(args.fallback_offset)
            else:
                offset_to_apply = args.manual_offset

            # plots per scan (use SAME corrections as metrics)
            if args.plot_spectra or args.per_sample_peaks:
                t_raw, I_raw = parse_xrdml_robust(f)
                t_plot = t_raw - (offset_to_apply or 0.0)
                if args.baseline_mode == "median":
                    I_plot = apply_rolling_median_baseline(t_plot, I_raw, bg_window_deg=args.bg_window_deg)
                elif args.baseline_mode == "asls":
                    I_plot, _ = baseline_asls(I_raw, lam=args.asls_lam, p=args.asls_p, niter=10)
                else:
                    I_plot = I_raw

                if args.plot_spectra:
                    save_corrected_spectrum_plot(
                        path=f,
                        two_theta=t_plot,
                        intensity_corr=I_plot,
                        spectra_outdir=args.spectra_outdir,
                        quartz_main_tt=args.quartz_main,
                        dolo_refs=DOLOMITE_REF
                    )
                if args.per_sample_peaks:
                    save_peak_check_panels(
                        path=f,
                        two_theta=t_plot,
                        intensity_corr=I_plot,
                        window_half_deg=args.peak_window_deg,
                        dolo_refs=DOLOMITE_REF,
                        quartz_refs=QUARTZ_REF,
                        outdir=args.peaks_outdir
                    )

            row = process_file(
                f,
                tol_deg=args.tol,
                min_rel_height=args.min_rel_height,
                min_distance_pts=args.min_distance_pts,
                apply_offset=offset_to_apply,
                quartz_main_tt=args.quartz_main,
                cal_window=args.cal_window,
                baseline_mode=args.baseline_mode,
                bg_window_deg=args.bg_window_deg,
                asls_lam=args.asls_lam,
                asls_p=args.asls_p
            )
            rows.append(row)
            print(f"[OK] {os.path.basename(f)}  peaks={row['n_peaks_detected']}")
        except Exception as e:
            print(f"[ERR] {f}: {e}")

    # write single CSV
    if rows:
        df = pd.DataFrame(rows)
        os.makedirs(os.path.dirname(args.single_csv_out) or ".", exist_ok=True)
        df.to_csv(args.single_csv_out, index=False)
        print(f"[WROTE] Single CSV: {args.single_csv_out}")

        # -------- Summary plots --------
        os.makedirs(args.plots_outdir, exist_ok=True)
        
        # Compute/repair the metrics we want
        def _num(col):
            return pd.to_numeric(df.get(col), errors="coerce")
        
        df["I015_over_I110"] = (_num("I_015") / _num("I_110")).replace([np.inf, -np.inf], np.nan)
        df["I021_over_I110"] = (_num("I_021") / _num("I_110")).replace([np.inf, -np.inf], np.nan)
        df["HWHM_015"]       = 0.5 * _num("FWHM_015")
        df["delta_tt_104"]   = _num("tt_104") - DOLOMITE_REF["(104)"]
        
        # 1) Grey histogram: I(015)/I(110)
        plt.figure()
        vals = df["I015_over_I110"].dropna().astype(float)
        plt.hist(vals, bins=30, color="0.7", edgecolor="k")  # grey bars
        plt.xlabel("I(015)/I(110)"); plt.ylabel("Count")
        
        # overlay label strip (use SCATTER so legend has circles only, no lines)
        if label_map:
            y0, y1 = plt.ylim()
            y_strip = y0 + 0.02 * (y1 - y0)
            for _, r in df.dropna(subset=["I015_over_I110"]).iterrows():
                disp, color = get_label_and_color(r["file"])
                if not disp:
                    continue
                x = float(r["I015_over_I110"])
                plt.scatter([x], [y_strip], s=60, facecolor=color, edgecolor="k", label=disp)
        
            # de-dupe legend
            handles, labels = plt.gca().get_legend_handles_labels()
            by = dict(zip(labels, handles))
            if by:
                plt.legend(by.values(), by.keys(), fontsize=8, loc="best")
        
        plt.tight_layout()
        plt.savefig(os.path.join(args.plots_outdir, "hist_I015_over_I110.png"), dpi=200)
        plt.close()
        
        # 2) Scatter: I(015)/I(110) vs HWHM(015)
        plt.figure()
        sub = df.dropna(subset=["I015_over_I110", "HWHM_015"]).copy()
        plt.scatter(sub["I015_over_I110"].astype(float), sub["HWHM_015"].astype(float),
                    s=16, alpha=0.7, color="0.6")
        if label_map:
            for _, r in sub.iterrows():
                disp, color = get_label_and_color(r["file"])
                if not disp: 
                    continue
                x = float(r["I015_over_I110"]); y = float(r["HWHM_015"])
                plt.scatter(x, y, s=60, edgecolor="k", facecolor=color, label=disp)
                plt.annotate(disp, (x, y), xytext=(5,5), textcoords="offset points", fontsize=8, color=color)
            h, l = plt.gca().get_legend_handles_labels()
            by = dict(zip(l, h))
            if by: plt.legend(by.values(), by.keys(), fontsize=8, loc="best")
        plt.xlabel("I(015)/I(110)"); plt.ylabel("HWHM(015) ($^\circ$)")
        plt.tight_layout()
        plt.savefig(os.path.join(args.plots_outdir, "scatter_I015over110_vs_HWHM015.png"), dpi=200)
        plt.close()
        
        # 3) Box plots with colored overlays (I015/110 and I021/110)
        plt.figure()
        vals_015_110 = df["I015_over_I110"].dropna().astype(float)
        vals_021_110 = df["I021_over_I110"].dropna().astype(float)
        plt.boxplot([vals_015_110, vals_021_110], labels=["I(015)/I(110)", "I(021)/I(110)"])
        if label_map:
            def _jit(xc): return xc + (np.random.rand() - 0.5) * 0.06
            for _, r in df.iterrows():
                disp, color = get_label_and_color(r["file"])
                if not disp: continue
                if pd.notna(r.get("I015_over_I110")):
                    plt.scatter(_jit(1.0), float(r["I015_over_I110"]), c=[color], edgecolor="k", label=disp)
                if pd.notna(r.get("I021_over_I110")):
                    plt.scatter(_jit(2.0), float(r["I021_over_I110"]), c=[color], edgecolor="k", label=disp)
            h, l = plt.gca().get_legend_handles_labels()
            by = dict(zip(l, h))
            if by: plt.legend(by.values(), by.keys(), fontsize=8, loc="best")
        plt.ylabel("Intensity ratio")
        plt.tight_layout()
        plt.savefig(os.path.join(args.plots_outdir, "box_I_ratios_over_110.png"), dpi=200)
        plt.close()
        
        # 4) Grey histogram of applied offsets (unchanged)
        plt.figure()
        vals_off = df["applied_offset_deg"].dropna().astype(float)
        plt.hist(vals_off, bins=30, color="0.7", edgecolor="k")
        plt.xlabel("Applied 2θ offset ($^\circ$)"); plt.ylabel("Count")
        # (optional) label strip
        if label_map:
            y0, y1 = plt.ylim()
            y_strip = y0 + 0.02 * (y1 - y0)
            for _, r in df.dropna(subset=["applied_offset_deg"]).iterrows():
                disp, color = get_label_and_color(r["file"])
                if not disp:
                    continue
                x = float(r["applied_offset_deg"])
                plt.scatter([x], [y_strip], s=60, facecolor=color, edgecolor="k", label=disp)
        
            handles, labels = plt.gca().get_legend_handles_labels()
            by = dict(zip(labels, handles))
            if by:
                plt.legend(by.values(), by.keys(), fontsize=8, loc="best")
        plt.tight_layout()
        plt.savefig(os.path.join(args.plots_outdir, "hist_applied_offsets.png"), dpi=200)
        plt.close()

        
        # 5) Scatter: I(015)/I(110) vs Δ2θ(104)  (tt_104 - seed_104)
        plt.figure()
        sub2 = df.dropna(subset=["I015_over_I110", "delta_tt_104"]).copy()
        plt.scatter(sub2["I015_over_I110"].astype(float), sub2["delta_tt_104"].astype(float),
                    s=16, alpha=0.7, color="0.6")
        if label_map:
            for _, r in sub2.iterrows():
                disp, color = get_label_and_color(r["file"])
                if not disp: continue
                x = float(r["I015_over_I110"]); y = float(r["delta_tt_104"])
                plt.scatter(x, y, s=60, edgecolor="k", facecolor=color, label=disp)
                plt.annotate(disp, (x, y), xytext=(5,5), textcoords="offset points", fontsize=8, color=color)
            h, l = plt.gca().get_legend_handles_labels()
            by = dict(zip(l, h))
            if by: plt.legend(by.values(), by.keys(), fontsize=8, loc="best")
        plt.axhline(0.0, ls="--", lw=1, color="0.3")
        plt.xlabel("I(015)/I(110)"); plt.ylabel("Δ2θ(104) ($^\circ$)")
        plt.tight_layout()
        plt.savefig(os.path.join(args.plots_outdir, "scatter_I015over110_vs_delta104.png"), dpi=200)
        plt.close()

        # Quartz(101) histogram (median-corrected, pre-offset), only where I > 50
        plt.figure()
        mask = (df["Q101_I_med_preOffset"].fillna(0) > 50) & df["Q101_tt_med_preOffset"].notna()
        vals_q = df.loc[mask, "Q101_tt_med_preOffset"].astype(float)
        
        plt.hist(vals_q, bins=30, color="0.7", edgecolor="k")
        plt.xlabel("Quartz (101) 2θ ($^\circ$)  [median baseline, pre-offset]")
        plt.ylabel("Count")
        
        # optional label strip
        if label_map and not vals_q.empty:
            y0, y1 = plt.ylim()
            y_strip = y0 + 0.02 * (y1 - y0)
            for _, r in df.loc[mask].iterrows():
                disp, color = get_label_and_color(r["file"])
                if not disp:
                    continue
                x = float(r["Q101_tt_med_preOffset"])
                plt.scatter([x], [y_strip], s=60, facecolor=color, edgecolor="k", label=disp)
            handles, labels = plt.gca().get_legend_handles_labels()
            by = dict(zip(labels, handles))
            if by:
                plt.legend(by.values(), by.keys(), fontsize=8, loc="best")
        
        plt.tight_layout()
        plt.savefig(os.path.join(args.plots_outdir, "hist_quartz101_tt_median_preOffset.png"), dpi=200)
        plt.close()

        # ---- Quartz (101) pre-offset stats (median baseline), I > 50 ----
        try:
            seed = float(args.quartz_main)  # your quartz(101) seed, e.g. 26.6737
            col_tt = "Q101_tt_med_preOffset"
            col_I  = "Q101_I_med_preOffset"
        
            if col_tt in df.columns and col_I in df.columns:
                mask = (df[col_I].fillna(0) > 50) & df[col_tt].notna()
                vals = df.loc[mask, col_tt].astype(float).to_numpy()
                n = int(vals.size)
        
                print("\n=== Quartz (101) pre-offset (median) stats, I > 50 ===")
                if n > 0:
                    mean_tt = float(vals.mean())
                    std_tt  = float(vals.std(ddof=1))  # sample std
                    mean_off = mean_tt - seed
                    print(f"N = {n}    mean = {mean_tt:.4f}°    std = {std_tt:.4f}°    "
                          f"meanΔ = {mean_off:+.4f}°  (vs seed {seed:.4f}°)")
                else:
                    print("No qualifying peaks.")
            else:
                print("\n[WARN] Q101_*_med_preOffset columns not found; stats not computed.")
        except Exception as e:
            print(f"\n[WARN] Failed to compute quartz stats: {e}")

        plt.figure()
        sub3 = df.dropna(subset=["I015_over_I110", "FWHM_104"]).copy()
        plt.scatter(sub3["I015_over_I110"].astype(float), sub3["FWHM_104"].astype(float),
                    s=16, alpha=0.7, color="0.6")
        
        if label_map:
            for _, r in sub3.iterrows():
                disp, color = get_label_and_color(r["file"])
                if not disp:
                    continue
                x = float(r["I015_over_I110"]); y = float(r["FWHM_104"])
                plt.scatter(x, y, s=60, edgecolor="k", facecolor=color, label=disp)
                plt.annotate(disp, (x, y), xytext=(5,5), textcoords="offset points",
                             fontsize=8, color=color)
            h, l = plt.gca().get_legend_handles_labels()
            by = dict(zip(l, h))
            if by:
                plt.legend(by.values(), by.keys(), fontsize=8, loc="best")
        
        plt.xlabel("I(015)/I(110)")
        plt.ylabel("FWHM(104) ($^\circ$)")
        plt.tight_layout()
        plt.savefig(os.path.join(args.plots_outdir, "scatter_I015over110_vs_FWHM104.png"), dpi=200)
        plt.close()

        # --- New scatter: Δ2θ(104) vs FWHM(104)
        plt.figure()
        sub5 = df.dropna(subset=["delta_tt_104", "FWHM_104"]).copy()
        plt.scatter(sub5["delta_tt_104"].astype(float), sub5["FWHM_104"].astype(float),
                    s=16, alpha=0.7, color="0.6")
        
        if label_map:
            for _, r in sub5.iterrows():
                disp, color = get_label_and_color(r["file"])
                if not disp:
                    continue
                x = float(r["delta_tt_104"]); y = float(r["FWHM_104"])
                plt.scatter(x, y, s=60, edgecolor="k", facecolor=color, label=disp)
                plt.annotate(disp, (x, y), xytext=(5,5), textcoords="offset points",
                             fontsize=8, color=color)
            h, l = plt.gca().get_legend_handles_labels()
            by = dict(zip(l, h))
            if by:
                plt.legend(by.values(), by.keys(), fontsize=8, loc="best")
        
        plt.axvline(0.0, ls="--", lw=1, color="0.3")
        plt.xlabel("Δ2θ(104) ($^\circ$)")
        plt.ylabel("FWHM(104) ($^\circ$)")
        plt.tight_layout()
        plt.savefig(os.path.join(args.plots_outdir, "scatter_delta104_vs_FWHM104.png"), dpi=200)
        plt.close()


if __name__ == "__main__":
    main()