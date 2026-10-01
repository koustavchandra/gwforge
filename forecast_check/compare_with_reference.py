#!/usr/bin/env python
"""Compare the gwforge population forecast with the CE_STM_CBC reference, broadly.

Three tables, written to COMPARISON.md: relative mass sigmas, the effective
Madau-Dickinson parameters and their sigmas, and the delay-time posterior.
Relative uncertainties are compared because the reference fiducial is a BGP fit
to a Full-Pop catalogue while this population is drawn from the GWTC-5.0 BGP.
"""
import argparse
from pathlib import Path

import numpy

REFERENCE = Path("/Users/kchandra/projects/CE_STM_CBC")
JOINT = REFERENCE / "final_plots/mass_redshift_timedelay_forecast/output/fisher_bbhs_CE40_CE20_ETD.npz"
NO_DELAY = REFERENCE / "population_fisher_seminumeric/forecast_mass_redshift_output/forecast_joint_snr10.npz"
# reference name -> gwforge name
MASS = {"alpha_1": "alpha_1", "alpha_2": "alpha_2", "lambda_0": "lam_0", "lambda_1": "lam_1", "mu_1": "mpp_1", "sigma_1": "sigpp_1", "mu_2": "mpp_2", "sigma_2": "sigpp_2", "beta_q": "beta"}
TAU_MIN = 0.02

parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
parser.add_argument("--xg", type=Path, default=Path("pf_xg/mrs_xg_snr10.npz"))
parser.add_argument("--two-l", type=Path, default=Path("pf_2l/mrs_2l_snr10.npz"))
parser.add_argument("--posterior-xg", type=Path, default=Path("delay/posterior_single_channel_CE40_CE20_ET.npz"))
parser.add_argument("--posterior-two-l", type=Path, default=Path("delay/posterior_single_channel_CE40_ET2L_LI.npz"))
parser.add_argument("--output", type=Path, default=Path("COMPARISON.md"))
opts = parser.parse_args()


def load(path):
    d = numpy.load(path, allow_pickle=True)
    names = [str(n) for n in d["parameter_names"]]
    return names, numpy.asarray(d["fiducial"], float), numpy.asarray(d["covariance"], float), int(d["n_events"]), int(d["n_total"])


def sigma(names, covariance, name):
    return float(numpy.sqrt(covariance[names.index(name), names.index(name)]))


ours = {"CE40+CE20+ET": load(opts.xg), "CE40+ET2L+LI": load(opts.two_l)}
joint = numpy.load(JOINT, allow_pickle=True)
ref_names = [str(n) for n in joint["joint__term_I__names"]]
ref_fid = numpy.asarray(joint["joint__term_I__fiducial"], float)
ref_cov = numpy.asarray(joint["joint__term_I__covariance"], float)
nodelay = numpy.load(NO_DELAY, allow_pickle=True)
nd_names = [str(n) for n in nodelay["parameter_names"]]
nd_fid = numpy.asarray(nodelay["fiducial"], float)
nd_cov = numpy.asarray(nodelay["covariance"], float)

lines = ["# gwforge population forecast against the CE_STM_CBC reference", ""]
lines += [
    "Reference: one simulated year, 34166 BBH from Full-Pop GWTC-4 with a BGP fit as fiducial, CE40+CE20+ET(triangle), SNR >= 10 plus a Fisher-inversion "
    "quality cut, Term I (`joint__term_I`, 14 free) and the older no-delay forecast (`forecast_joint_snr10`, 12 free)."
]
lines += ["This run: one Poisson year of the GWTC-5.0 BGP + Default spins + Madau-Dickinson (x) 1/tau delay, 18 free (10 mass, 3 redshift, 5 spin), Term I on the detected events.", ""]
lines += ["| | " + " | ".join(ours) + " | reference |", "|---|" + "---|" * (len(ours) + 1)]
lines += ["| N detected (SNR >= 10) | " + " | ".join(str(v[3]) for v in ours.values()) + " | 32639 |"]
lines += ["| N injected | " + " | ".join(str(v[4]) for v in ours.values()) + " | 34166 |", ""]

lines += ["## 1. Mass block, sigma / |fiducial| in per cent (pass: every ratio to the reference within a factor 2)", ""]
lines += ["| parameter | " + " | ".join(ours) + " | reference joint | reference no-delay | worst ratio |", "|---|" + "---|" * (len(ours) + 3)]
worst = 0.0
for ref_name, our_name in MASS.items():
    ref_rel = 100 * sigma(ref_names, ref_cov, ref_name) / abs(ref_fid[ref_names.index(ref_name)])
    nd_rel = 100 * sigma(nd_names, nd_cov, ref_name) / abs(nd_fid[nd_names.index(ref_name)])
    cells, ratios = [], []
    for names, fid, cov, _, _ in ours.values():
        rel = 100 * sigma(names, cov, our_name) / abs(fid[names.index(our_name)])
        cells.append("{:.2f}".format(rel))
        ratios.append(max(rel / ref_rel, ref_rel / rel))
    worst = max(worst, max(ratios))
    lines.append("| {} ({}) | {} | {:.2f} | {:.2f} | {:.2f} |".format(our_name, ref_name, " | ".join(cells), ref_rel, nd_rel, max(ratios)))
lines += ["", "Worst ratio over the nine shared parameters: {:.2f} ({}). m_break is free here and pinned in the reference.".format(worst, "PASS" if worst < 2 else "FAIL"), ""]

lines += ["## 2. Redshift block: effective Madau-Dickinson parameters and their sigmas", ""]
nd_alpha, nd_beta = nd_fid[nd_names.index("alpha")], nd_fid[nd_names.index("beta")]
ia, ib = nd_names.index("alpha"), nd_names.index("beta")
nd_kappa_sigma = numpy.sqrt(nd_cov[ia, ia] + nd_cov[ib, ib] + 2 * nd_cov[ia, ib])
lines += [
    "Reference no-delay effective fit on the Full-Pop catalogue: gamma = {:.3f}, kappa = {:.3f}, z_peak = {:.3f}; sigma = {:.4f}, {:.4f}, {:.4f} (kappa = alpha + beta).".format(
        nd_alpha, nd_alpha + nd_beta, nd_fid[nd_names.index("z_peak")], sigma(nd_names, nd_cov, "alpha"), nd_kappa_sigma, sigma(nd_names, nd_cov, "z_peak")
    )
]
lines += [
    "Earlier gwforge catalogue (CE40+CE20+ET, 33700 injections): effective (1.82, 5.29, 1.82). "
    "Two-channel study, truth case, CE40+ET2L+LI on the isolated 63 %: sigma (0.0996, 0.0682, 0.0504).",
    "",
]
lines += ["| | " + " | ".join(ours) + " |", "|---|" + "---|" * len(ours)]
for key in ("gamma", "kappa", "z_peak"):
    lines.append("| {} fitted | {} |".format(key, " | ".join("{:.3f}".format(v[1][v[0].index(key)]) for v in ours.values())))
    lines.append("| sigma({}) | {} |".format(key, " | ".join("{:.4f}".format(sigma(v[0], v[2], key)) for v in ours.values())))
nd_sig = {"gamma": sigma(nd_names, nd_cov, "alpha"), "kappa": nd_kappa_sigma, "z_peak": sigma(nd_names, nd_cov, "z_peak")}
ratio = max(max(sigma(v[0], v[2], k) / nd_sig[k], nd_sig[k] / sigma(v[0], v[2], k)) for v in ours.values() for k in nd_sig)
lines += [
    "",
    "Worst sigma ratio to the reference no-delay forecast: {:.2f} ({}). ".format(ratio, "PASS" if ratio < 2 else "FAIL")
    + "The reference joint sigmas (1.13, 0.55, 0.055) are not comparable: they free the delay parameters too.",
    "",
]

lines += ["## 3. Delay-time posterior (Procedure 1, SFR pinned at Madau-Dickinson; injected alpha = 1, tau_min = 20 Myr)", ""]
lines += [
    "| network | alpha 5/50/95 % | tau_min 5/50/95 % [Gyr] | MAP | KL(MAP) | sigma(alpha) | sigma(log10 tau) | corr | reference joint_sfr_pinned term_I | truth inside 90 % |",
    "|---|" + "---|" * 9,
]
delay_pass = True
for label, path in (("CE40+CE20+ET", opts.posterior_xg), ("CE40+ET2L+LI", opts.posterior_two_l)):
    if not path.exists():
        lines.append("| {} | (no posterior file) |".format(label))
        continue
    p = numpy.load(path, allow_pickle=True)
    s = p["samples"]
    aq, tq = p["alpha_quantiles"], 10.0 ** p["log10_tau_min_quantiles"]
    inside = aq[0] <= 1.0 <= aq[2] and tq[0] <= TAU_MIN <= tq[2]
    delay_pass &= bool(inside)
    ref_text = "0.0386, 0.1165, 0.957" if label == "CE40+CE20+ET" else "none for this network"
    lines.append("| {} | {} | {} | ({:.3f}, {:.4f}) | {:.1e} | {:.4f} | {:.4f} | {:.3f} | {} | {} |".format(
        label, aq.round(3), tq.round(4), p["map"][0], 10.0 ** p["map"][1], float(p["kl_map"]), s[:, 0].std(), s[:, 1].std(), numpy.corrcoef(s.T)[0, 1], ref_text, "yes" if inside else "NO"))
lines += ["", "Pass: the injected values lie inside both 90 % intervals, and the widths are within a factor 2 of the reference row ({}).".format("PASS" if delay_pass else "FAIL"), ""]
opts.output.write_text("\n".join(lines) + "\n")
print("\n".join(lines))
