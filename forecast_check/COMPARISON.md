# gwforge population forecast against the CE_STM_CBC reference

Reference: one simulated year, 34166 BBH from Full-Pop GWTC-4 with a BGP fit as fiducial, CE40+CE20+ET(triangle), SNR >= 10 plus a Fisher-inversion quality cut, Term I (`joint__term_I`, 14 free) and the older no-delay forecast (`forecast_joint_snr10`, 12 free).
This run: one Poisson year of the GWTC-5.0 BGP + Default spins + Madau-Dickinson (x) 1/tau delay, 18 free (10 mass, 3 redshift, 5 spin), Term I on the detected events.

| | CE40+CE20+ET | CE40+ET2L+LI | reference |
|---|---|---|---|
| N detected (SNR >= 10) | 33255 | 33113 | 32639 |
| N injected | 33519 | 33519 | 34166 |

## 1. Mass block, sigma / |fiducial| in per cent (pass: every ratio to the reference within a factor 2)

| parameter | CE40+CE20+ET | CE40+ET2L+LI | reference joint | reference no-delay | worst ratio |
|---|---|---|---|---|---|
| alpha_1 (alpha_1) | 1.96 | 1.97 | 2.12 | 2.11 | 1.08 |
| alpha_2 (alpha_2) | 3.80 | 3.80 | 1.29 | 1.28 | 2.95 |
| lam_0 (lambda_0) | 1.13 | 1.13 | 1.03 | 1.03 | 1.10 |
| lam_1 (lambda_1) | 0.71 | 0.71 | 0.66 | 0.66 | 1.08 |
| mpp_1 (mu_1) | 0.06 | 0.06 | 0.09 | 0.09 | 1.53 |
| sigpp_1 (sigma_1) | 0.73 | 0.73 | 0.74 | 0.74 | 1.02 |
| mpp_2 (mu_2) | 0.81 | 0.81 | 0.48 | 0.48 | 1.68 |
| sigpp_2 (sigma_2) | 5.01 | 5.01 | 3.72 | 3.71 | 1.35 |
| beta (beta_q) | 3.05 | 3.06 | 1.42 | 1.41 | 2.16 |

Worst ratio over the nine shared parameters: 2.95 (FAIL). m_break is free here and pinned in the reference.

## 2. Redshift block: effective Madau-Dickinson parameters and their sigmas

Reference no-delay effective fit on the Full-Pop catalogue: gamma = 1.881, kappa = 5.364, z_peak = 1.818; sigma = 0.0745, 0.0523, 0.0381 (kappa = alpha + beta).
Earlier gwforge catalogue (CE40+CE20+ET, 33700 injections): effective (1.82, 5.29, 1.82). Two-channel study, truth case, CE40+ET2L+LI on the isolated 63 %: sigma (0.0996, 0.0682, 0.0504).

| | CE40+CE20+ET | CE40+ET2L+LI |
|---|---|---|
| gamma fitted | 1.819 | 1.819 |
| sigma(gamma) | 0.0744 | 0.0745 |
| kappa fitted | 5.292 | 5.292 |
| sigma(kappa) | 0.0519 | 0.0520 |
| z_peak fitted | 1.824 | 1.824 |
| sigma(z_peak) | 0.0384 | 0.0385 |

Worst sigma ratio to the reference no-delay forecast: 1.01 (PASS). The reference joint sigmas (1.13, 0.55, 0.055) are not comparable: they free the delay parameters too.

## 3. Delay-time posterior (Procedure 1, SFR pinned at Madau-Dickinson; injected alpha = 1, tau_min = 20 Myr)

| network | alpha 5/50/95 % | tau_min 5/50/95 % [Gyr] | MAP | KL(MAP) | sigma(alpha) | sigma(log10 tau) | corr | reference joint_sfr_pinned term_I | truth inside 90 % |
|---|---|---|---|---|---|---|---|---|---|
| CE40+CE20+ET | [0.926 0.986 1.05 ] | [0.0091 0.0176 0.0281] | (0.983, 0.0175) | 1.6e-06 | 0.0380 | 0.1511 | 0.942 | 0.0386, 0.1165, 0.957 | yes |
| CE40+ET2L+LI | [0.925 0.986 1.05 ] | [0.0091 0.0176 0.0281] | (0.983, 0.0175) | 1.6e-06 | 0.0382 | 0.1521 | 0.941 | none for this network | yes |

Pass: the injected values lie inside both 90 % intervals, and the widths are within a factor 2 of the reference row (PASS).

