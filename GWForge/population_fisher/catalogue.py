r"""Assemble a detected catalogue from GWForge population and SNR files.

The pipeline this consumes is the ordinary GWForge one:

.. code-block:: bash

   gwforge_population   --config-file bgp.ini --output-file population.h5
   gwforge_optimal_snr  --injection-file population.h5 --output-file snr.h5 \
                        --ifos CE40 CE20 ET

``gwforge_population`` writes the source parameters, ``gwforge_optimal_snr``
writes one ``{IFO}_optimal_snr`` column per detector plus the ``index`` of the
injection it came from. This module joins the two on that index, forms the
network SNR in quadrature, applies a threshold, and hands back the events dict
the population models expect.

Joining on ``index`` rather than on row order matters as soon as the SNR stage
is split across jobs with ``--start-index``/``--end-index``, which is the normal
way to run 30,000 sources. Older SNR files predate the column; those are joined
positionally and the overlap is checked against the injection file, which will
notice a misalignment rather than silently forecasting for the wrong sources.

What ``n_total`` is for
-----------------------

The Term-I Fisher is a sum over *detected* events, so its scale is
:math:`N_{\rm det}` and nothing else -- the total number of injections never
enters the matrix. It is carried anyway, because the detected fraction is the
first thing to look at when a forecast comes out surprising: a threshold that
admits 99% of a population is telling you something different from one that
admits 3%.
"""

import logging

import h5py
import numpy

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)

# Source parameters the population models can ask for, pulled from the
# injection file when present.
CATALOGUE_KEYS = (
    "mass_1_source",
    "mass_2_source",
    "mass_1",
    "mass_2",
    "redshift",
    "luminosity_distance",
    "a_1",
    "a_2",
    "tilt_1",
    "tilt_2",
    "chi_1",
    "chi_2",
)

# Columns used to confirm that a positional join lined up. Compared with
# ``rtol`` below.
_ALIGNMENT_KEYS = ("mass_1", "luminosity_distance")
_ALIGNMENT_TOLERANCE = 1e-6

# ``gwforge_optimal_snr`` expands a triangular detector into three arms.
TRIANGULAR_ARMS = {"ET": ("ET1", "ET2", "ET3")}


class Catalogue:
    """A detected catalogue, and the counts that describe how it was cut.

    Attributes
    ----------
    events : dict
        Source parameters of the detected events, plus derived ``mass_ratio``
        and ``cos_tilt_1``/``cos_tilt_2`` where the inputs allow.
    network_snr : numpy.ndarray
        Network optimal SNR of the detected events.
    detectors : list of str
        Detector columns that went into the network SNR.
    snr_threshold : float
    n_total : int
        Injections in the population file.
    n_with_snr : int
        Injections the SNR stage actually covered.
    n_detected : int
        Injections above threshold. GWForge applies no duty cycle, so this is
        also the analysis sample.
    """

    def __init__(
        self,
        events,
        network_snr,
        detectors,
        snr_threshold,
        n_total,
        n_with_snr,
    ):
        self.events = events
        self.network_snr = network_snr
        self.detectors = list(detectors)
        self.snr_threshold = float(snr_threshold)
        self.n_total = int(n_total)
        self.n_with_snr = int(n_with_snr)
        self.n_detected = int(len(network_snr))

    @property
    def detected_fraction(self):
        """Fraction of the SNR-covered injections that passed the threshold."""
        return self.n_detected / self.n_with_snr if self.n_with_snr else float("nan")

    def summary(self):
        """A plain-text description of the cuts.

        Returns
        -------
        str
        """
        return "\n".join(
            [
                "  Detectors        : {}".format(", ".join(self.detectors)),
                "  SNR threshold    : {:g}".format(self.snr_threshold),
                "  Injections       : {}".format(self.n_total),
                "  With SNR         : {}".format(self.n_with_snr),
                "  Detected         : {} ({:.3f} of those with SNR)".format(
                    self.n_detected, self.detected_fraction
                ),
                "  Network SNR      : median {:.1f}, max {:.1f}".format(
                    numpy.median(self.network_snr) if self.n_detected else float("nan"),
                    self.network_snr.max() if self.n_detected else float("nan"),
                ),
            ]
        )

    def __repr__(self):
        return "Catalogue({} detected of {} injections)".format(
            self.n_detected, self.n_total
        )


def expand_detectors(detectors):
    """Expand triangular detectors into the arm names the SNR file uses.

    ``["CE40", "ET"]`` becomes ``["CE40", "ET1", "ET2", "ET3"]``, matching what
    ``gwforge_optimal_snr`` writes.

    Parameters
    ----------
    detectors : sequence of str

    Returns
    -------
    list of str
    """
    expanded = []
    for name in detectors:
        expanded.extend(TRIANGULAR_ARMS.get(name.upper(), (name,)))
    return expanded


def _read(path, keys=None):
    """Read one-dimensional datasets from an HDF5 file into a dict."""
    data = {}
    with h5py.File(path, "r") as handle:
        for key in handle.keys():
            if keys is not None and key not in keys:
                continue
            dataset = handle[key]
            if isinstance(dataset, h5py.Dataset) and dataset.ndim == 1:
                data[key] = dataset[:]
    return data


def network_snr_from(snr_data, detectors=None):
    """Quadrature sum of the per-detector optimal SNRs.

    Parameters
    ----------
    snr_data : dict
        As written by ``gwforge_optimal_snr``.
    detectors : sequence of str or None
        Detector names (triangles already expanded). Defaults to every
        ``*_optimal_snr`` column in the file.

    Returns
    -------
    tuple
        ``(network_snr, detectors_used)``.

    Raises
    ------
    KeyError
        If a requested detector has no column.
    """
    if detectors is None:
        detectors = sorted(
            key[: -len("_optimal_snr")]
            for key in snr_data
            if key.endswith("_optimal_snr")
        )
    missing = [
        name for name in detectors if "{}_optimal_snr".format(name) not in snr_data
    ]
    if missing:
        raise KeyError(
            "No SNR column for {}. The file has {}.".format(
                missing,
                sorted(
                    key[: -len("_optimal_snr")]
                    for key in snr_data
                    if key.endswith("_optimal_snr")
                ),
            )
        )
    total = sum(
        numpy.asarray(snr_data["{}_optimal_snr".format(name)], dtype=float) ** 2
        for name in detectors
    )
    return numpy.sqrt(total), list(detectors)


def derive_event_keys(events):
    """Add the parameters the population models want but the file does not store.

    ``mass_ratio`` from the component masses, and ``cos_tilt_i`` from the tilts.
    Both are cheap and both are what a model would otherwise recompute per call.

    Parameters
    ----------
    events : dict
        Modified in place and returned.

    Returns
    -------
    dict
    """
    if (
        "mass_ratio" not in events
        and "mass_1_source" in events
        and "mass_2_source" in events
    ):
        events["mass_ratio"] = events["mass_2_source"] / events["mass_1_source"]
    for index in (1, 2):
        tilt = "tilt_{}".format(index)
        cosine = "cos_tilt_{}".format(index)
        if cosine not in events and tilt in events:
            events[cosine] = numpy.cos(events[tilt])
    return events


def load_catalogue(
    population_file,
    snr_file,
    snr_threshold=10.0,
    detectors=None,
    keys=CATALOGUE_KEYS,
):
    """Join a population file to its SNR file and cut on network SNR.

    Parameters
    ----------
    population_file : str
        HDF5 from ``gwforge_population``.
    snr_file : str
        HDF5 from ``gwforge_optimal_snr``.
    snr_threshold : float
        Network optimal-SNR threshold for detection.
    detectors : sequence of str or None
        Detectors to sum. ``"ET"`` is expanded to its three arms. Defaults to
        every detector in the SNR file.
    keys : sequence of str
        Source parameters to carry over.

    Returns
    -------
    Catalogue

    Raises
    ------
    ValueError
        If the two files cannot be aligned, or nothing passes the threshold.
    """
    population = _read(population_file)
    snr_data = _read(snr_file)
    n_total = len(next(iter(population.values())))

    if "index" in snr_data:
        rows = numpy.asarray(snr_data["index"], dtype=int)
    else:
        logging.warning(
            "%s has no 'index' column (written by older gwforge_optimal_snr); "
            "joining positionally and verifying the overlap.",
            snr_file,
        )
        rows = numpy.arange(len(next(iter(snr_data.values()))))
    if rows.max() >= n_total:
        raise ValueError(
            "The SNR file refers to injection {} but {} has only {}. The two "
            "files are not from the same population.".format(
                rows.max(), population_file, n_total
            )
        )

    for name in _ALIGNMENT_KEYS:
        if name in snr_data and name in population:
            if not numpy.allclose(
                snr_data[name], population[name][rows], rtol=_ALIGNMENT_TOLERANCE
            ):
                raise ValueError(
                    "'{}' disagrees between {} and {} after the join. The SNR "
                    "file does not belong to this population.".format(
                        name, snr_file, population_file
                    )
                )

    if detectors is not None:
        detectors = expand_detectors(detectors)
    network_snr, detectors_used = network_snr_from(snr_data, detectors)
    detected = network_snr >= snr_threshold
    if not numpy.any(detected):
        raise ValueError(
            "No injection reaches a network SNR of {:g}; the loudest is "
            "{:.2f}.".format(snr_threshold, network_snr.max())
        )

    events = {
        name: numpy.asarray(population[name], dtype=float)[rows][detected]
        for name in keys
        if name in population
    }
    events = derive_event_keys(events)

    catalogue = Catalogue(
        events=events,
        network_snr=network_snr[detected],
        detectors=detectors_used,
        snr_threshold=snr_threshold,
        n_total=n_total,
        n_with_snr=len(network_snr),
    )
    logging.info("Loaded catalogue\n%s", catalogue.summary())
    return catalogue


def load_injections(population_file, keys=CATALOGUE_KEYS):
    """Every injection, detected or not, as an events dict.

    This is what the maximum-likelihood fits in
    :mod:`GWForge.population_fisher.fit` run on: the fiducial hyper-parameters
    describe the *astrophysical* population, so they are fitted to what was
    drawn, not to what was detected.

    Parameters
    ----------
    population_file : str
    keys : sequence of str

    Returns
    -------
    dict
    """
    population = _read(population_file)
    events = {
        name: numpy.asarray(population[name], dtype=float)
        for name in keys
        if name in population
    }
    return derive_event_keys(events)
