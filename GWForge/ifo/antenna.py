import logging

import numpy
from bilby.core.utils import ra_dec_to_theta_phi
from bilby_cython.geometry import greenwich_mean_sidereal_time
from lal import C_SI, MTSUN_SI

RESPONSE_PARAMETERS = (
    "ra",
    "dec",
    "psi",
    "geocent_time",
    "chirp_mass",
    "symmetric_mass_ratio",
    "mass_1",
    "mass_2",
    "mass_ratio",
)


def finite_size_factor(x, y):
    r"""Single-arm transfer function for finite light travel time.

    This is the frequency-domain response of one Fabry-Perot arm to a wave
    arriving from direction :math:`\hat\Omega`, from Rakhmanov, Romano & Whelan
    (`arXiv:0808.3805 <https://arxiv.org/abs/0808.3805>`_):

    .. math::

       D(x, y) = \tfrac{1}{2}\left[
            e^{-i\pi x(1+y)}\,\mathrm{sinc}\big(x(1-y)\big)
          + e^{ i\pi x(1-y)}\,\mathrm{sinc}\big(x(1+y)\big)\right]

    Parameters
    ----------
    x : array_like
        :math:`fL/c`, the arm length in wavelengths.
    y : array_like
        :math:`-\hat\Omega\cdot\hat{a}`, the cosine of the angle between the arm
        and the direction of propagation.

    Returns
    -------
    numpy.ndarray
        Complex transfer function. ``D(0, y) == 1`` for every ``y``, which is the
        long-wavelength limit -- so setting this to 1 recovers bilby exactly.

    Notes
    -----
    ``numpy.sinc(z)`` is :math:`\sin(\pi z)/(\pi z)`, which is what the formula
    above already assumes; do not insert an extra factor of :math:`\pi`.
    """
    x = numpy.asarray(x, dtype=float)
    y = numpy.asarray(y, dtype=float)
    return 0.5 * (
        numpy.exp(-numpy.pi * 1j * x * (1.0 + y)) * numpy.sinc(x * (1.0 - y))
        + numpy.exp(numpy.pi * 1j * x * (1.0 - y)) * numpy.sinc(x * (1.0 + y))
    )


def time_to_coalescence(frequency, chirp_mass, symmetric_mass_ratio, mode=2):
    r"""Time remaining to merger at a given frequency, to 3.5PN.

    This is the time-frequency relation that turns the frequency-domain response
    into a time-dependent one: the bin at frequency :math:`f` is emitted
    :math:`\tau(f)` seconds before merger, so the Earth's orientation there is
    the orientation at ``geocent_time - tau(f)``.

    Uses eq. (3.8b) of `arXiv:0907.0700 <https://arxiv.org/abs/0907.0700>`_.

    Parameters
    ----------
    frequency : array_like
        Frequency in Hz.
    chirp_mass : float
        Detector-frame chirp mass in solar masses.
    symmetric_mass_ratio : float
        :math:`\eta = m_1 m_2 / (m_1+m_2)^2`.
    mode : int
        Harmonic mode number. Mode ``m`` reaches a given frequency when the
        dominant quadrupole is at :math:`2f/m`, so ``tau_m(f) = tau_2(2f/m)``.

    Returns
    -------
    numpy.ndarray
        Seconds before coalescence, positive and decreasing with frequency.
    """
    frequency = numpy.asarray(frequency, dtype=float)
    if mode != 2:
        frequency = 2.0 * frequency / mode

    total_mass = chirp_mass * MTSUN_SI / (symmetric_mass_ratio ** (3.0 / 5.0))
    # Guard f = 0 (the DC bin): tau diverges there, and the bin is masked out of
    # every analysis anyway.
    with numpy.errstate(divide="ignore", invalid="ignore"):
        velocity = (numpy.pi * total_mass * frequency) ** (1.0 / 3.0)

        eta = symmetric_mass_ratio
        eta2 = eta * eta
        prefactor = 5.0 / 256 * total_mass / (eta * velocity**8)

        through_2p5 = (
            1.0
            + (743.0 / 252.0 + 11.0 / 3.0 * eta) * velocity**2
            - 32.0 / 5.0 * numpy.pi * velocity**3
            + (3058673.0 / 508032.0 + 5429.0 / 504.0 * eta + 617.0 / 72.0 * eta2)
            * velocity**4
            - (7729.0 / 252.0 - 13.0 / 3.0 * eta) * numpy.pi * velocity**5
        )
        order_3 = (
            -10052469856691.0 / 23471078400.0
            + 128.0 / 3.0 * numpy.pi**2
            + 6848.0 / 105.0 * numpy.euler_gamma
            + (3147553127.0 / 3048192.0 - 451.0 / 12.0 * numpy.pi**2) * eta
            - 15211.0 / 1728.0 * eta2
            + 25565.0 / 1296.0 * eta2 * eta
            + 3424.0 / 105.0 * numpy.log(16.0 * velocity**2)
        ) * velocity**6
        order_3p5 = (
            -15419335.0 / 127008.0 - 75703.0 / 756.0 * eta + 14809.0 / 378.0 * eta2
        ) * numpy.pi * velocity**7

        tau = prefactor * (through_2p5 + order_3 + order_3p5)

    return numpy.nan_to_num(tau, nan=0.0, posinf=0.0, neginf=0.0)


def chirp_mass_and_symmetric_mass_ratio(parameters):
    """Extract ``(chirp_mass, symmetric_mass_ratio)`` from a parameter dict.

    Accepts any of the mass parameterisations GWForge and bilby pass around:
    component masses, chirp mass with mass ratio, or chirp mass with symmetric
    mass ratio.

    Returns
    -------
    tuple of float, or None
        ``None`` when the masses cannot be determined -- e.g. a custom non-CBC
        source model -- in which case the caller falls back to a static response.
    """
    if "chirp_mass" in parameters and "symmetric_mass_ratio" in parameters:
        return float(parameters["chirp_mass"]), float(
            parameters["symmetric_mass_ratio"]
        )
    if "mass_1" in parameters and "mass_2" in parameters:
        mass_1 = float(parameters["mass_1"])
        mass_2 = float(parameters["mass_2"])
        total = mass_1 + mass_2
        eta = mass_1 * mass_2 / total**2
        return total * eta ** (3.0 / 5.0), eta
    if "chirp_mass" in parameters and "mass_ratio" in parameters:
        mass_ratio = float(parameters["mass_ratio"])
        eta = mass_ratio / (1.0 + mass_ratio) ** 2
        return float(parameters["chirp_mass"]), eta
    return None


def _arm_outer_products(interferometer):
    """Return ``(xx, yy)``, the half-outer-products of the arm unit vectors.

    ``xx - yy`` is exactly bilby's ``detector_tensor``; splitting it in two is
    what lets each arm carry its own frequency-dependent transfer function.
    """
    geometry = interferometer.geometry
    xx = 0.5 * numpy.einsum("i,j->ij", geometry.x, geometry.x)
    yy = 0.5 * numpy.einsum("i,j->ij", geometry.y, geometry.y)
    return xx, yy


def antenna_response(
    interferometer,
    ra,
    dec,
    geocent_time,
    psi,
    frequencies,
    start_time,
    times_to_coalescence=None,
    earth_rotation=True,
    finite_size=True,
):
    r"""Antenna pattern and arrival delay, per frequency bin.

    Parameters
    ----------
    interferometer : bilby.gw.detector.Interferometer
        Detector whose geometry is used. For ET this is one triangle arm; the
        three arms are separate ``Interferometer`` objects with their own
        vertices, so their responses differ correctly.
    ra, dec : float
        Sky position in radians.
    geocent_time : float
        GPS time of coalescence at the geocentre.
    psi : float
        Polarisation angle in radians, counter-clockwise about the propagation
        direction (bilby's convention).
    frequencies : array_like
        Frequencies in Hz at which to evaluate the response.
    start_time : float
        GPS start time of the data segment; the returned delay is measured from
        it, matching bilby's ``get_detector_response``.
    times_to_coalescence : array_like or None
        :math:`\tau(f)`, same length as ``frequencies``. Required when
        ``earth_rotation`` is True; ignored otherwise.
    earth_rotation : bool
        Evolve the antenna pattern and the arrival delay along :math:`t(f)`.
    finite_size : bool
        Apply the per-arm finite-light-travel-time transfer function.

    Returns
    -------
    tuple of (numpy.ndarray, numpy.ndarray, numpy.ndarray or float)
        ``(F_plus, F_cross, tau)``. The first two are complex when
        ``finite_size`` is True and real-valued otherwise; ``tau`` is the delay
        from ``start_time`` to arrival at this detector's vertex. The caller
        applies :math:`e^{-2\pi i f \tau}` -- keeping it separate leaves the
        rapidly oscillating factor available in closed form.
    """
    frequencies = numpy.asarray(frequencies, dtype=float)

    if earth_rotation:
        if times_to_coalescence is None:
            raise ValueError(
                "times_to_coalescence is required when earth_rotation is True"
            )
        # Linearise the sidereal rate over the signal: GMST is very nearly linear
        # in GPS time, and this avoids any 2*pi wrapping. The rate agrees with
        # 2*pi/sidereal_day to 3e-10 relative.
        gmst_at_coalescence = greenwich_mean_sidereal_time(geocent_time)
        day = 24.0 * 60.0 * 60.0
        rate = (greenwich_mean_sidereal_time(geocent_time + day) - gmst_at_coalescence) / day
        gmst = gmst_at_coalescence - rate * numpy.asarray(
            times_to_coalescence, dtype=float
        )
    else:
        gmst = numpy.full(len(frequencies), greenwich_mean_sidereal_time(geocent_time))

    gmst = numpy.mod(gmst, 2.0 * numpy.pi)
    theta, phi = ra_dec_to_theta_phi(ra, dec, gmst)
    theta = numpy.broadcast_to(theta, numpy.shape(phi))
    cos_phi, sin_phi = numpy.cos(phi), numpy.sin(phi)
    cos_theta, sin_theta = numpy.cos(theta), numpy.sin(theta)

    u = numpy.array([cos_phi * cos_theta, cos_theta * sin_phi, -sin_theta])
    v = numpy.array([-sin_phi, cos_phi, numpy.zeros_like(sin_phi)])
    m = -u * numpy.sin(psi) - v * numpy.cos(psi)
    n = -u * numpy.cos(psi) + v * numpy.sin(psi)
    omega = numpy.array([sin_theta * cos_phi, sin_theta * sin_phi, cos_theta])

    cross_term = numpy.einsum("ik,jk->ijk", m, n)
    polarization_plus = numpy.einsum("ik,jk->ijk", m, m) - numpy.einsum(
        "ik,jk->ijk", n, n
    )
    polarization_cross = cross_term + numpy.transpose(cross_term, axes=(1, 0, 2))

    if finite_size:
        arm_xx, arm_yy = _arm_outer_products(interferometer)
        plus_xx = numpy.einsum("ij,ijk->k", arm_xx, polarization_plus)
        plus_yy = numpy.einsum("ij,ijk->k", arm_yy, polarization_plus)
        cross_xx = numpy.einsum("ij,ijk->k", arm_xx, polarization_cross)
        cross_yy = numpy.einsum("ij,ijk->k", arm_yy, polarization_cross)

        # cos(angle) between each arm and the propagation direction.
        projection_x = -numpy.dot(omega.T, interferometer.geometry.x)
        projection_y = -numpy.dot(omega.T, interferometer.geometry.y)
        arm_length_in_wavelengths = (
            frequencies * interferometer.geometry.length * 1e3 / C_SI
        )
        transfer_x = finite_size_factor(arm_length_in_wavelengths, projection_x)
        transfer_y = finite_size_factor(arm_length_in_wavelengths, projection_y)

        response_plus = plus_xx * transfer_x - plus_yy * transfer_y
        response_cross = cross_xx * transfer_x - cross_yy * transfer_y
    else:
        detector_tensor = interferometer.geometry.detector_tensor
        response_plus = numpy.einsum("ij,ijk->k", detector_tensor, polarization_plus)
        response_cross = numpy.einsum("ij,ijk->k", detector_tensor, polarization_cross)

    vertex = interferometer.geometry.vertex
    if earth_rotation:
        vertex_delay = -numpy.dot(omega.T, vertex) / C_SI
    else:
        # Every column of omega is the same here, so the (len(frequencies), 3) @ (3,)
        # product would compute one value len(frequencies) times over just to keep the
        # first. The explicit three-term sum also mirrors bilby_cython's scalar
        # expression rather than going through a BLAS dgemv, whose accumulation order
        # varies between builds and shifts tau by an ulp.
        vertex_delay = (
            -(omega[0, 0] * vertex[0] + omega[1, 0] * vertex[1] + omega[2, 0] * vertex[2])
            / C_SI
        )
    tau = (geocent_time - start_time) + vertex_delay

    return response_plus, response_cross, tau


def detector_response(
    interferometer,
    waveform_polarizations,
    parameters,
    frequencies=None,
    earth_rotation=True,
    finite_size=True,):

    """Project polarizations onto one detector.

    Drop-in replacement for :meth:`bilby.gw.detector.Interferometer.get_detector_response`
    with the same masking and return contract: a full-length frequency-domain
    array, zero outside the interferometer's ``frequency_mask``.

    Parameters
    ----------
    interferometer : bilby.gw.detector.Interferometer
        Detector with strain data already set (its ``start_time`` is used).
    waveform_polarizations : dict
        ``{'plus': ..., 'cross': ...}`` frequency-domain polarizations.
    parameters : dict
        Source parameters. Needs ``ra``, ``dec``, ``psi``, ``geocent_time``, and
        -- when ``earth_rotation`` is True -- enough mass information to evaluate
        the time-frequency relation.
    frequencies : array_like or None
        Evaluate at these frequencies instead of the interferometer's own array.
        When given, no frequency masking is applied.
    earth_rotation, finite_size : bool
        Switch the two corrections off to recover bilby's response exactly.

    Returns
    -------
    numpy.ndarray
        Complex frequency-domain strain in this detector.
    """
    if frequencies is None:
        mask = interferometer.frequency_mask
        masked_frequencies = interferometer.frequency_array[mask]
    else:
        masked_frequencies = numpy.asarray(frequencies, dtype=float)
        mask = numpy.ones(len(masked_frequencies), dtype=bool)

    times_to_coalescence = None
    if earth_rotation:
        masses = chirp_mass_and_symmetric_mass_ratio(parameters)
        if masses is None:
            logging.warning(
                "%s: cannot determine the masses from the source parameters, so "
                "the time-frequency relation is unavailable; falling back to a "
                "static antenna pattern.",
                interferometer.name,
            )
            earth_rotation = False
        else:
            times_to_coalescence = time_to_coalescence(masked_frequencies, *masses)

    response_plus, response_cross, tau = antenna_response(
        interferometer,
        ra=parameters["ra"],
        dec=parameters["dec"],
        geocent_time=parameters["geocent_time"],
        psi=parameters["psi"],
        frequencies=masked_frequencies,
        start_time=interferometer.strain_data.start_time,
        times_to_coalescence=times_to_coalescence,
        earth_rotation=earth_rotation,
        finite_size=finite_size,
    )

    signal = numpy.zeros(len(mask), dtype=complex)
    signal[mask] = (
        waveform_polarizations["plus"][mask] * response_plus
        + waveform_polarizations["cross"][mask] * response_cross
    )
    signal[mask] = signal[mask] * numpy.exp(-1j * 2 * numpy.pi * tau * masked_frequencies)
    return signal


def inject_signal_with_response(
    interferometers,
    waveform_generator,
    parameters,
    earth_rotation=True,
    finite_size=True,
):
    """Inject a signal into a network using the frequency-dependent response.

    bilby's ``InterferometerList.inject_signal`` calls ``get_detector_response``
    internally and so cannot be reused. This mirrors
    ``Interferometer.inject_signal_from_waveform_polarizations``, including the
    ``meta_data`` it populates, so downstream code that reads
    ``ifo.meta_data['optimal_SNR']`` is unaffected.

    Parameters
    ----------
    interferometers : bilby.gw.detector.InterferometerList
        Network with strain data already set.
    waveform_generator : bilby.gw.WaveformGenerator
        Generator matched to the network's duration and sampling frequency.
    parameters : dict
        Injection parameters.
    earth_rotation, finite_size : bool
        Response options; see :func:`detector_response`.

    Returns
    -------
    dict
        The waveform polarizations that were injected.
    """
    polarizations = waveform_generator.frequency_domain_strain(parameters)

    for interferometer in interferometers:
        if not interferometer.strain_data.time_within_data(parameters["geocent_time"]):
            logging.warning(
                "Injecting signal outside segment: start_time=%s, merger time=%s.",
                interferometer.strain_data.start_time,
                parameters["geocent_time"],
            )
        signal = detector_response(
            interferometer,
            polarizations,
            parameters,
            earth_rotation=earth_rotation,
            finite_size=finite_size,
        )
        interferometer.strain_data.frequency_domain_strain += signal

        interferometer.meta_data["optimal_SNR"] = numpy.sqrt(
            interferometer.optimal_snr_squared(signal=signal)
        ).real
        interferometer.meta_data["matched_filter_SNR"] = (
            interferometer.matched_filter_snr(signal=signal)
        )
        interferometer.meta_data["parameters"] = parameters

    return polarizations
