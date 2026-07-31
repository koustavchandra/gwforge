r"""Validation of the frequency- and time-dependent antenna response.

Two reference implementations exist for this physics -- the ``bilby-x-g`` fork and
GWFast -- and they use different conventions, so neither can be treated as
ground truth on its own. The tests are therefore layered outward from claims that
cannot be wrong:

1. **Reduction to bilby.** Switching both corrections off must reproduce bilby's
   own ``get_detector_response`` to machine precision. This anchors every sign,
   the handedness of :math:`\psi`, and the delay convention against code GWForge
   already relies on. A convention error surfaces here, not in the GWFast test.
2. **Finite size, self-contained.** The closed-form arm transfer function is
   checked against a brute-force numerical integration of the round-trip light
   path, so it depends on no reference implementation at all.
3. **Earth rotation vs GWFast.** Only after 1 and 2 does the cross-code check
   mean anything, and it is set up to isolate geometry from the two codes'
   other differences (see :func:`test_earth_rotation_matches_gwfast`).
"""

import os

import numpy
import pytest

import bilby

from GWForge.ifo import antenna
from GWForge.ifo.detectors import Network

bilby.core.utils.setup_logger(log_level="error")

SPEED_OF_LIGHT = 299792458.0
GPS_TIME = 1893024018.0
NOISE_CURVES = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "GWForge",
    "ifo",
    "noise_curves",
)

#: A BNS-like source: long enough that the Earth turns appreciably during it.
CHIRP_MASS = 1.2
SYMMETRIC_MASS_RATIO = 0.2497


@pytest.fixture(scope="module")
def network():
    """CE40 plus the three ET arms, with strain data set."""
    interferometers = Network(ifos=["CE40", "ET"]).initialise_ifos()
    interferometers.set_strain_data_from_power_spectral_densities(
        sampling_frequency=2048.0, duration=8.0, start_time=GPS_TIME
    )
    return interferometers


# -- 1. the anchor: exact reduction to bilby ----------------------------------


@pytest.mark.parametrize(
    "ra, dec, psi",
    [
        (1.1, -0.4, 2.3),
        (0.0, 0.0, 0.0),
        (5.9, 1.3, 0.7),
        (3.2, -1.4, 3.0),
    ],
)
def test_reduces_to_bilby_when_both_corrections_are_off(network, ra, dec, psi):
    """With Earth rotation and finite size off, this *is* bilby's response.

    The tolerance is machine precision, not a loose 1e-6: anything looser would
    let a genuine convention error hide.
    """
    generator = numpy.random.default_rng(0)
    length = len(network[0].frequency_array)
    polarizations = {
        mode: generator.normal(size=length) + 1j * generator.normal(size=length)
        for mode in ("plus", "cross")
    }
    parameters = dict(
        ra=ra,
        dec=dec,
        psi=psi,
        geocent_time=GPS_TIME + 6.0,
        mass_1=1.4,
        mass_2=1.35,
    )

    for interferometer in network:
        reference = interferometer.get_detector_response(polarizations, parameters)
        mine = antenna.detector_response(
            interferometer,
            polarizations,
            parameters,
            earth_rotation=False,
            finite_size=False,
        )
        mask = interferometer.frequency_mask
        scale = numpy.abs(reference[mask]).max()
        assert numpy.abs(mine[mask] - reference[mask]).max() / scale < 1e-12
        # Outside the band bilby returns exactly zero; so must we.
        assert numpy.count_nonzero(mine[~mask]) == 0


def test_triangle_arms_have_distinct_responses(network):
    """ET's three arms must not be degenerate.

    bilby's ``TriangularInterferometer`` walks the vertices apart, so each arm has
    its own antenna pattern and its own arrival delay. If they collapsed to one
    response the null stream would be wrong.
    """
    arms = [ifo for ifo in network if ifo.name.startswith("ET")]
    assert len(arms) == 3
    responses = []
    for arm in arms:
        plus, cross, tau = antenna.antenna_response(
            arm,
            ra=1.1,
            dec=-0.4,
            geocent_time=GPS_TIME + 6.0,
            psi=2.3,
            frequencies=numpy.array([100.0]),
            start_time=GPS_TIME,
            earth_rotation=False,
            finite_size=False,
        )
        responses.append((plus[0], cross[0]))
    for first in range(3):
        for second in range(first + 1, 3):
            assert not numpy.isclose(responses[first][0], responses[second][0])


# -- 2. finite size, checked against first principles -------------------------


@pytest.mark.parametrize("y", numpy.linspace(-1.0, 1.0, 9))
def test_finite_size_factor_long_wavelength_limit(y):
    """``D(0, y) == 1``: the long-wavelength limit, exactly, for every direction."""
    assert abs(complex(antenna.finite_size_factor(0.0, y)) - 1.0) < 1e-15


def _round_trip_light_path(x, y, samples=400001):
    r"""Brute-force one-arm response, from the light path alone.

    Averages :math:`h_{aa}` over both legs of the round trip and references the
    result to the round-trip midpoint. Convention:
    :math:`h \sim e^{-2\pi i f(t - \hat\Omega\cdot\vec{x}/c)}`,
    :math:`y = -\hat\Omega\cdot\hat{a}`, with the photon leaving the vertex at
    :math:`t - 2L/c` and returning at :math:`t`.

    This shares no code with :func:`~GWForge.ifo.antenna.finite_size_factor`, so
    agreement validates the closed form independently of any reference package.
    """
    fraction = numpy.linspace(0.0, 1.0, samples)
    outbound = numpy.exp(-2j * numpy.pi * x * (-2.0 + fraction + fraction * y))
    inbound = numpy.exp(-2j * numpy.pi * x * (-fraction + fraction * y))
    average = 0.5 * (
        numpy.trapezoid(outbound, fraction) + numpy.trapezoid(inbound, fraction)
    )
    return average * numpy.exp(-2j * numpy.pi * x)


@pytest.mark.parametrize("x", [0.0, 0.01, 0.1, 0.5, 1.0, 2.0])
@pytest.mark.parametrize("y", [-1.0, -0.37, 0.0, 0.63, 1.0])
def test_finite_size_factor_matches_the_light_path_integral(x, y):
    closed_form = complex(antenna.finite_size_factor(x, y))
    brute_force = _round_trip_light_path(x, y)
    assert abs(closed_form - brute_force) < 1e-9


def test_finite_size_leading_order_is_a_light_travel_phase():
    r"""To first order :math:`D \approx 1 - i\pi x y`.

    Worth pinning explicitly, because it says the leading finite-size effect is a
    *phase* -- the extra light travel to the arm midpoint -- and not an amplitude
    change. That is why the correction shows up as linear in :math:`f` rather
    than quadratic, and why much of it is degenerate with the coalescence time.
    """
    y = 0.7
    # The imaginary part is the first-order term itself, so it is accurate to
    # the relative size of the next order.
    for x in (1e-5, 1e-4):
        value = complex(antenna.finite_size_factor(x, y))
        assert value.imag == pytest.approx(-numpy.pi * x * y, rel=1e-6)

    # The residual after removing the first-order term must be O(x^2): doubling
    # x must quadruple it. This is what shows the expansion is right, rather
    # than merely small.
    def residual(x):
        return abs(complex(antenna.finite_size_factor(x, y)) - (1.0 - 1j * numpy.pi * x * y))

    assert residual(2e-4) / residual(1e-4) == pytest.approx(4.0, rel=1e-3)


@pytest.mark.parametrize("arm_length_km, name", [(40.0, "CE40"), (10.0, "ET")])
def test_finite_size_null_at_the_free_spectral_range(arm_length_km, name):
    r"""The arm response vanishes at :math:`f = c/2L` for a transverse source.

    At :math:`y = 0`, :math:`D = \mathrm{sinc}(x)\cos(\pi x)`, whose first zero is
    at :math:`x = 1/2`. That is 3.75 kHz for CE's 40 km arms -- inside the
    analysis band, which is precisely why this correction matters for XG.
    """
    free_spectral_range = SPEED_OF_LIGHT / (2.0 * arm_length_km * 1e3)
    x = free_spectral_range * arm_length_km * 1e3 / SPEED_OF_LIGHT
    assert abs(complex(antenna.finite_size_factor(x, 0.0))) < 1e-12
    # ... and it is non-trivial well below that.
    assert abs(complex(antenna.finite_size_factor(0.5 * x, 0.0))) < 0.75


def test_finite_size_correction_grows_with_frequency(network):
    """Small in the low band, O(10%) near the free spectral range."""
    interferometer = network[0]
    frequencies = numpy.array([10.0, 100.0, 1000.0])
    static = antenna.antenna_response(
        interferometer,
        1.1,
        -0.4,
        GPS_TIME + 6.0,
        2.3,
        frequencies,
        GPS_TIME,
        earth_rotation=False,
        finite_size=False,
    )[0]
    corrected = antenna.antenna_response(
        interferometer,
        1.1,
        -0.4,
        GPS_TIME + 6.0,
        2.3,
        frequencies,
        GPS_TIME,
        earth_rotation=False,
        finite_size=True,
    )[0]
    departure = numpy.abs(corrected - static) / numpy.abs(static)
    assert departure[0] < 1e-3
    assert numpy.all(numpy.diff(departure) > 0)


# -- 3. Earth rotation, cross-checked against GWFast --------------------------


def test_earth_rotation_uses_the_sidereal_rate():
    """The Earth turns once per sidereal day, not per solar day.

    GWFast advances its pattern by ``2*pi * t/86400`` -- a solar day -- which over
    a one-hour inspiral is ~7e-4 rad adrift of the true sidereal rotation. This
    asserts GWForge uses the sidereal rate, and is the reason
    :func:`test_earth_rotation_matches_gwfast` equalises the two codes' Earth
    orientations rather than comparing them at face value.
    """
    from bilby_cython.geometry import greenwich_mean_sidereal_time

    day = 86400.0
    rate = (
        greenwich_mean_sidereal_time(GPS_TIME + day)
        - greenwich_mean_sidereal_time(GPS_TIME)
    ) / day
    sidereal_day = 86164.0905
    assert rate == pytest.approx(2 * numpy.pi / sidereal_day, rel=1e-8)
    assert abs(rate - 2 * numpy.pi / day) / rate > 1e-3


def test_time_to_coalescence_matches_gwfast_tau_star():
    """Same 3.5PN time-frequency relation as GWFast (arXiv:0907.0700 eq. 3.8b).

    Using an identical ``tau(f)`` is what makes the pattern comparison a test of
    *geometry*, rather than of two different time-frequency relations.
    """
    gwfast_waveforms = pytest.importorskip("gwfast.waveforms")
    frequencies = numpy.array([5.0, 20.0, 100.0, 1000.0])
    reference = gwfast_waveforms.IMRPhenomD().tau_star(
        frequencies, Mc=CHIRP_MASS, eta=SYMMETRIC_MASS_RATIO
    )
    mine = antenna.time_to_coalescence(
        frequencies, CHIRP_MASS, SYMMETRIC_MASS_RATIO
    )
    numpy.testing.assert_allclose(mine, reference, rtol=1e-13)


def test_time_to_coalescence_mode_scaling():
    """Mode ``m`` reaches frequency ``f`` when the quadrupole is at ``2f/m``."""
    frequencies = numpy.array([20.0, 100.0])
    quadrupole = antenna.time_to_coalescence(
        frequencies, CHIRP_MASS, SYMMETRIC_MASS_RATIO, mode=2
    )
    fourth = antenna.time_to_coalescence(
        frequencies, CHIRP_MASS, SYMMETRIC_MASS_RATIO, mode=4
    )
    expected = antenna.time_to_coalescence(
        frequencies / 2.0, CHIRP_MASS, SYMMETRIC_MASS_RATIO, mode=2
    )
    numpy.testing.assert_allclose(fourth, expected, rtol=1e-13)
    assert numpy.all(fourth > quadrupole)


def test_earth_rotation_matches_gwfast(network):
    r"""Time-varying :math:`F_{+,\times}` against GWFast, geometry isolated.

    Convention mapping (the ``xax`` offset was verified against H1, L1 and Virgo,
    where GWFast's ``xax`` minus bilby's ``xarm_azimuth`` is exactly 45 degrees):

    * ``det_xax = xarm_azimuth + 45`` -- GWFast orients by the arm *bisector*.
    * ``theta = pi/2 - dec``, ``phi = ra``.
    * GWFast's ``t`` is a rotation angle in units of days, entering its pattern
      functions only as ``2*pi*t``. Passing ``gmst/(2*pi)`` therefore places both
      codes at *identical* Earth orientations, which removes GWFast's solar-day
      approximation (see :func:`test_earth_rotation_uses_the_sidereal_rate`) and
      leaves a pure geometry comparison.
    * Only the real :math:`F_{+,\times}` amplitudes are compared: GWFast writes
      :math:`h \sim e^{+i\Psi}` where bilby writes :math:`e^{-i\ldots}`, so the
      complex responses differ by a conjugation that says nothing about geometry.

    CE40 is used rather than ET because GWFast models a spherical Earth and
    ignores detector elevation (CE40's ``.ifo`` has ``elevation = 0``), and
    because GWFast places its triangle's three arms at a single vertex whereas
    bilby walks them apart.
    """
    pytest.importorskip("gwfast")
    # gwfast 1.1.2 calls scipy.integrate.cumtrapz, removed in modern scipy.
    import scipy.integrate

    if not hasattr(scipy.integrate, "cumtrapz"):
        scipy.integrate.cumtrapz = scipy.integrate.cumulative_trapezoid
    from bilby_cython.geometry import greenwich_mean_sidereal_time
    from gwfast.signal import GWSignal
    from gwfast.waveforms import IMRPhenomD

    interferometer = network[0]
    assert interferometer.name == "CE40"

    signal = GWSignal(
        IMRPhenomD(),
        psd_path=os.path.join(NOISE_CURVES, "CE40-asd.txt"),
        detector_shape="L",
        det_lat=interferometer.latitude,
        det_long=interferometer.longitude,
        det_xax=interferometer.xarm_azimuth + 45.0,
        is_ASD=True,
        useEarthMotion=True,
        fmin=5.0,
        fmax=2048.0,
        verbose=False,
    )

    ra, dec, psi = 1.1, -0.4, 2.3
    frequencies = numpy.array([7.0, 20.0, 50.0, 200.0, 1000.0])
    times_to_coalescence = antenna.time_to_coalescence(
        frequencies, CHIRP_MASS, SYMMETRIC_MASS_RATIO
    )
    # The lowest bin sits ~45 minutes before merger: the Earth really has turned.
    assert times_to_coalescence[0] > 40 * 60

    sidereal = numpy.array(
        [greenwich_mean_sidereal_time(GPS_TIME - t) for t in times_to_coalescence]
    )
    reference_plus, reference_cross = signal._PatternFunction(
        numpy.pi / 2 - dec, ra, sidereal / (2 * numpy.pi), psi
    )

    plus, cross, _ = antenna.antenna_response(
        interferometer,
        ra=ra,
        dec=dec,
        geocent_time=GPS_TIME,
        psi=psi,
        frequencies=frequencies,
        start_time=GPS_TIME - 1000.0,
        times_to_coalescence=times_to_coalescence,
        earth_rotation=True,
        finite_size=False,
    )

    numpy.testing.assert_allclose(plus, reference_plus, rtol=0, atol=1e-8)
    numpy.testing.assert_allclose(cross, reference_cross, rtol=0, atol=1e-8)
    # The pattern must genuinely vary across the band, or the comparison is
    # vacuous -- it would just be re-checking the static limit five times.
    assert numpy.ptp(plus) > 1e-2


def test_earth_rotation_changes_the_response_over_a_long_signal(network):
    """A BNS spans enough Earth rotation to move the pattern measurably."""
    interferometer = network[0]
    frequencies = numpy.array([7.0, 1000.0])
    times_to_coalescence = antenna.time_to_coalescence(
        frequencies, CHIRP_MASS, SYMMETRIC_MASS_RATIO
    )
    rotating = antenna.antenna_response(
        interferometer,
        1.1,
        -0.4,
        GPS_TIME,
        2.3,
        frequencies,
        GPS_TIME - 1000.0,
        times_to_coalescence=times_to_coalescence,
        earth_rotation=True,
        finite_size=False,
    )[0]
    static = antenna.antenna_response(
        interferometer,
        1.1,
        -0.4,
        GPS_TIME,
        2.3,
        frequencies,
        GPS_TIME - 1000.0,
        earth_rotation=False,
        finite_size=False,
    )[0]
    # The highest bin is essentially at merger, so it must match the static value.
    assert rotating[-1] == pytest.approx(static[-1], rel=1e-6)
    # The lowest bin is ~45 min earlier, and must not.
    assert abs(rotating[0] - static[0]) > 1e-2


# -- fallbacks and wiring -----------------------------------------------------


def test_falls_back_to_static_without_masses(network, caplog):
    """A source model with no masses cannot supply t(f); say so and degrade."""
    interferometer = network[0]
    length = len(interferometer.frequency_array)
    polarizations = {
        mode: numpy.ones(length, dtype=complex) for mode in ("plus", "cross")
    }
    parameters = dict(ra=1.1, dec=-0.4, psi=2.3, geocent_time=GPS_TIME + 6.0)

    fallback = antenna.detector_response(
        interferometer, polarizations, parameters, earth_rotation=True, finite_size=False
    )
    static = antenna.detector_response(
        interferometer,
        polarizations,
        parameters,
        earth_rotation=False,
        finite_size=False,
    )
    numpy.testing.assert_allclose(fallback, static, rtol=0, atol=0)


@pytest.mark.parametrize(
    "parameters, expected_chirp_mass",
    [
        ({"chirp_mass": 1.2, "symmetric_mass_ratio": 0.25}, 1.2),
        ({"mass_1": 1.4, "mass_2": 1.4}, 2.8 * 0.25**0.6),
        ({"chirp_mass": 1.2, "mass_ratio": 1.0}, 1.2),
    ],
)
def test_mass_resolver(parameters, expected_chirp_mass):
    resolved = antenna.chirp_mass_and_symmetric_mass_ratio(parameters)
    assert resolved is not None
    assert resolved[0] == pytest.approx(expected_chirp_mass, rel=1e-12)
    assert 0 < resolved[1] <= 0.25


def test_mass_resolver_returns_none_when_unavailable():
    assert antenna.chirp_mass_and_symmetric_mass_ratio({"ra": 1.0}) is None


def test_inject_signal_with_response_populates_meta_data():
    """The injection helper must be a drop-in for bilby's, ``meta_data`` included."""
    interferometers = Network(ifos=["CE40"]).initialise_ifos()
    duration, sampling_frequency = 8.0, 2048.0
    interferometers.set_strain_data_from_power_spectral_densities(
        sampling_frequency=sampling_frequency,
        duration=duration,
        start_time=GPS_TIME,
    )
    waveform_generator = bilby.gw.WaveformGenerator(
        duration=duration,
        sampling_frequency=sampling_frequency,
        frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
        parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
        waveform_arguments={
            "waveform_approximant": "IMRPhenomD",
            "reference_frequency": 20.0,
            "minimum_frequency": 7.0,
        },
    )
    parameters = dict(
        mass_1=36.0,
        mass_2=29.0,
        a_1=0.0,
        a_2=0.0,
        tilt_1=0.0,
        tilt_2=0.0,
        phi_12=0.0,
        phi_jl=0.0,
        luminosity_distance=1000.0,
        theta_jn=0.4,
        psi=2.659,
        phase=1.3,
        geocent_time=GPS_TIME + 6.0,
        ra=1.375,
        dec=-1.2108,
    )

    antenna.inject_signal_with_response(
        interferometers, waveform_generator, parameters, earth_rotation=True
    )
    for interferometer in interferometers:
        assert interferometer.meta_data["optimal_SNR"] > 1.0
        assert numpy.isfinite(interferometer.meta_data["matched_filter_SNR"])
        assert interferometer.meta_data["parameters"] is parameters


def test_response_flags_shift_the_snr_by_a_small_amount():
    """Turning the corrections on must move the SNR -- but not wildly.

    A BBH at these masses is short and band-limited, so Earth rotation is
    negligible and finite size is a percent-level effect. A large shift would
    mean a normalisation error rather than new physics.
    """
    parameters = dict(
        mass_1=36.0,
        mass_2=29.0,
        a_1=0.0,
        a_2=0.0,
        tilt_1=0.0,
        tilt_2=0.0,
        phi_12=0.0,
        phi_jl=0.0,
        luminosity_distance=1000.0,
        theta_jn=0.4,
        psi=2.659,
        phase=1.3,
        geocent_time=GPS_TIME + 6.0,
        ra=1.375,
        dec=-1.2108,
    )

    def optimal_snr(**flags):
        interferometers = Network(ifos=["CE40"]).initialise_ifos()
        interferometers.set_strain_data_from_power_spectral_densities(
            sampling_frequency=2048.0, duration=8.0, start_time=GPS_TIME
        )
        waveform_generator = bilby.gw.WaveformGenerator(
            duration=8.0,
            sampling_frequency=2048.0,
            frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
            parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
            waveform_arguments={
                "waveform_approximant": "IMRPhenomD",
                "reference_frequency": 20.0,
                "minimum_frequency": 7.0,
            },
        )
        antenna.inject_signal_with_response(
            interferometers, waveform_generator, parameters, **flags
        )
        return interferometers[0].meta_data["optimal_SNR"]

    baseline = optimal_snr(earth_rotation=False, finite_size=False)
    full = optimal_snr(earth_rotation=True, finite_size=True)
    assert baseline > 1.0
    assert full == pytest.approx(baseline, rel=0.05)
