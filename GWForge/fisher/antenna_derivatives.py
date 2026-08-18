r"""Closed-form derivatives of the detector response geometry.

``ra``, ``dec`` and ``psi`` reach the strain only through the detector
projection, never through the waveform, so their derivatives can be written down
exactly. This module differentiates the same expressions
:func:`GWForge.ifo.antenna.antenna_response` evaluates -- the ``u``, ``v``, ``m``
and ``n`` triad, the propagation direction, the finite-size arm transfer
function and the vertex delay -- with respect to the internal angles
:math:`(\phi, \theta, \psi)`, and chain-rules them onto the source parameters.

That holds under the frequency-dependent, Earth-rotating response too, which
changes only *how* :math:`F_{+,\times}` and :math:`\tau` are built. Freezing the
polarizations and differencing the response alone reproduces the full derivative
bit-identically; ``tests/test_fisher.py`` asserts exactly that.

The work is done on scalars rather than on :math:`3\times3\times n`
polarization tensors: contracting an arm with the polarization basis reduces to
:math:`(\hat a\cdot m)` and :math:`(\hat a\cdot n)`, so only those and their
derivatives are ever needed.
"""

import numpy
from bilby.core.utils import speed_of_light
from bilby_cython.geometry import greenwich_mean_sidereal_time

_SMALL = 1e-4


def _sinc_derivatives(z):
    r"""Return ``(sinc, sinc', sinc'')`` for :math:`\mathrm{sinc}(z)`.

    Uses numpy's convention :math:`\mathrm{sinc}(z) = \sin(\pi z)/(\pi z)`, to
    match :func:`GWForge.ifo.antenna.finite_size_factor`.

    Both derivatives are singular-looking at ``z = 0`` (they contain ``1/z`` and
    ``1/z**2``) while being perfectly finite there, so small arguments are
    evaluated from the Taylor series instead.

    Parameters
    ----------
    z : array_like

    Returns
    -------
    tuple of numpy.ndarray
    """
    z = numpy.asarray(z, dtype=float)
    pi = numpy.pi
    value = numpy.sinc(z)

    small = numpy.abs(z) < _SMALL
    safe = numpy.where(small, 1.0, z)

    cosine = numpy.cos(pi * safe)
    first = (cosine - numpy.sinc(safe)) / safe
    second = (
        -pi * safe * numpy.sin(pi * safe) + 2.0 * numpy.sinc(safe) - 2.0 * cosine
    ) / safe**2

    series_first = -(pi**2) * z / 3.0 + pi**4 * z**3 / 30.0
    series_second = -(pi**2) / 3.0 + pi**4 * z**2 / 10.0

    return (
        value,
        numpy.where(small, series_first, first),
        numpy.where(small, series_second, second),
    )


def _transfer_derivatives(x, y):
    r"""Arm transfer function :math:`D(x, y)` and its first two ``y`` derivatives.

    :math:`D` is :func:`GWForge.ifo.antenna.finite_size_factor`; ``y`` is the
    cosine between the arm and the propagation direction, which is where the sky
    position enters. ``x = fL/c`` does not depend on any source parameter, so no
    ``x`` derivatives are needed.

    Parameters
    ----------
    x : array_like
        Arm length in wavelengths.
    y : array_like
        :math:`-\hat\Omega\cdot\hat a`.

    Returns
    -------
    tuple of numpy.ndarray
        ``(D, dD/dy, d2D/dy2)``, complex.
    """
    x = numpy.asarray(x, dtype=float)
    y = numpy.asarray(y, dtype=float)
    pi = numpy.pi

    minus, plus = x * (1.0 - y), x * (1.0 + y)
    sinc_minus, dsinc_minus, ddsinc_minus = _sinc_derivatives(minus)
    sinc_plus, dsinc_plus, ddsinc_plus = _sinc_derivatives(plus)

    exp_minus = numpy.exp(-pi * 1j * x * (1.0 + y))
    exp_plus = numpy.exp(pi * 1j * x * (1.0 - y))

    value = 0.5 * (exp_minus * sinc_minus + exp_plus * sinc_plus)

    first = 0.5 * (
        -1j * pi * x * (exp_minus * sinc_minus + exp_plus * sinc_plus)
        - x * exp_minus * dsinc_minus
        + x * exp_plus * dsinc_plus
    )

    second = 0.5 * (
        -(pi**2) * x**2 * (exp_minus * sinc_minus + exp_plus * sinc_plus)
        + 2j * pi * x**2 * exp_minus * dsinc_minus
        - 2j * pi * x**2 * exp_plus * dsinc_plus
        + x**2 * (exp_minus * ddsinc_minus + exp_plus * ddsinc_plus)
    )

    return value, first, second


def _basis_derivatives(theta, phi):
    r"""Sky basis triad and propagation direction, with first and second derivatives.

    Mirrors the ``u``, ``v`` and ``omega`` built in
    :func:`GWForge.ifo.antenna.antenna_response`, so that differentiating here is
    differentiating exactly what the response evaluates.

    Parameters
    ----------
    theta, phi : array_like
        :math:`\theta = \pi/2 - \delta` and :math:`\phi = \alpha - \mathrm{gmst}`,
        broadcast to a common shape.

    Returns
    -------
    dict
        Maps ``"u"``, ``"v"``, ``"omega"`` and the derivative keys ``"_phi"``,
        ``"_theta"``, ``"_phiphi"``, ``"_phitheta"``, ``"_thetatheta"`` to
        ``(3, n)`` arrays.
    """
    theta, phi = numpy.broadcast_arrays(
        numpy.asarray(theta, dtype=float), numpy.asarray(phi, dtype=float)
    )
    cos_phi, sin_phi = numpy.cos(phi), numpy.sin(phi)
    cos_theta, sin_theta = numpy.cos(theta), numpy.sin(theta)
    zero = numpy.zeros_like(phi)

    def stack(first, second, third):
        return numpy.array([first, second, third])

    basis = {
        "u": stack(cos_theta * cos_phi, cos_theta * sin_phi, -sin_theta),
        "u_phi": stack(-cos_theta * sin_phi, cos_theta * cos_phi, zero),
        "u_theta": stack(-sin_theta * cos_phi, -sin_theta * sin_phi, -cos_theta),
        "u_phiphi": stack(-cos_theta * cos_phi, -cos_theta * sin_phi, zero),
        "u_phitheta": stack(sin_theta * sin_phi, -sin_theta * cos_phi, zero),
        "u_thetatheta": stack(-cos_theta * cos_phi, -cos_theta * sin_phi, sin_theta),
        "v": stack(-sin_phi, cos_phi, zero),
        "v_phi": stack(-cos_phi, -sin_phi, zero),
        "v_theta": stack(zero, zero, zero),
        "v_phiphi": stack(sin_phi, -cos_phi, zero),
        "v_phitheta": stack(zero, zero, zero),
        "v_thetatheta": stack(zero, zero, zero),
        "omega": stack(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta),
        "omega_phi": stack(-sin_theta * sin_phi, sin_theta * cos_phi, zero),
        "omega_theta": stack(cos_theta * cos_phi, cos_theta * sin_phi, -sin_theta),
        "omega_phiphi": stack(-sin_theta * cos_phi, -sin_theta * sin_phi, zero),
        "omega_phitheta": stack(-cos_theta * sin_phi, cos_theta * cos_phi, zero),
        "omega_thetatheta": stack(
            -sin_theta * cos_phi, -sin_theta * sin_phi, -cos_theta
        ),
    }
    return basis


#: The nine derivative labels tracked through the antenna calculation: the
#: value, three first derivatives and five independent second derivatives with
#: respect to the internal angles ``(phi, theta, psi)``.
_LABELS = (
    "",
    "_phi",
    "_theta",
    "_psi",
    "_phiphi",
    "_phitheta",
    "_thetatheta",
    "_phipsi",
    "_thetapsi",
    "_psipsi",
)


def _polarization_projections(basis, psi, arm):
    r"""Project the polarization triad onto one arm, with all derivatives.

    Rather than forming :math:`3\times3\times n` polarization tensors and
    contracting them, note that

    .. math::

       \tfrac{1}{2}\hat a\otimes\hat a : e^+ = \tfrac{1}{2}\big[(\hat a\cdot m)^2
       - (\hat a\cdot n)^2\big],
       \qquad
       \tfrac{1}{2}\hat a\otimes\hat a : e^\times = (\hat a\cdot m)(\hat a\cdot n),

    so everything reduces to the two scalars :math:`\hat a\cdot m` and
    :math:`\hat a\cdot n` and their derivatives.

    Parameters
    ----------
    basis : dict
        Output of :func:`_basis_derivatives`.
    psi : float
        Polarisation angle.
    arm : numpy.ndarray
        Unit vector along the arm.

    Returns
    -------
    tuple of dict
        ``(plus, cross)``, each keyed by the labels in ``_LABELS``.
    """
    sin_psi, cos_psi = numpy.sin(psi), numpy.cos(psi)

    def m_of(u_part, v_part):
        return -arm @ u_part * sin_psi - arm @ v_part * cos_psi

    def n_of(u_part, v_part):
        return -arm @ u_part * cos_psi + arm @ v_part * sin_psi

    # m = -u sin(psi) - v cos(psi) and n = -u cos(psi) + v sin(psi), so
    # d/d(psi) turns m into n and n into -m; the (phi, theta) derivatives just
    # carry through the same linear combination.
    m = {
        "": m_of(basis["u"], basis["v"]),
        "_phi": m_of(basis["u_phi"], basis["v_phi"]),
        "_theta": m_of(basis["u_theta"], basis["v_theta"]),
        "_psi": n_of(basis["u"], basis["v"]),
        "_phiphi": m_of(basis["u_phiphi"], basis["v_phiphi"]),
        "_phitheta": m_of(basis["u_phitheta"], basis["v_phitheta"]),
        "_thetatheta": m_of(basis["u_thetatheta"], basis["v_thetatheta"]),
        "_phipsi": n_of(basis["u_phi"], basis["v_phi"]),
        "_thetapsi": n_of(basis["u_theta"], basis["v_theta"]),
    }
    m["_psipsi"] = -m[""]

    n = {
        "": n_of(basis["u"], basis["v"]),
        "_phi": n_of(basis["u_phi"], basis["v_phi"]),
        "_theta": n_of(basis["u_theta"], basis["v_theta"]),
        "_psi": -m[""],
        "_phiphi": n_of(basis["u_phiphi"], basis["v_phiphi"]),
        "_phitheta": n_of(basis["u_phitheta"], basis["v_phitheta"]),
        "_thetatheta": n_of(basis["u_thetatheta"], basis["v_thetatheta"]),
        "_phipsi": -m["_phi"],
        "_thetapsi": -m["_theta"],
    }
    n["_psipsi"] = -n[""]

    plus = _product_rule(m, m, 0.5)
    minus = _product_rule(n, n, 0.5)
    plus = {key: plus[key] - minus[key] for key in _LABELS}
    cross = _product_rule(m, n, 1.0)
    return plus, cross


def _product_rule(left, right, factor):
    """Symmetrised product of two derivative dictionaries.

    Returns ``factor * left * right`` differentiated up to second order in the
    internal angles, which is the only place the product rule appears in this
    module.

    Parameters
    ----------
    left, right : dict
        Keyed by the labels in ``_LABELS``.
    factor : float
        Overall multiplier.

    Returns
    -------
    dict
    """
    product = {
        key: factor * (left[""] * right[key] + left[key] * right[""]) for key in _LABELS
    }
    product[""] = factor * left[""] * right[""]

    for key, first, second in (
        ("_phiphi", "_phi", "_phi"),
        ("_phitheta", "_phi", "_theta"),
        ("_thetatheta", "_theta", "_theta"),
        ("_phipsi", "_phi", "_psi"),
        ("_thetapsi", "_theta", "_psi"),
        ("_psipsi", "_psi", "_psi"),
    ):
        product[key] = factor * (
            left[""] * right[key]
            + left[key] * right[""]
            + left[first] * right[second]
            + left[second] * right[first]
        )
    return product


def _scale(mapping, factor):
    """Multiply every entry of a derivative dictionary by a scalar or array."""
    return {key: mapping[key] * factor for key in _LABELS}


def _combine(left, right, sign=1.0):
    """Add or subtract two derivative dictionaries entry by entry."""
    return {key: left[key] + sign * right[key] for key in _LABELS}


class AntennaDerivatives:
    r"""Exact derivatives of :math:`F_+`, :math:`F_\times` and :math:`\tau`.

    Differentiates the same expressions
    :func:`GWForge.ifo.antenna.antenna_response` evaluates, with respect to the
    internal angles :math:`(\phi, \theta, \psi)`, then chain-rules those onto the
    source parameters ``ra``, ``dec``, ``psi`` and ``geocent_time``:

    .. math::

       \frac{\partial}{\partial\alpha} = \frac{\partial}{\partial\phi},
       \qquad
       \frac{\partial}{\partial\delta} = -\frac{\partial}{\partial\theta},
       \qquad
       \frac{\partial}{\partial t_c} \supset
           -\dot{\mathrm{gmst}}\,\frac{\partial}{\partial\phi},

    the last because :math:`\phi = \alpha - \mathrm{gmst}(t_c)`. So the
    Earth-rotation part of the ``geocent_time`` derivative comes free from the
    sky partials. ``geocent_time`` additionally enters :math:`\tau` explicitly
    with unit derivative, which is the dominant :math:`-2\pi i f h` term.

    Attributes
    ----------
    plus, cross : dict
        :math:`F_+` and :math:`F_\times` and their derivatives, keyed by source
        parameter name (``"ra"``) or sorted pair (``("dec", "ra")``).
    delay : dict
        Same for :math:`\tau`.
    """

    def __init__(
        self,
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
        """
        Parameters
        ----------
        interferometer : bilby.gw.detector.Interferometer
            Detector geometry. For ET this is one arm of the triangle.
        ra, dec, geocent_time, psi : float
            Source parameters, radians and GPS seconds.
        frequencies : array_like
            Frequencies at which to evaluate, in Hz.
        start_time : float
            Segment start time; the delay is measured from it, matching bilby.
        times_to_coalescence : array_like or None
            :math:`\\tau(f)`; required when ``earth_rotation`` is True.
        earth_rotation, finite_size : bool
            Must match the settings used for the response itself.
        """
        frequencies = numpy.asarray(frequencies, dtype=float)
        day = 24.0 * 60.0 * 60.0
        gmst_at_coalescence = greenwich_mean_sidereal_time(geocent_time)
        # d(gmst)/d(t_c). antenna.py linearises the sidereal rate over the
        # signal for the same reason: gmst is very nearly linear in GPS time and
        # this sidesteps 2*pi wrapping.
        self.sidereal_rate = (
            greenwich_mean_sidereal_time(geocent_time + day) - gmst_at_coalescence
        ) / day

        if earth_rotation:
            if times_to_coalescence is None:
                raise ValueError(
                    "times_to_coalescence is required when earth_rotation is True"
                )
            gmst = gmst_at_coalescence - self.sidereal_rate * numpy.asarray(
                times_to_coalescence, dtype=float
            )
        else:
            gmst = numpy.full(len(frequencies), gmst_at_coalescence)

        # Wrap gmst *before* subtracting, exactly as antenna.py does.
        # greenwich_mean_sidereal_time returns an unwrapped angle -- of order
        # 3.6e4 rad here -- so forming ``ra - gmst`` first and wrapping after
        # quantises phi at ~7e-12 rad and destroys five digits of any
        # finite-difference check against this class.
        gmst = numpy.mod(gmst, 2.0 * numpy.pi)
        theta = numpy.pi / 2.0 - dec
        phi = ra - gmst
        basis = _basis_derivatives(numpy.full_like(phi, theta), phi)

        geometry = interferometer.geometry
        plus, cross, delay = self._responses(
            basis, psi, geometry, frequencies, finite_size
        )
        # antenna.py measures the delay from the segment start, not the vertex.
        delay[""] = delay[""] + (geocent_time - start_time)
        self._store(plus, cross, delay)

    @staticmethod
    def _responses(basis, psi, geometry, frequencies, finite_size):
        """Assemble the two pattern functions and the delay, with derivatives."""
        arm_x, arm_y = geometry.x, geometry.y
        plus_x, cross_x = _polarization_projections(basis, psi, arm_x)
        plus_y, cross_y = _polarization_projections(basis, psi, arm_y)

        if finite_size:
            arm_length_in_wavelengths = frequencies * geometry.length * 1e3 / speed_of_light
            transfer_x = _arm_transfer(basis, arm_x, arm_length_in_wavelengths)
            transfer_y = _arm_transfer(basis, arm_y, arm_length_in_wavelengths)
        else:
            # The long-wavelength limit is D = 1 with no angular dependence, so
            # the same product rule below collapses to bilby's detector tensor.
            transfer_x = transfer_y = _unit_derivatives(basis["u"].shape[1:])

        plus = _combine(
            _product_rule(plus_x, transfer_x, 1.0),
            _product_rule(plus_y, transfer_y, 1.0),
            sign=-1.0,
        )
        cross = _combine(
            _product_rule(cross_x, transfer_x, 1.0),
            _product_rule(cross_y, transfer_y, 1.0),
            sign=-1.0,
        )

        vertex = geometry.vertex
        delay = {
            key: -(vertex @ basis["omega" + key]) / speed_of_light
            if ("omega" + key) in basis
            else numpy.zeros_like(basis["omega"][0])
            for key in _LABELS
        }
        return plus, cross, delay

    def _store(self, plus, cross, delay):
        """Chain-rule the internal-angle derivatives onto source parameters."""
        rate = self.sidereal_rate
        # d/d(ra) = d/d(phi); d/d(dec) = -d/d(theta); d/d(t_c) picks up
        # -rate * d/d(phi) because phi = ra - gmst(t_c).
        first = {
            "ra": "_phi",
            "dec": "_theta",
            "psi": "_psi",
            "geocent_time": "_phi",
        }
        sign = {"ra": 1.0, "dec": -1.0, "psi": 1.0, "geocent_time": -rate}
        second = {
            ("ra", "ra"): ("_phiphi", 1.0),
            ("dec", "ra"): ("_phitheta", -1.0),
            ("dec", "dec"): ("_thetatheta", 1.0),
            ("psi", "ra"): ("_phipsi", 1.0),
            ("dec", "psi"): ("_thetapsi", -1.0),
            ("psi", "psi"): ("_psipsi", 1.0),
            ("geocent_time", "ra"): ("_phiphi", -rate),
            ("dec", "geocent_time"): ("_phitheta", rate),
            ("geocent_time", "psi"): ("_phipsi", -rate),
            ("geocent_time", "geocent_time"): ("_phiphi", rate**2),
        }

        self.plus, self.cross, self.delay = {}, {}, {}
        for target, source in ((self.plus, plus), (self.cross, cross), (self.delay, delay)):
            target["value"] = source[""]
            for name, label in first.items():
                target[name] = sign[name] * source[label]
            for pair, (label, factor) in second.items():
                target[tuple(sorted(pair))] = factor * source[label]

        # t_c enters tau explicitly as well: tau = (t_c - start) + vertex_delay.
        self.delay["geocent_time"] = self.delay["geocent_time"] + 1.0


def _arm_transfer(basis, arm, arm_length_in_wavelengths):
    r"""Finite-size transfer function for one arm, with angular derivatives.

    The arm enters only through :math:`y = -\hat\Omega\cdot\hat a`, so the chain
    rule is one-dimensional: :math:`\partial_a D = D_y\,\partial_a y` and
    :math:`\partial_{ab} D = D_{yy}\,\partial_a y\,\partial_b y +
    D_y\,\partial_{ab} y`.

    Parameters
    ----------
    basis : dict
        Output of :func:`_basis_derivatives`.
    arm : numpy.ndarray
        Unit vector along the arm.
    arm_length_in_wavelengths : numpy.ndarray
        :math:`fL/c`.

    Returns
    -------
    dict
        Keyed by the labels in ``_LABELS``.
    """
    projection = {
        key: -(arm @ basis["omega" + key]) if ("omega" + key) in basis else None
        for key in _LABELS
    }
    # psi does not rotate the propagation direction, so every psi derivative of
    # the projection vanishes.
    for key in ("_psi", "_phipsi", "_thetapsi", "_psipsi"):
        projection[key] = numpy.zeros_like(projection[""])

    value, first, second = _transfer_derivatives(arm_length_in_wavelengths, projection[""])

    transfer = {"": value}
    for key in ("_phi", "_theta", "_psi"):
        transfer[key] = first * projection[key]
    for key, left, right in (
        ("_phiphi", "_phi", "_phi"),
        ("_phitheta", "_phi", "_theta"),
        ("_thetatheta", "_theta", "_theta"),
        ("_phipsi", "_phi", "_psi"),
        ("_thetapsi", "_theta", "_psi"),
        ("_psipsi", "_psi", "_psi"),
    ):
        transfer[key] = (
            second * projection[left] * projection[right] + first * projection[key]
        )
    return transfer


def _unit_derivatives(shape):
    """A derivative dictionary for the constant function 1."""
    ones = numpy.ones(shape)
    zeros = numpy.zeros(shape)
    return {key: (ones if key == "" else zeros) for key in _LABELS}
