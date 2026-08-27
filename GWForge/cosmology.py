import logging
from lal import C_SI

import numpy
from astropy import units
from scipy.special import roots_legendre

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)
# Speed of light in km/s, so that ``c / H0`` is Mpc for ``H0`` in km/s/Mpc.
SPEED_OF_LIGHT = C_SI / 1000.0

# Gauss-Legendre order for the comoving-distance integral. 1/E is smooth on
# [0, z], so convergence is geometric: measured against ``scipy.integrate.quad``
# this reaches ~1e-14 relative by order 32 and is flat thereafter. 64 is used
# because the cost is one matrix product either way.
QUADRATURE_ORDER = 64

# Newton iterations for :meth:`FlatwCDM.redshift_of_distance`. The bracketing
# table gives ~1e-3 and Newton doubles the digits each step, so 4 is already
# past double precision; the residual is checked rather than assumed.
NEWTON_ITERATIONS = 6

# Redshift nodes for the initial guess in the distance inversion. Log-spaced
# because d_L(z) is close to linear in z at low z and the low-z end is where
# the events are.
_GUESS_NODES = numpy.concatenate(
    [[0.0], numpy.logspace(-4, numpy.log10(100.0), 512)]
)

# Parameters :meth:`FlatwCDM.derivatives` knows how to differentiate.
COSMOLOGY_PARAMETERS = ("H0", "Om0", "w0")


def astropy_cosmology(
    name="Planck18", H0=None, Om0=None, Ode0=None, Tcmb0=2.725, Ob0=None
):
    """Resolve a GWForge cosmology specification to an ``astropy`` cosmology.

    Setting ``Tcmb0 = 0`` removes radiation and massive neutrinos, which makes
    the result *identical* to ``FlatwCDM(H0, Om0, w0=-1)`` when
    ``Ode0 = 1 - Om0``.

    Parameters
    ----------
    name : str
        An ``astropy.cosmology`` realization name, or anything else (e.g.
        ``"custom"``) to request an explicit ``LambdaCDM``.
    H0 : float or None
        Hubble constant in km/s/Mpc. Required if ``name`` is not a realization.
    Om0, Ode0 : float or None
        Present-day matter and dark-energy density parameters.
    Tcmb0 : float
        CMB temperature in K. ``0`` disables radiation and neutrinos.
    Ob0 : float or None
        Baryon density parameter.

    Returns
    -------
    astropy.cosmology.FLRW

    Raises
    ------
    ValueError
        If ``name`` is not a realization and the explicit parameters are
        incomplete.
    """
    from astropy import cosmology as astropy_module

    realization = getattr(astropy_module, str(name), None)
    if isinstance(realization, astropy_module.FLRW):
        return realization
    if H0 is None or Om0 is None or Ode0 is None:
        raise ValueError(
            "Cosmology '{}' is not an astropy realization, and building a "
            "custom LambdaCDM needs H0, Om0 and Ode0 (got {}, {}, {}).".format(
                name, H0, Om0, Ode0
            )
        )
    logging.info("Building a custom LambdaCDM cosmology for '%s'", name)
    settings = {"H0": H0, "Om0": Om0, "Ode0": Ode0, "Tcmb0": Tcmb0}
    # astropy 8 rejects an explicit Ob0=None
    if Ob0 is not None:
        settings["Ob0"] = Ob0
    return astropy_module.LambdaCDM(**settings)


class FlatwCDM:
    r"""Flat dark-energy background with analytic parameter derivatives.

    .. math::

       E(z)^2 = \Omega_{m,0}(1 + z)^3
                + (1 - \Omega_{m,0})(1 + z)^{3(1 + w_0)}

    Usage
    -----
    >>> cosmology = FlatwCDM(H0=67.66, Om0=0.3111)
    >>> cosmology.luminosity_distance(1.0)            # Mpc
    >>> cosmology.derivatives(1.0)["Om0"]["luminosity_distance"]

    Attributes
    ----------
    H0 : float
        Hubble constant, km/s/Mpc.
    Om0 : float
        Present-day matter density parameter.
    w0 : float
        Dark-energy equation of state. ``-1`` is flat LCDM.
    hubble_distance : float
        :math:`D_H = c / H_0`, Mpc.
    """

    def __init__(self, H0=67.66, Om0=0.3111, w0=-1.0):
        """
        Parameters
        ----------
        H0 : float
            Hubble constant in km/s/Mpc.
        Om0 : float
            Present-day matter density parameter, in ``(0, 1)``.
        w0 : float
            Dark-energy equation of state.
        """
        if not 0.0 < Om0 < 1.0:
            raise ValueError("Om0 must lie in (0, 1), got {}".format(Om0))
        if H0 <= 0.0:
            raise ValueError("H0 must be positive, got {}".format(H0))
        self.H0 = float(H0)
        self.Om0 = float(Om0)
        self.w0 = float(w0)
        self.hubble_distance = SPEED_OF_LIGHT / self.H0
        self._nodes, self._weights = roots_legendre(QUADRATURE_ORDER)

    def __repr__(self):
        return "FlatwCDM(H0={}, Om0={}, w0={})".format(self.H0, self.Om0, self.w0)

    @classmethod
    def from_astropy(cls, cosmology, w0=-1.0, fold_neutrinos=True):
        """Build a :class:`FlatwCDM` matching an ``astropy`` cosmology.

        ``astropy`` realizations carry radiation and massive neutrinos, which
        this two-component background does not have. ``fold_neutrinos`` adds
        ``Onu0`` to ``Om0``, treating neutrinos as cold matter -- a good
        approximation at the redshifts gravitational-wave catalogues reach, and
        the convention this module documents its residuals against.

        Parameters
        ----------
        cosmology : astropy.cosmology.FLRW
        w0 : float
            Equation of state to give the result.
        fold_neutrinos : bool
            Add ``Onu0`` to ``Om0``.

        Returns
        -------
        FlatwCDM
        """
        matter = float(cosmology.Om0)
        if fold_neutrinos:
            matter += float(cosmology.Onu0)
        return cls(H0=float(cosmology.H0.value), Om0=matter, w0=w0)

    def to_astropy(self):
        """The exactly equivalent ``astropy`` cosmology (no radiation, no neutrinos).

        Returns
        -------
        astropy.cosmology.FlatLambdaCDM or astropy.cosmology.FlatwCDM
        """
        from astropy import cosmology as astropy_module

        if self.w0 == -1.0:
            return astropy_module.FlatLambdaCDM(H0=self.H0, Om0=self.Om0, Tcmb0=0.0)
        return astropy_module.FlatwCDM(
            H0=self.H0, Om0=self.Om0, w0=self.w0, Tcmb0=0.0
        )

    # -- background ------------------------------------------------------

    def efunc(self, redshift):
        r""":math:`E(z) = H(z) / H_0`."""
        return numpy.sqrt(self._efunc_squared(numpy.asarray(redshift, dtype=float)))

    def _efunc_squared(self, redshift):
        one_plus_z = 1.0 + redshift
        return self.Om0 * one_plus_z**3 + (1.0 - self.Om0) * one_plus_z ** (
            3.0 * (1.0 + self.w0)
        )

    def defunc_dz(self, redshift):
        r""":math:`dE/dz`.

        Carries the dark-energy term. The LCDM expression
        :math:`1.5\,\Omega_{m,0}(1+z)^2 / E` is the :math:`w_0 = -1` special
        case only, and using it under wCDM is wrong by 26% at :math:`w_0 = -0.8`.
        """
        redshift = numpy.asarray(redshift, dtype=float)
        one_plus_z = 1.0 + redshift
        matter = 3.0 * self.Om0 * one_plus_z**2
        dark_energy = (
            3.0
            * (1.0 + self.w0)
            * (1.0 - self.Om0)
            * one_plus_z ** (3.0 * (1.0 + self.w0) - 1.0)
        )
        return (matter + dark_energy) / (2.0 * self.efunc(redshift))

    def _quadrature_nodes(self, redshift):
        """Gauss-Legendre nodes and half-width for ``[0, z]``, broadcast over ``z``.

        Returns ``(sample_redshifts, half)`` with shapes ``(..., order)`` and
        ``(..., 1)``, so the integral of ``f`` is ``(half * (weights * f).sum(-1))``.
        """
        redshift = numpy.asarray(redshift, dtype=float)
        half = redshift[..., numpy.newaxis] / 2.0
        return half * (self._nodes + 1.0), half

    def comoving_distance(self, redshift):
        r""":math:`D_C(z) = D_H \int_0^z dz'/E(z')`, Mpc."""
        sample, half = self._quadrature_nodes(redshift)
        integral = half[..., 0] * numpy.sum(
            self._weights / self.efunc(sample), axis=-1
        )
        return self.hubble_distance * integral

    def luminosity_distance(self, redshift):
        r""":math:`d_L = (1 + z) D_C(z)`, Mpc."""
        redshift = numpy.asarray(redshift, dtype=float)
        return (1.0 + redshift) * self.comoving_distance(redshift)

    def ddL_dz(self, redshift):
        r""":math:`dd_L/dz = D_C(z) + (1 + z) D_H / E(z)`, Mpc."""
        redshift = numpy.asarray(redshift, dtype=float)
        return self.comoving_distance(redshift) + (
            1.0 + redshift
        ) * self.hubble_distance / self.efunc(redshift)

    def d2dL_dz2(self, redshift):
        r""":math:`d^2 d_L/dz^2 = 2 D_H/E - (1 + z) D_H E'/E^2`, Mpc."""
        redshift = numpy.asarray(redshift, dtype=float)
        efunc = self.efunc(redshift)
        return (
            2.0 * self.hubble_distance / efunc
            - (1.0 + redshift)
            * self.hubble_distance
            * self.defunc_dz(redshift)
            / efunc**2
        )

    def differential_comoving_volume(self, redshift):
        r"""Full-sky :math:`dV_c/dz = 4\pi D_C^2 D_H / E`, Mpc^3.

        Note the :math:`4\pi`: ``astropy``'s method of the same name is *per
        steradian*.
        """
        redshift = numpy.asarray(redshift, dtype=float)
        comoving = self.comoving_distance(redshift)
        return (
            4.0
            * numpy.pi
            * comoving**2
            * self.hubble_distance
            / self.efunc(redshift)
        )

    def redshift_of_distance(self, distance):
        r"""Invert :math:`d_L(z)`, by table-seeded Newton on the exact derivative.

        Parameters
        ----------
        distance : array_like
            Luminosity distance in Mpc.

        Returns
        -------
        numpy.ndarray
        """
        distance = numpy.asarray(distance, dtype=float)
        table = self.luminosity_distance(_GUESS_NODES)
        redshift = numpy.interp(distance, table, _GUESS_NODES)
        for _ in range(NEWTON_ITERATIONS):
            residual = self.luminosity_distance(redshift) - distance
            step = residual / self.ddL_dz(redshift)
            redshift = numpy.maximum(redshift - step, 0.0)
        return redshift

    # -- derivatives -----------------------------------------------------

    def _dinv_efunc(self, redshift):
        """``d(1/E)/dLambda`` for each cosmology parameter, as a dict of arrays.

        These are the closed-form integrands that make
        :meth:`derivatives` exact: the same Gauss-Legendre weights that build
        :math:`D_C` also build :math:`\\partial D_C/\\partial\\Lambda`.
        """
        one_plus_z = 1.0 + redshift
        efunc = self.efunc(redshift)
        dark_energy_power = one_plus_z ** (3.0 * (1.0 + self.w0))
        prefactor = -0.5 / efunc**3
        return {
            # E carries no H0; the H0 dependence is entirely in D_H.
            "H0": numpy.zeros_like(efunc),
            "Om0": prefactor * (one_plus_z**3 - dark_energy_power),
            "w0": prefactor
            * 3.0
            * (1.0 - self.Om0)
            * dark_energy_power
            * numpy.log(one_plus_z),
        }

    def _defunc(self, redshift):
        """``dE/dLambda`` at fixed ``z``, as a dict of arrays."""
        one_plus_z = 1.0 + redshift
        efunc = self.efunc(redshift)
        dark_energy_power = one_plus_z ** (3.0 * (1.0 + self.w0))
        return {
            "H0": numpy.zeros_like(efunc),
            "Om0": (one_plus_z**3 - dark_energy_power) / (2.0 * efunc),
            "w0": 3.0
            * (1.0 - self.Om0)
            * dark_energy_power
            * numpy.log(one_plus_z)
            / (2.0 * efunc),
        }

    def derivatives(self, redshift, parameters=COSMOLOGY_PARAMETERS):
        r"""Analytic derivatives of the background quantities at fixed ``z``.

        Parameters
        ----------
        redshift : array_like
        parameters : sequence of str
            Any subset of :data:`COSMOLOGY_PARAMETERS`.

        Returns
        -------
        dict
            ``{parameter: {"comoving_distance": ..., "luminosity_distance": ...,
            "ddL_dz": ..., "differential_comoving_volume": ...}}``, every entry
            an array shaped like ``redshift``.

        Notes
        -----
        ``H0`` enters only through :math:`D_H = c/H_0`, because :math:`E(z)` is
        dimensionless and carries no :math:`H_0`. Hence
        :math:`\partial_{H_0} D_C = -D_C/H_0`,
        :math:`\partial_{H_0}(dV_c/dz) = -3 (dV_c/dz)/H_0` and
        :math:`\partial_{H_0}(dd_L/dz) = -(dd_L/dz)/H_0` -- the
        :math:`-3/H_0 + 1/H_0 = -2/H_0` that the spectral-siren measure term
        reduces to.
        """
        unknown = [name for name in parameters if name not in COSMOLOGY_PARAMETERS]
        if unknown:
            raise ValueError(
                "Unknown cosmology parameter(s) {}; expected a subset of {}.".format(
                    unknown, list(COSMOLOGY_PARAMETERS)
                )
            )
        redshift = numpy.asarray(redshift, dtype=float)
        one_plus_z = 1.0 + redshift
        efunc = self.efunc(redshift)
        comoving = self.comoving_distance(redshift)
        volume = 4.0 * numpy.pi * comoving**2 * self.hubble_distance / efunc

        sample, half = self._quadrature_nodes(redshift)
        dinv_efunc = self._dinv_efunc(sample)
        defunc = self._defunc(redshift)

        results = {}
        for name in parameters:
            if name == "H0":
                dcomoving = -comoving / self.H0
                dhubble_over_efunc = -self.hubble_distance / (efunc * self.H0)
                dvolume = -3.0 * volume / self.H0
            else:
                dcomoving = self.hubble_distance * half[..., 0] * numpy.sum(
                    self._weights * dinv_efunc[name], axis=-1
                )
                dhubble_over_efunc = (
                    -self.hubble_distance * defunc[name] / efunc**2
                )
                dvolume = 4.0 * numpy.pi * self.hubble_distance * (
                    2.0 * comoving * dcomoving / efunc
                    - comoving**2 * defunc[name] / efunc**2
                )
            results[name] = {
                "comoving_distance": dcomoving,
                "luminosity_distance": one_plus_z * dcomoving,
                "ddL_dz": dcomoving + one_plus_z * dhubble_over_efunc,
                "differential_comoving_volume": dvolume,
            }
        return results

    def dredshift_dparameter(self, redshift, parameters=COSMOLOGY_PARAMETERS):
        r"""``dz/dLambda`` at **fixed luminosity distance**.

        Implicit differentiation of :math:`d_L(z;\Lambda) = \mathrm{const}`:

        .. math::

           \frac{\partial z}{\partial\Lambda}
             = -\left(\frac{\partial d_L}{\partial\Lambda}\right)_{\!z}
                \Big/ \frac{d d_L}{d z}.

        This is the chain that carries cosmology into the spectral-siren score:
        a trial cosmology moves the redshift assigned to a measured distance,
        which moves the source-frame masses and the merger-rate shape.

        Parameters
        ----------
        redshift : array_like
            Redshift at the fiducial cosmology.
        parameters : sequence of str

        Returns
        -------
        dict
            ``{parameter: array}``.
        """
        derivatives = self.derivatives(redshift, parameters=parameters)
        ddL_dz = self.ddL_dz(redshift)
        return {
            name: -derivatives[name]["luminosity_distance"] / ddL_dz
            for name in parameters
        }

    def parameters(self):
        """The cosmology parameters as a plain dict."""
        return {"H0": self.H0, "Om0": self.Om0, "w0": self.w0}


def differential_comoving_volume(cosmology, redshift):
    """Full-sky ``dV_c/dz`` in Mpc^3 for either cosmology flavour.

    ``astropy``'s ``differential_comoving_volume`` is per steradian and carries
    units; :class:`FlatwCDM`'s is already full-sky and unitless. Population code
    wants one number either way, so it goes through here.

    Parameters
    ----------
    cosmology : FlatwCDM or astropy.cosmology.FLRW
    redshift : array_like

    Returns
    -------
    numpy.ndarray
    """
    if isinstance(cosmology, FlatwCDM):
        return cosmology.differential_comoving_volume(redshift)
    return (
        (cosmology.differential_comoving_volume(redshift) * 4.0 * numpy.pi * units.sr)
        .to(units.Mpc**3)
        .value
    )


def luminosity_distance(cosmology, redshift):
    """Luminosity distance in Mpc for either cosmology flavour.

    Parameters
    ----------
    cosmology : FlatwCDM or astropy.cosmology.FLRW
    redshift : array_like

    Returns
    -------
    numpy.ndarray
    """
    if isinstance(cosmology, FlatwCDM):
        return cosmology.luminosity_distance(redshift)
    return cosmology.luminosity_distance(redshift).to(units.Mpc).value
