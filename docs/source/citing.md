# Citing GWForge

If GWForge contributes to a publication, please cite the paper that describes it:

Chandra, K. (2025), *gwforge: a user-friendly package to generate gravitational-wave mock data*,
Class. Quant. Grav. 42, 025003, [doi:10.1088/1361-6382/ad9b68](https://doi.org/10.1088/1361-6382/ad9b68),
[arXiv:2407.21109](https://arxiv.org/abs/2407.21109).

```bibtex
@article{Chandra:2024dhf,
    author = "Chandra, Koustav",
    title = "{gwforge: a user-friendly package to generate gravitational-wave mock data}",
    eprint = "2407.21109",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    doi = "10.1088/1361-6382/ad9b68",
    journal = "Class. Quant. Grav.",
    volume = "42",
    number = "2",
    pages = "025003",
    year = "2025"
}
```

The repository also ships a `CITATION.cff`, so GitHub's "Cite this repository" button gives the same entry.

## Acknowledging the software GWForge builds on

GWForge is a thin layer over the community's gravitational-wave software. Please also cite the packages that did the work in your analysis:

- [bilby](https://bilby-dev.github.io/bilby/citing-bilby.html) -- priors, waveform generation, interferometers and PSDs, and the sampler behind the DALI posteriors.
- [PyCBC](https://pycbc.org/pycbc/latest/html/credit.html) -- time-domain injections, PSD estimation and FFTs.
- [GWpy](https://gwpy.github.io/docs/) -- reading and writing the strain files.
- [Astropy](https://www.astropy.org/acknowledging.html) -- cosmology and units.

GWForge is released under the [MIT License](https://opensource.org/license/mit/).
