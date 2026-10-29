A lightweight code to generate mock gravitational wave detector data using user-defined population and arbitrary detector sensitivity. It ships six tools -- `gwforge_population`, `gwforge_noise`, `gwforge_inject`, `gwforge_glitch`, `gwforge_fisher` and `gwforge_population_fisher` -- plus `gwforge_workflow`, which chains noise, injections and glitches into an HTCondor DAG, and `gwforge_optimal_snr` and `gwforge_estimate_psd` for the bookkeeping in between.

# Installation
To install this library for development or use, follow the instructions [here](https://koustavchandra.github.io/gwforge/install.html)

# Documentation

Heavily under-development documentation is accessible [here](https://koustavchandra.github.io/gwforge/index.html)

# Citation
If you use `gwforge`, please cite:

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

# License
`gwforge` is released under the MIT License. (see [LICENSE](LICENSE)).
