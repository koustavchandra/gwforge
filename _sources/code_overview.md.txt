# Code overview

Below is a flowchart showing the layout of the packages and modules within GWForge.

```{mermaid}
flowchart TB
    pkg(GWForge)
    population(.population)
    pkg --> population
    population_0(mass)
    population --> population_0
    population_1(_smoothed_mass)
    population --> population_1
    population_2(pairing)
    population --> population_2
    population_3(pdb_external)
    population --> population_3
    population_4(pdb_mass_sampler)
    population --> population_4
    population_5(spin)
    population --> population_5
    population_6(redshift)
    population --> population_6
    population_7(extrinsic)
    population --> population_7
    population_8(user_defined)
    population --> population_8
    population_0 ~~~ population_1
    population_1 ~~~ population_2
    population_2 ~~~ population_3
    population_3 ~~~ population_4
    population_4 ~~~ population_5
    population_5 ~~~ population_6
    population_6 ~~~ population_7
    population_7 ~~~ population_8
    ifo(.ifo)
    pkg --> ifo
    ifo_0(detectors)
    ifo --> ifo_0
    ifo_1(antenna)
    ifo --> ifo_1
    ifo_2(reproducible_noise)
    ifo --> ifo_2
    ifo_0 ~~~ ifo_1
    ifo_1 ~~~ ifo_2
    inject(.inject)
    pkg --> inject
    inject_0(bilby_inject)
    inject --> inject_0
    inject_1(pycbc_inject)
    inject --> inject_1
    inject_0 ~~~ inject_1
    glitch(.glitch)
    pkg --> glitch
    glitch_0(sine_gaussian)
    glitch --> glitch_0
    fisher(.fisher)
    pkg --> fisher
    fisher_0(parameters)
    fisher --> fisher_0
    fisher_1(stepsize)
    fisher --> fisher_1
    fisher_2(inner_product)
    fisher --> fisher_2
    fisher_3(antenna_derivatives)
    fisher --> fisher_3
    fisher_4(derivatives)
    fisher --> fisher_4
    fisher_5(matrix)
    fisher --> fisher_5
    fisher_6(dali)
    fisher --> fisher_6
    fisher_7(config)
    fisher --> fisher_7
    fisher_8(plot)
    fisher --> fisher_8
    fisher_0 ~~~ fisher_1
    fisher_1 ~~~ fisher_2
    fisher_2 ~~~ fisher_3
    fisher_3 ~~~ fisher_4
    fisher_4 ~~~ fisher_5
    fisher_5 ~~~ fisher_6
    fisher_6 ~~~ fisher_7
    fisher_7 ~~~ fisher_8
    popfisher(.population_fisher)
    pkg --> popfisher
    popfisher_0(base)
    popfisher --> popfisher_0
    popfisher_1(model)
    popfisher --> popfisher_1
    popfisher_2(catalogue)
    popfisher --> popfisher_2
    popfisher_3(fit)
    popfisher --> popfisher_3
    popfisher_4(compare)
    popfisher --> popfisher_4
    popfisher_5(population_fisher_term_I)
    popfisher --> popfisher_5
    popfisher_6(config)
    popfisher --> popfisher_6
    popfisher_7(plot)
    popfisher --> popfisher_7
    popfisher_0 ~~~ popfisher_1
    popfisher_1 ~~~ popfisher_2
    popfisher_2 ~~~ popfisher_3
    popfisher_3 ~~~ popfisher_4
    popfisher_4 ~~~ popfisher_5
    popfisher_5 ~~~ popfisher_6
    popfisher_6 ~~~ popfisher_7
    classDef pkg fill:#dae8fc,stroke:#6c8ebf,color:#000
    classDef mod fill:#ffe6cc,stroke:#d79b00,color:#000
    class pkg,population,ifo,inject,glitch,fisher,popfisher pkg
    class population_0,population_1,population_2,population_3,population_4,population_5,population_6,population_7,population_8,ifo_0,ifo_1,ifo_2,inject_0,inject_1,glitch_0,fisher_0,fisher_1,fisher_2,fisher_3,fisher_4,fisher_5,fisher_6,fisher_7,fisher_8,popfisher_0,popfisher_1,popfisher_2,popfisher_3,popfisher_4,popfisher_5,popfisher_6,popfisher_7 mod
```

The scripts are described in {doc}`usage`; every module is documented in the {doc}`API reference <autoapi/GWForge/index>`.
