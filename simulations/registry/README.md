# Simulation Registry

The production fixture registry is defined in small family-specific R files and
rendered to CSV manifests for review.

Tiers:

- `short`: regular stored fixture cases on the short grid.
- `long`: regular stored fixture cases on the long grid.
- `adversarial`: numerically stressful or deliberately misspecified cases.
- `regression`: cases tied to known bugs or fragile behavior.

`smoke` cases are intentionally not a registry tier. Developer smoke helpers
live in `simulations/R/smoke.R`; package smoke tests should generate tiny
temporary examples directly.

Source files:

```text
registry/designs/*.R
registry/analyses/*.R
```

Review files generated from the source registry:

```text
registry/manifests/design-cases.csv
registry/manifests/analysis-cases.csv
registry/manifests/run-plan.csv
```

Render manifests from the repository root:

```sh
Rscript simulations/scripts/render_registry_manifests.R
```

Validate registry metadata:

```sh
Rscript simulations/scripts/validate_registry.R
```

## Adding A Design Case

Add cases to the relevant family file under `registry/designs/` using the
factory helpers. Include `tier`, `tags`, and `rationale` so later reviewers can
understand why the case exists.

```r
bfpwr_sim_design_z(
    id = "z-dnorm-m0p5-s0p1-short",
    tier = "short",
    look_grid = "short",
    nsim = 10000,
    chunk_size = 1000,
    seed = 20261001,
    dpm = 0.5,
    dpsd = 0.1,
    tags = c("medium-effect", "normal-design"),
    rationale = "Baseline uncertain medium-effect design for z-family analyses."
)
```

## Adding An Analysis Case

Add cases to the relevant family file under `registry/analyses/`. Analysis cases
usually loop over selected designs and call an analysis factory.

```r
lapply(designs, function(design) {
    bfpwr_sim_analysis_z_normal(
        design = design,
        pm = 0,
        psd = 1,
        k1 = 1 / 10,
        k0 = 10,
        tags = c("normal-prior", "standard-thresholds"),
        rationale = "Standard normal-prior z analysis applied to selected designs."
    )
})
```

For local dry runs, keep the registry case unchanged and use
`bfpwr_sim_design_case_override(design, nsim = 100, chunk_size = 100)` in a
temporary script. Keep production registry tiers limited to `short`, `long`,
`adversarial`, and `regression`.
