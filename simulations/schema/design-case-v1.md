# Design Case Schema v1

A design case describes one data-generating simulation setting.

Required fields:

- `schema_version`: integer, currently `1`
- `corpus_version`: string, currently `v1`
- `design_case_id`: stable ASCII identifier
- `test_family`: one of `z`, `t`, `binomial`
- `tier`: one of `short`, `long`, `adversarial`, `regression`
- `tags`: character vector for filtering and review
- `rationale`: one-sentence reason the case exists
- `look_grid_name`: one of `short`, `long`
- `look_grid`: numeric vector of looks
- `nsim`: number of Monte Carlo replicates
- `chunk_size`: number of replicates per chunk
- `chunks`: data frame with `chunk_id`, `replicate_start`, `replicate_end`, `seed`
- `rng_kind`: RNG kind used for chunks
- `master_seed`: integer seed used to derive chunk seeds
- `design_prior`: list describing the data-generating prior
- `generation`: list describing family-specific generation settings

Continuous design priors:

- point: `list(family = "point", mean = <number>, sd = 0)`
- normal: `list(family = "normal", mean = <number>, sd = <positive number>)`

Binomial design priors:

- point: `list(family = "point", prob = <number>)`
- beta: `list(family = "beta", shape1 = <number>, shape2 = <number>, lower = 0, upper = 1)`
- truncated beta: same as beta with `0 <= lower < upper <= 1`
