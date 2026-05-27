# Trajectory Schema v1

Design trajectories store sufficient statistics by replicate and look. Bayes
factors are derived later by analysis cases. Production corpora store
trajectories as chunk files under `chunks/chunk-*.rds`; a finalized design adds
`chunk-index.rds`, `chunk-index.csv`, and `sha256.txt`. A monolithic
`trajectories.rds` file is optional and intended only for local/debug exports.

Common columns:

- `design_case_id`
- `replicate_id`
- `chunk_id`
- `look`
- `n`

z-family columns:

- `se`
- `true_effect`
- `estimate`

t-family columns:

- `n1`
- `n2`
- `true_effect`
- `estimate`
- `t`

binomial-family columns:

- `true_p`
- `x`

The key `(replicate_id, look)` must be unique within a design case. The true
effect or true probability must be constant across looks within a replicate.
Each chunk must contain the full look grid for exactly the replicate range
declared in the design case metadata.
