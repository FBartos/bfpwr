# Fixture Summary Schema v1

Fixture summaries are derived artifacts built from cached BF trajectory chunks.
They do not replace the raw chunks; they provide compact probability curves and
search targets for package tests.

## Fixed summaries

`tail-summary.rds` stores one row per BF prior, design, fixed `n`, evidence
threshold, and tail. `tail = "H1"` means `BF01 <= k`, with `k = 1 / evidence`.
`tail = "H0"` means `BF01 >= k`, with `k = evidence`.

`decision-summary.rds` combines the two tails for symmetric threshold pairs and
stores `pH1`, `pH0`, and `pInc`.

`search-summary.rds` records the first grid `n` whose empirical probability is
at least the requested target probability.

## Sequential summaries

`cumulative-summary.rds` stores cumulative stopping probabilities by look for a
registered schedule and threshold pair.

`final-summary.rds` stores final stopping probabilities and sample-size
summaries for the whole schedule. Inconclusive sequences are assigned the
schedule maximum `n` when computing `EN`, `VarN`, and quantiles.

`search-summary.rds` records the first look within each schedule whose
cumulative probability reaches the requested target probability. It may be a
zero-row table with the standard search-summary columns when sequential search
targets are intentionally disabled.

## Storage

Fixture directories contain `spec.rds`, one or more summary `.rds` files,
`overview.csv`, `manifest.csv`, and `sha256.txt`. Full summary CSV files are
intentionally not written by default because long IDs can push CSV artifacts
above GitHub's single-file size limit.
