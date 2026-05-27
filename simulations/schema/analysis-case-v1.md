# Analysis Case Schema v1

An analysis case applies an analysis prior and BF thresholds to one design case.

Required fields:

- `schema_version`: integer, currently `1`
- `corpus_version`: string, currently `v1`
- `analysis_case_id`: stable ASCII identifier
- `design_case_id`: identifier of the source design case
- `test_family`: one of `z`, `t`, `binomial`
- `tier`: one of `short`, `long`, `adversarial`, `regression`
- `tags`: character vector for filtering and review
- `rationale`: one-sentence reason the case exists
- `bf_type`: BF type for the family
- `k1`: BF01 threshold for stopping in favor of H1
- `k0`: BF01 threshold for stopping in favor of H0
- `analysis_prior`: list of family-specific prior parameters
- `strict`: logical for sequential integration where applicable
- `drange`: root-search range for t-test functions where applicable

Supported `bf_type` values:

- z family: `normal`, `directional`, `moment`
- t family: `t`
- binomial family: `point`, `direction`
