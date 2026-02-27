# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.0] - 2026-02-27

### Added
- Added `fit!(A::IndependentStatistic, B::IndependentStatistic)` to merge accumulated statistics in-place.
- Added `merge`/`merge!` support for `IndependentStatistic`.
- Added README Usage example checker script under `test/check_readme_examples.jl` and integrated execution in CI.
- Added richer generated docs examples for `UnivariateStatistic` and `IndependentStatistic`.

### Changed
- **Breaking**: standardized `UnivariateStatistic` constructors around the canonical order (`T`, `K`, `I`, then sample/weights when present).
- **Breaking**: standardized `IndependentStatistic` constructors to `K`-first forms (e.g. `IndependentStatistic(K, x; dims=...)`, `IndependentStatistic(K, x, w; dims=...)`).
- Refactored constructor internals to centralize raw-moment object construction via private helpers.
- Updated tests and docs to use canonical constructor APIs and one-line summary/show expectations.

### Deprecated
- Deprecated legacy `UnivariateStatistic` constructor call orders in favor of canonical signatures.
- Deprecated legacy `IndependentStatistic(x, K; dims=...)` and `IndependentStatistic(x, w, K; dims=...)` in favor of `K`-first signatures.

### Fixed
- Fixed variance warning/error path wording for `UnivariateStatistic`.
- Ensured `fit!` methods consistently return the mutated statistic object.
- Improved type stability for sliced `IndependentStatistic` fitting paths.

### CI/Docs
- Added/updated CI checks for ExplicitImports, docs doctests (including manual doctests), and README example execution.
- Updated README and docs pages to reflect constructor API changes and weighted usage examples.

## [0.2] - 2026-02-26

### Added
- Added pretty-printing support for `UnivariateStatistic` and `IndependentStatistic`.
- Added compact and plain-text display methods (`show`) with one-line `summary` behavior.
- Added broader test coverage for display behavior and edge cases.

### Changed
- Refined display formatting and summary output for statistics objects.
- Added configurable significant digits in univariate display output.
- Refactored import style to explicit `import` usage in `src/`.
- Updated CI to include `ExplicitImports` checks.

### Fixed
- Fixed weighted-moment regression in the `K = 2` path.
- Fixed error message paths in `kurtosis` and `merge!` for `UnivariateStatistic`.
- Enabled `skewness` and `kurtosis` usage for weighted `UnivariateStatistic`.
- Improved type stability/performance in `IndependentStatistic` internals.

### Dependencies
- Removed stale dependency usage (`StaticArrays`).
- Updated dependency constraints and CI workflow dependencies.

## [0.1] - 2024-12-07
- Initial 0.1 series release.

[0.3.0]: https://github.com/FerreolS/OnlineSampleStatistics.jl/compare/v0.2...v0.3.0
[0.2]: https://github.com/FerreolS/OnlineSampleStatistics.jl/compare/v0.1...v0.2
[0.1]: https://github.com/FerreolS/OnlineSampleStatistics.jl/releases/tag/v0.1
