# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[Unreleased]: https://github.com/FerreolS/OnlineSampleStatistics.jl/compare/v0.2...HEAD
[0.2]: https://github.com/FerreolS/OnlineSampleStatistics.jl/compare/v0.1...v0.2
[0.1]: https://github.com/FerreolS/OnlineSampleStatistics.jl/releases/tag/v0.1
