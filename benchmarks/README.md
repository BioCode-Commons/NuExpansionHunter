# Publishing benchmarks

`published.json` contains the comparisons displayed by the BioCode Commons
website. Record the results and methods in the [benchmark report](../docs/08_Benchmark.md).
Commit the report first, then pin `evidence_url` and the report link in
`caption` to the full SHA of the report commit.

Each row records the dataset, resource allocation, elapsed seconds, and peak
RSS for ExpansionHunter v5.0.0 and NuExpansionHunter. Set `rss_unit`,
`rss_scope`, and `time_scope` to match the measurements. The website computes
speedup from the timings. Use `display` to set the precision and speedup-column
visibility. Use `caption` and `note` for context in text segments with optional
HTTPS GitHub links.

The [publication workflow](../.github/workflows/published-benchmark.yml)
validates the record with the
[shared validator](https://github.com/BioCode-Commons/BioCode-Commons.github.io/blob/main/scripts/benchmarks.py).
On `main` it requests a site refresh using the organization’s GitHub App.
The organization configures `BENCHMARK_APP_ID` as an Actions variable and
`BENCHMARK_APP_PRIVATE_KEY` as an Actions secret. The site also checks hourly
and deploys after validating data from all four project repositories.
