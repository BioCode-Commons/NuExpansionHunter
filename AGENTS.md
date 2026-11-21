# Publishing benchmarks

Update `benchmarks/published.json` when changing the benchmark table on the BioCode Commons website. Put the measurements and methods in [docs/08_Benchmark.md](docs/08_Benchmark.md), commit the report, and use a full commit SHA in the data file's `evidence_url`. Changes to the report or README alone do not publish the website table.

Keep the baseline and candidate revision, inputs, hardware, thread allocation, timing and memory scope, repetitions, output checks, and excluded runs explicit. Use seconds in the data file. Preserve the memory metric and units used in the report. Keep index construction separate when reporting per-sample runs. Do not estimate missing measurements.

The **Published benchmark** workflow validates the data on a pull request. After a merge to `main`, it triggers [the website publisher](https://github.com/BioCode-Commons/BioCode-Commons.github.io/actions/workflows/benchmarks.yml). The website calculates speedups from the underlying values. The organization manages the GitHub App and its credentials.
