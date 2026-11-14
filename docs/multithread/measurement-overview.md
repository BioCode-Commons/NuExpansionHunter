# Frozen v6: complete selected measurement matrix

14 fixture/settings cases (fixed and automatic I/O share the genotyped10k fixture); 27 main configurations; 81 selected formal pairs plus 27 warmup pairs. Including the preserved original T8 and its supplement: 84 formal pairs plus 28 warmup pairs, 224 analysis processes. Cache preparation processes are excluded from these counts.

The original genotyped10k-fixed T8 configuration has an objective timing inconsistency in baseline formal round 2: process wall and application/UTC elapsed times disagree. Its entire original result file and all runs remain in evidence. Main tables use the original T1/T4 configurations and the complete independent T8 recheck (three formal pairs plus one warmup), rather than selecting individual favorable repetitions. No measurements are changed. Original T8 is excluded only from selected main statistics.

Read-phase boundaries differ: b00650d logs the phase before unmatched-read and input-stream destruction; v6 logs it after the input scope closes. These diagnostic elapsed times are not strictly equivalent intervals. Pending-mate cleanup may contribute to an observed difference; no entire delta is attributed to cleanup. End-to-end process wall time is the primary comparison.

B = frozen RAMOPT baseline; C = frozen v6. Medians/maxima use only formal runs. Paired time change is the median C/B ratio minus one. RSS maxima change compares the two formal maxima. Largest paired RSS delta reports separate maxima of absolute and percentage differences over matched pairs; these two maxima can come from different pairs. Negative values are retained. Memory reference means the original baseline peak RSS plus 20,000,000 bytes per compute worker. The user later allowed a modest disclosed excess; this reference is not an acceptance gate. Every observed slowdown and memory-reference exceedance remains visible. Complete outputs do not establish performance acceptance. Three pairs on one ARM host are not a server-wide claim.

| Case/settings | T | E2E median B/C (s) | Paired E2E change | Read median B/C (s) | Paired read change | Peak RSS max B/C (MiB) | RSS maxima change | Largest paired RSS delta | Memory reference |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| primary100k | 8 | 184.026 / 90.113 | -51.03% | 169.532 / 79.725 | -53.30% | 3839.406 / 4142.922 | +7.91% | +381.703 MiB / +11.43% | Exceeded |
| genotyped10k-fixed | 1 | 14.269 / 13.013 | -8.80% | 12.056 / 10.730 | -11.00% | 433.172 / 435.469 | +0.53% | +3.984 MiB / +0.92% | Within |
| genotyped10k-fixed | 4 | 11.879 / 9.112 | -22.22% | 10.699 / 7.927 | -23.49% | 435.984 / 445.500 | +2.18% | +9.516 MiB / +2.18% | Within |
| genotyped10k-auto-io | 8 | 21.835 / 10.900 | -51.23% | 20.458 / 9.418 | -55.05% | 442.750 / 449.516 | +1.53% | +6.766 MiB / +1.53% | Within |
| repository-example | 1 | 0.150 / 0.111 | -26.35% | 0.098 / 0.059 | -40.03% | 16.344 / 17.984 | +10.04% | +1.641 MiB / +10.08% | Within |
| repository-example | 8 | 0.151 / 0.074 | -51.07% | 0.099 / 0.021 | -78.33% | 16.438 / 23.000 | +39.92% | +6.672 MiB / +40.86% | Within |
| longread1000 | 1 | 0.791 / 0.735 | +0.61% | 0.682 / 0.625 | -0.62% | 24.328 / 26.031 | +7.00% | +1.781 MiB / +7.35% | Within |
| longread1000 | 8 | 0.910 / 0.454 | -48.98% | 0.334 / 0.327 | -1.81% | 61.641 / 44.641 | -27.58% | -16.969 MiB / -27.58% | Within |
| higherror-indels150 | 1 | 2.145 / 2.124 | +0.14% | 2.036 / 2.013 | +0.23% | 23.891 / 25.000 | +4.64% | +1.219 MiB / +5.14% | Within |
| higherror-indels150 | 8 | 0.851 / 0.913 | +7.20% | 0.775 / 0.825 | +6.63% | 29.641 / 33.344 | +12.49% | +3.953 MiB / +13.46% | Within |
| iupac-missingquality150 | 1 | 0.632 / 0.580 | -8.80% | 0.513 / 0.467 | -9.17% | 24.406 / 25.266 | +3.52% | +0.984 MiB / +4.07% | Within |
| iupac-missingquality150 | 8 | 0.281 / 0.256 | -9.94% | 0.195 / 0.173 | -11.25% | 28.141 / 31.328 | +11.33% | +3.562 MiB / +12.83% | Within |
| expandedrepeat150-genotyped | 1 | 5.195 / 4.809 | -7.42% | 5.007 / 4.592 | -8.28% | 46.125 / 49.469 | +7.25% | +3.578 MiB / +7.80% | Within |
| expandedrepeat150-genotyped | 8 | 4.496 / 2.253 | -49.01% | 4.409 / 2.150 | -50.36% | 50.938 / 66.734 | +31.01% | +15.797 MiB / +31.01% | Within |
| dense-fanout-large | 1 | 10.882 / 8.877 | -19.27% | 10.010 / 8.105 | -20.20% | 136.609 / 142.078 | +4.00% | +5.469 MiB / +4.00% | Within |
| dense-fanout-large | 8 | 10.766 / 6.296 | -42.12% | 10.512 / 5.999 | -43.49% | 167.203 / 185.078 | +10.69% | +20.406 MiB / +12.39% | Within |
| hot-locus | 1 | 0.480 / 0.325 | -33.00% | 0.386 / 0.228 | -41.68% | 35.328 / 32.062 | -9.24% | -3.234 MiB / -9.16% | Within |
| hot-locus | 8 | 0.455 / 0.237 | -50.40% | 0.374 / 0.157 | -63.42% | 36.391 / 43.016 | +18.21% | +7.250 MiB / +20.45% | Within |
| empty-input | 1 | 0.051 / 0.050 | -2.59% | 0.000105 / 0.000179 | +38.52% | 15.859 / 15.781 | -0.49% | -0.031 MiB / -0.20% | Within |
| empty-input | 8 | 0.041 / 0.039 | -3.86% | 0.000158 / 0.000171 | +9.87% | 15.922 / 15.844 | -0.49% | -0.031 MiB / -0.20% | Within |
| pending-million-missing-mates | 1 | 0.622 / 0.616 | +0.83% | 0.444 / 0.550 | +25.07% | 858.609 / 795.219 | -7.38% | -63.359 MiB / -7.38% | Within |
| pending-million-missing-mates | 8 | 0.595 / 0.572 | -2.41% | 0.438 / 0.521 | +20.02% | 858.922 / 795.422 | -7.39% | -63.344 MiB / -7.38% | Within |
| pending-large-far-apart | 1 | 84.627 / 85.454 | +0.98% | 82.568 / 82.913 | +0.42% | 117.078 / 116.484 | -0.51% | +15.656 MiB / +15.53% | Within |
| pending-large-far-apart | 8 | 34.038 / 36.142 | +6.18% | 33.118 / 35.077 | +5.92% | 211.859 / 237.422 | +12.07% | +25.562 MiB / +12.07% | Within |
| mixed100k | 1 | 17.658 / 16.652 | -5.70% | 1.936 / 1.663 | -14.08% | 2381.750 / 2384.062 | +0.10% | +2.312 MiB / +0.10% | Within |
| mixed100k | 8 | 10.998 / 9.844 | -10.88% | 1.992 / 0.901 | -55.05% | 2381.906 / 2397.094 | +0.64% | +19.484 MiB / +0.82% | Within |
| genotyped10k-t8-recheck | 8 | 21.084 / 9.772 | -53.65% | 19.717 / 8.559 | -56.59% | 441.812 / 449.250 | +1.68% | +9.250 MiB / +2.10% | Within |
