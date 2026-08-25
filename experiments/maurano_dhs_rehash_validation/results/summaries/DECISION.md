# Decision report

**Classification: exact but estimator-sensitive.**

All prespecified exact, p=18 seed, and precision-follow-up gates passed. The historical `k=10,w=30` cell is an exact-feature leader. The exploratory seed extension is complete (372/372 new runs; 101 total seeds per precision). The missing-precision interpolation is complete (909/909 runs).

## Historical-cell precision evidence

| p | seeds | ARI median [IQR] | ARI range | NMI median [IQR] | partitions / ranked hierarchies / topologies | median clade distance | exact partition (95% CI) | median exact MAE |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 12 | 101 | 0.766154 [0.766154, 0.910180] | 0.590924–0.910180 | 0.911405 [0.911405, 0.960988] | 9 / 101 / 90 | 0.222 | 44/101 (0.343–0.533) | 0.006067553 |
| 13 | 101 | 0.910180 [0.766154, 0.910180] | 0.693152–0.910180 | 0.960988 [0.911405, 0.960988] | 5 / 101 / 60 | 0.167 | 52/101 (0.419–0.610) | 0.0041579418 |
| 14 | 101 | 0.910180 [0.766154, 0.910180] | 0.766154–0.910180 | 0.960988 [0.911405, 0.960988] | 3 / 98 / 32 | 0.056 | 56/101 (0.457–0.648) | 0.0028850169 |
| 15 | 101 | 0.766154 [0.766154, 0.910180] | 0.766154–0.910180 | 0.911405 [0.911405, 0.960988] | 2 / 82 / 17 | 0.056 | 50/101 (0.400–0.591) | 0.0020443047 |
| 16 | 101 | 0.910180 [0.766154, 0.910180] | 0.766154–0.910180 | 0.960988 [0.911405, 0.960988] | 2 / 36 / 3 | 0.056 | 54/101 (0.438–0.629) | 0.0013264078 |
| 17 | 101 | 0.910180 [0.766154, 0.910180] | 0.766154–0.910180 | 0.960988 [0.911405, 0.960988] | 2 / 18 / 4 | 0.000 | 61/101 (0.506–0.694) | 0.00080842638 |
| 18 | 101 | 0.910180 [0.766154, 0.910180] | 0.766154–0.910180 | 0.960988 [0.911405, 0.960988] | 2 / 7 / 2 | 0.000 | 52/101 (0.419–0.610) | 0.00053106263 |
| 19 | 101 | 0.910180 [0.766154, 0.910180] | 0.766154–0.910180 | 0.960988 [0.911405, 0.960988] | 2 / 4 / 2 | 0.000 | 52/101 (0.419–0.610) | 0.00035367915 |
| 20 | 101 | 0.910180 [0.766154, 0.910180] | 0.766154–0.910180 | 0.960988 [0.911405, 0.960988] | 2 / 4 / 2 | 0.000 | 61/101 (0.506–0.694) | 0.00024467689 |
| 21 | 101 | 0.910180 [0.766154, 0.910180] | 0.766154–0.910180 | 0.960988 [0.911405, 0.960988] | 2 / 4 / 2 | 0.000 | 58/101 (0.477–0.666) | 0.0001747532 |
| 22 | 101 | 0.910180 [0.766154, 0.910180] | 0.766154–0.910180 | 0.960988 [0.911405, 0.960988] | 2 / 4 / 2 | 0.000 | 54/101 (0.438–0.629) | 0.00011989787 |
| 23 | 101 | 0.766154 [0.766154, 0.910180] | 0.766154–0.910180 | 0.911405 [0.911405, 0.960988] | 2 / 3 / 2 | 0.000 | 41/101 (0.315–0.503) | 8.5429165e-05 |
| 24 | 101 | 0.766154 [0.766154, 0.910180] | 0.766154–0.910180 | 0.911405 [0.911405, 0.960988] | 2 / 2 / 1 | 0.000 | 44/101 (0.343–0.533) | 6.1937514e-05 |

## Precision dispersion

Exact-partition recovery differs across paired precisions by Cochran's Q (Q=21.788, df=12, p=0.03997), but the recovery frequency is nonmonotonic and no adjacent exact McNemar comparison passes Holm correction. Coarse ARI/NMI states therefore do not track the strong, monotonic improvement in matrix error and cophenetic agreement.

## Eight-class muscle-merged sensitivity

| p | median eight-class ARI | ARI range | partitions | exact eight-cluster partition |
|---:|---:|---:|---:|---:|
| 12 | 0.861952 | 0.498818–1.000000 | 11 | 48/101 |
| 13 | 0.861952 | 0.498818–1.000000 | 8 | 60/101 |
| 14 | 0.861952 | 0.733434–0.904642 | 4 | 64/101 |
| 15 | 0.861952 | 0.733434–0.861952 | 2 | 75/101 |
| 16 | 0.861952 | 0.733434–0.861952 | 2 | 86/101 |
| 17 | 0.861952 | 0.733434–0.861952 | 2 | 95/101 |
| 18 | 0.861952 | 0.861952–0.861952 | 1 | 101/101 |
| 19 | 0.861952 | 0.861952–0.861952 | 1 | 101/101 |
| 20 | 0.861952 | 0.861952–0.861952 | 1 | 101/101 |
| 21 | 0.861952 | 0.861952–0.861952 | 1 | 101/101 |
| 22 | 0.861952 | 0.861952–0.861952 | 1 | 101/101 |
| 23 | 0.861952 | 0.861952–0.861952 | 1 | 101/101 |
| 24 | 0.861952 | 0.861952–0.861952 | 1 | 101/101 |

The eight-cluster cut is stable earlier: all 101 seeds reproduce the exact eight-cluster partition and identical ARI/NMI at every precision from p=18 through p=24. Thus the persistent high-precision ten-class split is specific to merge ranking at its exceptionally narrow cut, not a failure of the coarser biological grouping or unranked hierarchical convergence. Cochran's Q for the eight-cluster recovery series is 393.500 (df=12, p=8.997e-77); no adjacent exact McNemar comparison passes Holm correction.

Tissue rankings are exploratory because selection and recovery use the same 20 labelled samples. No compatibility default or manuscript text was changed.
