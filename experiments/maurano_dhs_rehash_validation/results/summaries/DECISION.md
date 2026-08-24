# Decision report

**Classification: exact but estimator-sensitive.**

All prespecified exact, p=18 seed, and precision-follow-up gates passed. The historical `k=10,w=30` cell is an exact-feature leader. The exploratory seed extension is complete (372/372 new runs; 101 total seeds per precision).

## Historical-cell precision evidence

| p | seeds | ARI min | ARI median | ARI max | NMI min | NMI median | NMI max | exact partition (95% CI) | median exact MAE |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 12 | 101 | 0.590924 | 0.766154 | 0.910180 | 0.866958 | 0.911405 | 0.960988 | 44/101 (0.343–0.533) | 0.006067553 |
| 18 | 101 | 0.766154 | 0.910180 | 0.910180 | 0.911405 | 0.960988 | 0.960988 | 52/101 (0.419–0.610) | 0.00053106263 |
| 22 | 101 | 0.766154 | 0.910180 | 0.910180 | 0.911405 | 0.960988 | 0.960988 | 54/101 (0.438–0.629) | 0.00011989787 |
| 24 | 101 | 0.766154 | 0.766154 | 0.910180 | 0.911405 | 0.911405 | 0.960988 | 44/101 (0.343–0.533) | 6.1937514e-05 |

Tissue rankings are exploratory because selection and recovery use the same 20 labelled samples. No compatibility default or manuscript text was changed.
