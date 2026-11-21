# Targeted bootstrap endpoint stability check

This supplement reuses the protocol lognormal generator and deterministic sample seeds. It reruns only the two n=30 conditions with the largest observed ROS undercoverage. It does not recompute the main simulation metrics or modify their results.

- Samples per condition: 30; conditions: n=30, F(LOQ)=0.10 and 0.40.
- Paired bootstrap sizes: B=99 and B=499; common random-number seeds make the smaller run a prefix of the larger run.
- Additional bootstrap work: 35880 refits per family (MLE or ROS); 71760 total refits across MLE and ROS.
- The deterministic comparison was executed twice: the first pass aggregated endpoint changes, and this second pass reports lower and upper endpoints separately. Total compute spent across both passes was 143520 refits; only the endpoint-separated second-pass result is reported below.
- Endpoint differences are reported separately for lower and upper endpoints; coverage counts are descriptive for these 30 paired samples and are not fresh validation estimates.

| Condition | Procedure / endpoint | Median absolute change | 90th percentile | Maximum |
|---|---|---:|---:|---:|
| n30_cp0.1 | MLE parametric percentile lower endpoint | 0.0202 | 0.0599 | 0.1789 |
| n30_cp0.1 | MLE parametric percentile upper endpoint | 0.0768 | 0.1950 | 0.5175 |
| n30_cp0.1 | ROS percentile lower endpoint | 0.0213 | 0.0962 | 0.1024 |
| n30_cp0.1 | ROS percentile upper endpoint | 0.0600 | 0.2114 | 0.3114 |
| n30_cp0.1 | ROS basic lower endpoint | 0.0600 | 0.2114 | 0.3114 |
| n30_cp0.1 | ROS basic upper endpoint | 0.0213 | 0.0962 | 0.1024 |
| n30_cp0.4 | MLE parametric percentile lower endpoint | 0.0366 | 0.0679 | 0.1368 |
| n30_cp0.4 | MLE parametric percentile upper endpoint | 0.0870 | 0.2515 | 0.2800 |
| n30_cp0.4 | ROS percentile lower endpoint | 0.0204 | 0.0746 | 0.1416 |
| n30_cp0.4 | ROS percentile upper endpoint | 0.0368 | 0.2388 | 0.3485 |
| n30_cp0.4 | ROS basic lower endpoint | 0.0368 | 0.2388 | 0.3485 |
| n30_cp0.4 | ROS basic upper endpoint | 0.0204 | 0.0746 | 0.1416 |

## Empirical inclusion counts in the paired subset

| Condition | Procedure | B=99 | B=499 |
|---|---|---:|---:|
| n30_cp0.1 | MLE parametric percentile | 29/30 | 28/30 |
| n30_cp0.1 | ROS percentile | 26/30 | 27/30 |
| n30_cp0.1 | ROS basic | 26/30 | 25/30 |
| n30_cp0.4 | MLE parametric percentile | 27/30 | 28/30 |
| n30_cp0.4 | ROS percentile | 27/30 | 27/30 |
| n30_cp0.4 | ROS basic | 26/30 | 25/30 |
