## R CMD check results

### Local check
`devtools::check()` result:

  0 errors ✓ | 0 warnings ✓ | 0 notes ✓

### Online check
`devtools::check_rhub()` result:

0 errors ✓ | 0 warnings ✓ | 0 notes x

Notes:
  - We are submitting 1.0.1 with only minor updates to trigger a re-run on M1 mac check. We received notice that it was failing check on M1 but we have been unable to recreate.
  - We believe the issue with M1 was a "one-off" issue since it appeared to be a floating point precision error that isn't able to be recreated, so we hope this run will re-trigger the flow and confirm that.

You can also see the results of R CMD check on Windows, Linux, and MacOS by going to the GitHub Actions run labeled `check-deploy` [here](https://github.com/hugheylab/simphony/actions).

## Downstream dependencies
There are no downstream dependencies for simphony.

## Tests
When checking locally, tests take about 2.7s total. If this is too long for tests let us know and we can skip some tests on CRAN.
