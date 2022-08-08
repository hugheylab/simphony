## R CMD check results

### Local check
`devtools::check()` result:

  0 errors ✓ | 0 warnings ✓ | 0 notes ✓

Notes:
  - We have been unable to reproduce the issue with a test failing on M1mac. Nonetheless, we have revised the test in an attempt to resolve the issue, including by setting the seed on the RNG (as requested by Uwe Ligges).

You can also see the results of R CMD check on Windows, Linux, and MacOS by going to the GitHub Actions run labeled `check-deploy` [here](https://github.com/hugheylab/simphony/actions).

## Downstream dependencies
There are no downstream dependencies for simphony.

## Tests
When checking locally, tests take about 2.7s total. If this is too long for tests let us know and we can skip some tests on CRAN.
