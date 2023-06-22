---
name: Bug report
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behaviour, e.g.
```r
PlotSignature(PhyloExpressionSetExample, 
              measure       = "TAI", 
              permutations  = 100,
              TestStatistic = "FlatLineTest",
              ylab = "Transcriptome Age Index")
```
If the bug relates to the input data (e.g. `MyPhyloExpressionSet`), please provide some information too, i.e.
```r
utils::head(MyPhyloExpressionSet)
```
and/or
```r
utils::str(MyPhyloExpressionSet)
```

**Expected behaviour**
A clear and concise description of what you expected to happen.

**Screenshots or code**
If applicable, add screenshots or
```
code
```
to help explain your problem.

Please also feel free to [refer to the `myTAI` code in the issue](https://docs.github.com/en/issues/tracking-your-work-with-issues/creating-an-issue#creating-an-issue-from-code)

**Session info:**
- Please note session info in R
```r
utils::sessionInfo()
```
**Additional context**
Add any other context about the problem here.
