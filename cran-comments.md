## Submission

This release prepares myTAI for S7 1.0.0, which the S7 maintainers plan to submit to CRAN on 2026-11-16 (https://github.com/RConsortium/S7/issues/786). It also adds a new constructor and fixes bugs.

## R CMD check results

0 errors | 1 warning | 0 notes

* checking for code/documentation mismatches ... WARNING
  Codoc mismatches for `BulkPhyloExpressionSet` and `ScPhyloExpressionSet` (`identities_label` default).

  This warning only occurs with S7 0.2.2. With S7 0.2.2, generated constructors use the parent class's property default instead of the subclass override. S7 1.0.0 fixes this, and our documentation matches the S7 1.0.0 constructors. The S7 maintainer contributed this update (https://github.com/drostlab/myTAI/pull/80). With S7 1.0.0, the check gives 0 warnings.
