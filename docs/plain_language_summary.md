# Plain-Language Summary

For centuries, the Inka Empire administered millions of people across the Andes without alphabetic writing, recording censuses, tribute, and accounts on knotted-string devices called *khipus*. While scholars know how decimal numbers were knotted onto cords, how whole khipus organized multi-column tables, subsidiary cords, and linked museum collections remains only partially understood.

In this open-science study of 611 khipus from the Open Khipu Repository (OKR v2.1.0) and 651 Khipu Field Guide (KFG) Excel workbooks, we separated internally **pre-specified** tests from exploratory corpus-wide discovery. First, controlling for archive size and column count, tabular khipus from the same archaeological site share column-level structure under our primary **modal color match** metric (`p = 0.0128`, though sensitive to omitting Pachacámac, `p = 0.2853`; full-distribution color JSD robustness `p = 0.0034`) and primary numeric column rank (`p = 0.0010`), while geometric outliers selectively capture known calendar and paired khipus. Second, after fixing a KFG summary bug that silently dropped 2,313 subsidiary cords, five exploratory candidate subsidiary regularities replicated in direction on 36 disjoint KFG-only khipus. Third, a descriptive corpus-wide sweep surfaced two **putative** unflagged links in Berlin collections: a three-khipu split-ledger match (`VA 42518` versus `VA 63038 b` and `VA 63038(B)`) and a **putative** reversed physical duplicate (`VA 24370(B)` and `VA 66834`). Finally, testing coastal khipus against five highland colonial texts yielded zero BH-FDR matches; under the actual BH-FDR (`m = 30`) rule at 20% noise, detection power is 70.8% (`k = 5`), 59.2% (`k = 7`), and 0% for an isolated `k = 6` hit due to a `2/720` permutation floor, pointing to coastal-versus-highland geographic mismatch.

## What This Does NOT Claim

We do **not** claim to decipher Inka language, decode non-decimal narrative khipus, prove copying direction among the Berlin khipus, or establish formal external pre-registration (protocol freeze and execution occurred on the same day).

## What Remains to Be Verified

Physical museum inspection in Berlin is required to verify provenance and duplicate status, and externally pre-registered tests against coastal colonial census records are needed to evaluate cross-modal alignment directly.
