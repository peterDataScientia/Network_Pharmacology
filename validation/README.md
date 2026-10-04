# EGCG–RISI validation fixture

This folder contains the 32 shared EGCG–RISI targets used by the uploaded reference workflow.

## Validation goals

Run the publication enrichment engine twice:

1. `background_mode=default` to reproduce the older systematic notebook design.
2. `background_mode=annotated` to test the newer explicitly defined all-human OrgDb universe.

The two result sets must be treated as different statistical analyses and must never be mixed.

## Older executed notebook reference

The older systematic notebook reported, after BH-adjusted P < 0.05:

- GO-BP significant: 1432; reduced: 49
- GO-CC significant: 23; reduced: 17
- GO-MF significant: 54; reduced: 26
- Reactome significant: 276; reduced: 91

These values are a historical validation reference for the notebook environment, not universal constants. Annotation and package-version changes can alter counts. Validation should therefore inspect:
- identifier mapping,
- background definition,
- raw term identities and statistics,
- reduced representative terms,
- package/database versions,
not just total counts.

## Run locally

```bash
Rscript scripts/install_publication_enrichment.R

python validation/run_egcg_validation.py
```

The validator writes separate result directories for default and annotated backgrounds.
