# Colonial Khipu Transcription Schemas (`data/colonial_schemas/`)

Each `.json` file in this directory encodes **one primary-source verified colonial Andean khipu transcription** for internally pre-specified Step 3 (`src/step3_colonial_alignment.py`). See [`DATA_LICENSE.md`](../../DATA_LICENSE.md) for the explicit licensing exclusion notice pending written permission from the *Textos Andinos* copyright holders.

## Required JSON Contract

```json
{
  "schema_id": "santa_valley_1670",
  "title": "Santa Valley Tributary Census (1670)",
  "source_citation": "Medrano & Urton (2018), Ethnohistory 65(1): 1-23",
  "region": "Santa Valley, Ancash, Peru",
  "year": 1670,
  "domain": "census",
  "role": "positive_control",
  "editorial_notes": "Document any unresolved edition discrepancies here.",
  "column_names": ["col_1", "col_2", "..."],
  "row_names": ["row_1", "row_2", "..."],
  "value_matrix": [
    [0, 0, 0]
  ],
  "structural_flags": {
    "has_moiety_split": true,
    "has_summary_row": false,
    "has_round_100_totals": false
  }
}
```

## Verification Rule
No schema enters this directory until every cell of `value_matrix` has been verified against a primary scholarly edition and any editorial discrepancies (e.g., Murra 1982 vs. Julien 1988 for Huánuco 1549) are explicitly documented in `editorial_notes`.
