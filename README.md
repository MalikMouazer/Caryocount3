# Caryocount

Caryocount is an open-source, rule-based application for parsing cytogenetic
formulas and applying two anomaly-counting systems: ISCN 2024 and Jondreville
2020. Unlike a statistical or machine-learning model, it follows an explicit,
ordered catalogue of rules. The same input and software version therefore
produce the same result, and each reported count can be traced to the rule that
was applied.

- **Online application:** https://caryocount3.streamlit.app/
- **Source code:** https://github.com/MalikMouazer/Caryocount3
- **Test spreadsheet:** https://docs.google.com/spreadsheets/d/1_QPAzEr3PNHaNVu8Qxuv_hWCUoBQqqRcNSnX-Q0Egm8/edit?gid=98233212#gid=98233212
- **Official ISCN 2024 reference:** https://doi.org/10.1159/000538512
- **ISCN 2024 errata:** https://karger.com/cgr/article/165/2/99/920397/Erratum and https://karger.com/cgr/article/166/1/64/939149/Erratum-ISCN-2024-An-International-System-for

The example formulas in the test spreadsheet are derived from examples in the
official ISCN guide and are provided for software testing and demonstration.
Users should consult the licensed ISCN publication and its errata for the
authoritative nomenclature. Caryocount supports counting and review; it does not
replace cytogenetic interpretation or clinical validation.

## Quick start

### Use the hosted application

1. Open https://caryocount3.streamlit.app/.
2. Select English or French.
3. Enter one ISCN formula for an immediate analysis, or upload a CSV/XLS/XLSX
   file containing a `Formule` column for batch analysis.
4. Review the ISCN 2024 and Jondreville 2020 totals and use the `(?)` control to
   inspect the ordered rule path and the selected rule.
5. Download the batch results as an Excel file when required.

Simple formulas for a first check:

```text
46,XX[20]
47,XX,+8[20]
46,XY,t(9;22)(q34;q11.2)[20]
46,XX,del(5)(q13q33)[20]
```

More complex, multi-clone examples are available in the public test
spreadsheet linked above.

### Run locally

Python 3.10 or later is recommended.

```bash
git clone https://github.com/MalikMouazer/Caryocount3.git
cd Caryocount3
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
streamlit run app.py
```

The application will print a local URL, usually `http://localhost:8501`.

## How the method works

The input formula is normalized and divided into clones and individual
abnormalities. Condensed notations, multiplicities, constitutional findings,
repeated-clone notation and relationships between structural abnormalities are
then identified. Each abnormality is evaluated against a fixed priority list.
The first applicable rule determines its score; rules assigning zero prevent
events already represented elsewhere from being counted twice. The application
returns anomaly-level decisions and totals for both supported scoring systems.

The rule engine is deterministic: it contains no trained model, probabilistic
classification or generative component. Rule identifiers, priorities, default
scores, explanations and technical checks are exposed in the application. This
makes a result auditable and reproducible for a given input and code version.

## Maintaining the scoring rules

The business rules are centralized in `My_expert_karyo_functions.py`:

- `RULE_CATALOG` is the canonical source for stable rule IDs, scoring systems,
  default scores, labels and explanations.
- `RULE_PRIORITY` records the actual evaluation order for each scoring system.
- `RULE_TECHNICAL_CHECKS` documents the condition implemented for every rule.
- `rules_catalog_reference.csv` is the version-controlled snapshot displayed
  and checked by the application.
- `RuleDecision` records the rule actually applied during calculation.
- `analyser_formule(..., debug=True)` exposes rule IDs and technical
  explanations in its output.

This separation keeps the software maintainable: clinical wording can be
reviewed independently, while changes to scoring logic remain explicit,
version-controlled and reviewable in the open-source repository.

After adding, modifying or removing a rule:

1. Update its condition and score in the relevant calculation function.
2. Update `RULE_CATALOG`, `RULE_TECHNICAL_CHECKS` and `RULE_PRIORITY` together.
3. Add or update examples covering the changed behaviour.
4. Regenerate the internal catalogue:

   ```bash
   python scripts/sync_rule_catalog_reference.py
   ```

5. Run `validate_rule_catalog_integrity()` and compare the generated catalogue
   with the public catalogue before committing the change.

### Public rule catalogue

A published Google Sheet may provide reviewed French and English wording. Its
supported columns are:

- `Rule ID`
- `Libellé` and `Explication`
- `Label EN`, `Explanation EN` and `Technical criterion EN`
- `Référentiel`, `Score par défaut` and `Critère technique` for reference

To use a different published sheet in a deployment, set
`RULE_CATALOG_SHEET_URL` in `app.py` to its published URL:

```python
RULE_CATALOG_SHEET_URL = "https://docs.google.com/spreadsheets/d/.../edit?gid=0"
```

The sheet cannot silently add or remove rules, change IDs, alter scores or
replace executable criteria. Unknown IDs or structural differences are
reported, and the application falls back to the version-controlled catalogue
when the remote content cannot be safely used.

## Reproducibility and scope

For reproducible reporting, retain the exact input formula, the Caryocount Git
commit, the selected scoring system and the output file. Results depend on the
implemented scope of the parser and rule catalogue. New or unusual ISCN
constructs should be reviewed by a qualified cytogeneticist and added through a
documented rule and regression example when appropriate.
