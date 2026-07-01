## Dépendances

Ce projet utilise plusieurs bibliothèques open source. Les bibliothèques suivantes sont incluses dans le code source de ce projet, mais n'ont pas été modifiées.

### Bibliothèques utilisées :

- **openpyxl**
  - Lien vers le projet : https://github.com/soxhub/openpyxl 
- **et_xmlfile**
  - Lien vers le projet : https://github.com/biydnd/et_xmlfile 

## Maintenir les règles de scoring

Les règles métier sont centralisées dans `My_expert_karyo_functions.py`.

- `RULE_CATALOG` décrit les règles visibles dans l'application : identifiant, référentiel, score par défaut, libellé et explication.
- `RuleDecision` est la décision effectivement appliquée pendant le calcul.
- `analyser_formule(..., debug=True)` ajoute les colonnes `RuleID_ISCN`, `RuleID_Jon` et les explications techniques associées.
- L'application affiche un bouton `(?)` devant chaque comptage ISCN/Jon pour montrer le parcours des règles jusqu'à la règle retenue.

Pour modifier une règle, garder le même réflexe :

1. Modifier la condition ou le score dans la fonction de calcul concernée.
2. Mettre à jour l'entrée correspondante dans `RULE_CATALOG`.
3. Lancer une comparaison avec `compare_sheet.py` sur le fichier de référence.
