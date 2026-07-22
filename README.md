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
- `RULE_CATALOG` reste la source canonique des IDs, du nombre de règles, des référentiels et des scores.
- `RULE_TECHNICAL_CHECKS` décrit ce que le code vérifie concrètement pour chaque règle : expression régulière, préfixe, fonction métier ou relation entre clones.
- `rules_catalog_reference.csv` est la photographie interne suivie par Git. Sa colonne `N°` matérialise l'ordre de priorité propre à chaque référentiel. Le fichier doit être régénéré dans le même commit que tout ajout, modification ou suppression de règle.
- Les lignes du CSV et leurs numéros dans l'interface suivent `RULE_PRIORITY`, donc le numéro affiché correspond exactement à l'ordre réel d'évaluation montré dans le parcours des règles.
- Une feuille Google Sheets publique peut surcharger uniquement les colonnes `Libellé` et `Explication`.
- `RuleDecision` est la décision effectivement appliquée pendant le calcul.
- `analyser_formule(..., debug=True)` ajoute les colonnes `RuleID_ISCN`, `RuleID_Jon` et les explications techniques associées.
- L'application affiche un bouton `(?)` devant chaque comptage ISCN/Jon pour montrer le parcours des règles jusqu'à la règle retenue.

### Catalogue public Google Sheets

Créer une feuille Google Sheets publiée avec ces colonnes :

- `Rule ID` : obligatoire, doit correspondre exactement à un ID existant dans `RULE_CATALOG`.
- `Critère technique` : aide à la rédaction, générée par le code, à ne pas modifier dans la feuille.
- `Libellé` : modifiable par les utilisateurs métier.
- `Explication` : modifiable par les utilisateurs métier.
- `Référentiel` et `Score par défaut` peuvent être gardés dans la feuille pour lecture, mais ils sont ignorés par l'application.

Configurer ensuite l'application avec la variable d'environnement :

```bash
RULE_CATALOG_SHEET_URL="https://docs.google.com/spreadsheets/d/.../edit?gid=0"
```

La feuille ne peut pas ajouter de règle, supprimer de règle, changer un score, changer un ID ou changer le critère technique. Si elle contient un `Rule ID` inconnu, le chargement est refusé et l'application revient au catalogue interne.

### Bonne conduite quand une règle change

Si la modification est uniquement clinique ou rédactionnelle :

1. Modifier `Libellé` et/ou `Explication` dans le Google Sheet.
2. Ne pas modifier `Rule ID`, `Référentiel`, `Score par défaut` ni `Critère technique`.

Si la logique de calcul change :

1. Modifier la condition ou le score dans la fonction de calcul concernée.
2. Mettre à jour l'entrée correspondante dans `RULE_CATALOG` si le libellé interne, le score par défaut ou le référentiel changent.
3. Mettre à jour `RULE_TECHNICAL_CHECKS` avec le test concret réellement codé : regex utilisée, préfixe testé, normalisation appliquée, compteurs, relation entre clones ou règle d'exclusion.
4. Vérifier que l'ID existe aussi dans `RULE_PRIORITY`, dans l'ordre réel d'évaluation.
5. Lancer `validate_rule_catalog_integrity()` pour contrôler la cohérence catalogue/priorités.
6. Vérifier qu'il existe un critère technique pour chaque `Rule ID`.
7. Lancer une comparaison avec `compare_sheet.py` sur le fichier de référence.

Si une nouvelle règle est ajoutée :

1. Créer un `Rule ID` stable et explicite.
2. Ajouter le `RuleSpec` dans `RULE_CATALOG`.
3. Ajouter le critère dans `RULE_TECHNICAL_CHECKS`.
4. Ajouter l'ID dans `RULE_PRIORITY` à l'endroit où la règle est réellement testée.
5. Ajouter ou mettre à jour des exemples de référence couvrant cette règle.

Après toute création, modification ou suppression de règle, régénérer le CSV interne :

```bash
python scripts/sync_rule_catalog_reference.py
```

L'application compare automatiquement ce CSV avec le Google Sheet dans le tableau
principal, numéroté séparément (`1_ISCN`, `1_JON`, etc.) pour chaque référentiel.
L'en-tête et la colonne `Rule ID` restent figés pendant le défilement. Une règle
absente du fichier distant apparaît en bleu à sa
position canonique; une règle présente uniquement dans le distant apparaît en
rouge. Les cellules distantes `Libellé` ou `Explication` différentes apparaissent
en jaune et affichent directement les versions distante et interne à harmoniser.

Si deux situations ont le même score mais nécessitent des explications cliniques
différentes, elles doivent avoir deux `Rule ID` différents. Exemple :
`CONSTITUTIONAL_GAIN`, `CONSTITUTIONAL_SUSPECT` et
`CONSTITUTIONAL_CERTAIN` valent tous 0, mais ils restent séparés pour que les
cliniciens puissent rédiger un libellé et une explication adaptés dans le
Google Sheet.
