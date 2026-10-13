# Site de cours — STT-7360

Site du cours STT-7360 Méthodes biostatistiques (hiver 2027) : accueil, plan de cours, calendrier, projet de recherche et études de cas, générés avec [Quarto](https://quarto.org).

## Développement local

```bash
quarto preview          # prévisualise le site
quarto render --to html # génère le site dans _site/
```

## Source du contenu

Le contenu des pages est tiré de `plan_STT-7360_H27.tex`. Le PDF compilé (`plan_STT-7360_H27.pdf`) est versionné et offert en téléchargement sur le site. Après une modification du `.tex`, recompiler le PDF localement et reporter les changements dans les pages `.qmd` concernées.

## Publication

Chaque push sur `main` déclenche `.github/workflows/publish.yml`, qui rend le site et le publie sur GitHub Pages. Dans les réglages du dépôt, *Settings → Pages → Source* doit être réglé sur « GitHub Actions ».

## Structure

```
.
├── _quarto.yml            # configuration du site
├── index.qmd              # accueil
├── plan-de-cours.qmd      # objectifs, principes, évaluation, compétences, références
├── calendrier.qmd         # les 15 semaines
├── projet.qmd             # mini-projet de recherche : format, jalons, pistes, ADEMP
├── etudes-de-cas.qmd      # format et description des trois études de cas
├── styles.css             # thème commun aux sites de cours (copié de STT-4300)
├── date-modification-bas-page.html
├── images/simulation.svg  # illustration de l'accueil et favicon
└── plan_STT-7360_H27.{tex,pdf}
```
