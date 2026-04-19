# bioRxiv Paper Project

This folder contains a clean LaTeX starter for a bioRxiv preprint.

## Structure

- `main.tex` — manuscript source
- `references.bib` — bibliography database
- `figures/` — figure files
- `sections/` — optional section-level `.tex` files
- `supplement/` — optional supplementary material

## Build

Use one of the following commands from this folder:

```bash
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

or

```bash
latexmk -pdf main.tex
```

## Notes for bioRxiv submission

- Upload a single PDF generated from `main.tex`.
- Keep figure captions descriptive and include statistical details.
- Include data/code links in the "Data and Code Availability" section.
- Confirm author list, affiliations, and competing interests before submission.