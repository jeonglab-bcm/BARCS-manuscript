# Generated figures

Everything in this directory is output. The scripts under `examples/` write
here, and none of them creates the directory first, so this file exists to keep
it in the repository.

The figure PDFs themselves are **not** versioned here. Three of them are cited
by the manuscript and are versioned in the `overleaf/` submodule, which is the
Overleaf project:

| Figure | Written by |
|---|---|
| `crispulator_facs_moi_10k_f1_by_fdr.pdf` | `examples/manuscript_crispulator_figure.R` |
| `liang_longitudinal_volcano_trajectories.pdf` | `examples/manuscript_liang_figure.R` |
| `simcrispr_interaction.pdf` | `examples/simcrispr_interaction_figure.R` |

After regenerating any of those, publish it to Overleaf:

```sh
scripts/publish_figures.sh
```

That copies the three cited figures into `overleaf/figures/` and tells you what
changed. Committing and pushing inside `overleaf/` is what makes the new
version appear in Overleaf.

Every other PDF written here is a diagnostic or an exploratory benchmark plot
that the manuscript does not include. Those stay local.
