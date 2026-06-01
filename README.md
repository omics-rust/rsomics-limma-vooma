# rsomics-limma-vooma

vooma-style mean-variance modelling for log-expression (microarray-like) data.

`vooma` ("mean-variance modelling at the observational level for arrays") is the
counts-free analogue of `voom`. It estimates the mean-variance trend of a
log-expression matrix and turns it into per-observation precision (inverse
variance) weights, which downstream `lmFit`/`eBayes` use to model
heteroscedasticity. Unlike `voom`, there are no counts and no library sizes —
the input is already on the log-expression scale.

```
rsomics-limma-vooma expr.tsv --design design.tsv -o weights.tsv
```

- `expr.tsv` — log-expression matrix, header = sample ids, column 1 = gene ids.
- `design.tsv` — model matrix, header = coefficient names, column 1 = sample ids
  (one row per sample, in the column order of `expr.tsv`).
- `weights.tsv` — precision-weights matrix (same shape as the expression
  matrix). `-o -` writes to stdout.
- `--trend trend.tsv` — optionally also write the fitted mean-variance line
  (average log-expression vs sqrt residual SD).

## Method

`lmFit(expr, design)` fits each gene by ordinary least squares, giving a
residual standard deviation `sigma` and a row mean `Amean`. The trend is a
LOWESS of `sqrt(sigma)` against `Amean` (span from `chooseLowessSpan` with
`small.n=50, min.span=0.3, power=1/3`). Each observation's predicted mean
`X·beta` is mapped through that trend to a predicted SD; the precision weight is
`1 / SD^4`.

## Origin

This crate is an independent Rust reimplementation of `vooma` from the **limma**
package, based on:

- The published method: Law CW (2013), *Precision weights for gene expression
  analysis*, PhD Thesis, University of Melbourne
  (<http://hdl.handle.net/11343/38150>); Law et al. (2014), *Genome Biology*
  15:R29 (the `voom` companion method), DOI 10.1186/gb-2014-15-2-r29.
- The public limma documentation for `vooma` and `chooseLowessSpan`.
- Black-box behaviour testing against the limma binary (weights matched to
  ~1e-15 relative on the test fixtures).

No source code from the GPL upstream was used as reference during
implementation. Test fixtures are independently generated.

License: MIT OR Apache-2.0.
Upstream credit: limma <https://bioconductor.org/packages/limma/> (GPL >= 2).
