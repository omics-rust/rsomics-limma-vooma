//! vooma: mean-variance modelling at the observational level for arrays.
//!
//! Clean-room reimplementation of limma's `vooma()` for a log-expression
//! matrix and design. Reference: Law CW (2013), "Precision weights for gene
//! expression analysis", PhD Thesis, University of Melbourne
//! (hdl.handle.net/11343/38150). No limma (GPL) source was consulted; the
//! method follows the published documentation and is validated black-box
//! against the binary.
//!
//! lmFit(y, design) gives per-gene residual sd (sigma) and row mean (Amean);
//! the trend is a lowess of sqrt(sigma) on Amean; each observation's precision
//! weight is `1 / trend(fitted)^4` where fitted is its model-predicted mean.

mod fit;
mod lowess;
mod matrix;

use std::io::{BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

pub use matrix::{Design, Expr, read_design, read_expr};

/// limma `chooseLowessSpan(n, small.n=50, min.span=0.3, power=1/3)`: the span
/// vooma uses by default. (The man page's "power=1.3" is a documentation slip;
/// the binary uses the chooseLowessSpan default 1/3 — verified black-box.)
fn lowess_span(n: usize) -> f64 {
    if n <= 50 {
        return 1.0;
    }
    (0.3 + 0.7 * (50.0 / n as f64).powf(1.0 / 3.0)).min(1.0)
}

pub struct Vooma {
    pub samples: Vec<String>,
    pub genes: Vec<String>,
    /// precision weights, row-major [gene][sample]
    pub weights: Vec<Vec<f64>>,
    /// the fitted mean-variance trend, as (Amean, sqrt-sd) line points
    pub trend_x: Vec<f64>,
    pub trend_y: Vec<f64>,
}

pub fn vooma(expr: &Expr, design: &Design) -> Result<Vooma> {
    let n = expr.samples.len();
    let ng = expr.genes.len();
    if design.x.len() != n {
        return Err(RsomicsError::InvalidInput(format!(
            "design has {} rows, expression has {n} samples",
            design.x.len()
        )));
    }

    let fit = fit::lm_fit(&expr.y, &design.x)?;

    let sx = &fit.amean;
    let sy: Vec<f64> = fit.sigma.iter().map(|&s| s.sqrt()).collect();

    let mut order: Vec<usize> = (0..ng).collect();
    order.sort_by(|&a, &b| sx[a].partial_cmp(&sx[b]).unwrap());
    let lx: Vec<f64> = order.iter().map(|&i| sx[i]).collect();
    let ly: Vec<f64> = order.iter().map(|&i| sy[i]).collect();
    let delta = 0.01 * (lx[ng - 1] - lx[0]);
    let fitted_line = lowess::lowess(&lx, &ly, lowess_span(ng), 3, delta);
    let trend = lowess::Trend::new(&lx, &fitted_line);

    // weight = 1 / trend(fitted)^4, fitted = X beta per observation
    let mut weights = vec![vec![0.0; n]; ng];
    for (gi, beta) in fit.coefficients.iter().enumerate() {
        for (j, w) in weights[gi].iter_mut().enumerate() {
            let xrow = &design.x[j];
            let predicted_mean: f64 = beta.iter().zip(xrow).map(|(&b, &xij)| b * xij).sum();
            let sd = trend.eval(predicted_mean);
            *w = 1.0 / sd.powi(4);
        }
    }

    Ok(Vooma {
        samples: expr.samples.clone(),
        genes: expr.genes.clone(),
        weights,
        trend_x: lx,
        trend_y: fitted_line,
    })
}

pub fn write_weights(v: &Vooma, out: &mut dyn Write) -> Result<()> {
    let mut w = BufWriter::with_capacity(1 << 20, out);
    write!(w, "gene").map_err(RsomicsError::Io)?;
    for s in &v.samples {
        write!(w, "\t{s}").map_err(RsomicsError::Io)?;
    }
    writeln!(w).map_err(RsomicsError::Io)?;

    let mut fmt = ryu::Buffer::new();
    let mut line = String::with_capacity(v.samples.len() * 16);
    for (gene, row) in v.genes.iter().zip(&v.weights) {
        line.clear();
        line.push_str(gene);
        for &val in row {
            line.push('\t');
            line.push_str(fmt.format(val));
        }
        line.push('\n');
        w.write_all(line.as_bytes()).map_err(RsomicsError::Io)?;
    }
    w.flush().map_err(RsomicsError::Io)?;
    Ok(())
}

pub fn write_trend(v: &Vooma, path: &Path) -> Result<()> {
    let f = std::fs::File::create(path).map_err(RsomicsError::Io)?;
    let mut w = BufWriter::new(f);
    writeln!(w, "AveLogExpr\tsqrtSD").map_err(RsomicsError::Io)?;
    let mut xb = ryu::Buffer::new();
    let mut yb = ryu::Buffer::new();
    for (&x, &y) in v.trend_x.iter().zip(&v.trend_y) {
        writeln!(w, "{}\t{}", xb.format(x), yb.format(y)).map_err(RsomicsError::Io)?;
    }
    w.flush().map_err(RsomicsError::Io)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_span_caps_at_one() {
        assert_eq!(lowess_span(40), 1.0);
        assert_eq!(lowess_span(50), 1.0);
        assert!(lowess_span(1000) < 1.0 && lowess_span(1000) > 0.3);
    }

    #[test]
    fn intercept_only_weights_positive() {
        let expr = Expr {
            samples: vec!["s1".into(), "s2".into(), "s3".into(), "s4".into()],
            genes: (0..60).map(|i| format!("g{i}")).collect(),
            y: (0..60)
                .map(|i| {
                    let m = (i % 10) as f64;
                    vec![m + 0.1, m - 0.1, m + 0.2, m - 0.2]
                })
                .collect(),
        };
        let design = Design {
            coef_names: vec!["Intercept".into()],
            x: vec![vec![1.0]; 4],
        };
        let v = vooma(&expr, &design).unwrap();
        assert_eq!(v.weights.len(), 60);
        assert!(v.weights.iter().all(|r| r.iter().all(|&w| w > 0.0)));
    }
}
