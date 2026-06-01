//! Per-gene ordinary least squares (limma lmFit, method="ls").
//!
//! Householder QR of the design once, reused for every gene. Yields per-gene
//! coefficients, residual sd (sigma) and row mean (Amean).

use rsomics_common::{Result, RsomicsError};

pub struct Qr {
    n: usize,
    p: usize,
    qr: Vec<Vec<f64>>,
    rdiag: Vec<f64>,
}

impl Qr {
    fn new(x: &[Vec<f64>]) -> Result<Qr> {
        let n = x.len();
        let p = x[0].len();
        let mut qr: Vec<Vec<f64>> = x.to_vec();
        let mut rdiag = vec![0.0; p];
        for k in 0..p {
            let mut nrm = 0.0f64;
            for row in qr.iter().take(n).skip(k) {
                nrm = nrm.hypot(row[k]);
            }
            if nrm == 0.0 {
                return Err(RsomicsError::InvalidInput(
                    "design matrix is rank-deficient".into(),
                ));
            }
            if qr[k][k] < 0.0 {
                nrm = -nrm;
            }
            for row in qr.iter_mut().take(n).skip(k) {
                row[k] /= nrm;
            }
            qr[k][k] += 1.0;
            for j in (k + 1)..p {
                let mut s = 0.0;
                for row in qr.iter().take(n).skip(k) {
                    s += row[k] * row[j];
                }
                s = -s / qr[k][k];
                for row in qr.iter_mut().take(n).skip(k) {
                    let add = s * row[k];
                    row[j] += add;
                }
            }
            rdiag[k] = -nrm;
        }
        Ok(Qr { n, p, qr, rdiag })
    }

    #[allow(clippy::needless_range_loop)]
    fn qty(&self, y: &mut [f64]) {
        for k in 0..self.p {
            let mut s = 0.0;
            for i in k..self.n {
                s += self.qr[i][k] * y[i];
            }
            s = -s / self.qr[k][k];
            for i in k..self.n {
                y[i] += s * self.qr[i][k];
            }
        }
    }

    /// (beta[p], residual sum of squares).
    fn solve(&self, y: &[f64]) -> (Vec<f64>, f64) {
        let mut qty = y.to_vec();
        self.qty(&mut qty);
        let rss: f64 = qty[self.p..].iter().map(|&e| e * e).sum();
        let mut beta = vec![0.0; self.p];
        for j in (0..self.p).rev() {
            beta[j] = qty[j];
            for k in (j + 1)..self.p {
                beta[j] -= self.qr[j][k] * beta[k];
            }
            beta[j] /= self.rdiag[j];
        }
        (beta, rss)
    }
}

pub struct Fit {
    /// [gene][coef]
    pub coefficients: Vec<Vec<f64>>,
    /// residual sd per gene
    pub sigma: Vec<f64>,
    pub amean: Vec<f64>,
}

pub fn lm_fit(y: &[Vec<f64>], x: &[Vec<f64>]) -> Result<Fit> {
    let n = x.len();
    let p = x[0].len();
    if n < p {
        return Err(RsomicsError::InvalidInput(format!(
            "design has {n} samples < {p} coefficients (rank-deficient)"
        )));
    }
    if y.iter().any(|row| row.len() != n) {
        return Err(RsomicsError::InvalidInput(
            "expression samples do not match design rows".into(),
        ));
    }
    let df = (n - p) as f64;
    if df < 1.0 {
        return Err(RsomicsError::InvalidInput(
            "residual degrees of freedom < 1 (need more samples than coefficients)".into(),
        ));
    }
    let qr = Qr::new(x)?;

    let mut coefficients = Vec::with_capacity(y.len());
    let mut sigma = Vec::with_capacity(y.len());
    let mut amean = Vec::with_capacity(y.len());
    for row in y {
        let (beta, rss) = qr.solve(row);
        coefficients.push(beta);
        sigma.push((rss / df).sqrt());
        amean.push(row.iter().sum::<f64>() / n as f64);
    }

    Ok(Fit {
        coefficients,
        sigma,
        amean,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_group_means() {
        let x = vec![
            vec![1.0, 0.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![1.0, 1.0],
        ];
        let y = vec![vec![1.0, 3.0, 5.0, 7.0]];
        let f = lm_fit(&y, &x).unwrap();
        assert!((f.coefficients[0][0] - 2.0).abs() < 1e-9);
        assert!((f.coefficients[0][1] - 4.0).abs() < 1e-9);
        assert!((f.amean[0] - 4.0).abs() < 1e-9);
    }
}
