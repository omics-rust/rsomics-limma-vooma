use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

fn open(path: &Path) -> Result<BufReader<File>> {
    let f = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    Ok(BufReader::new(f))
}

fn parse_f64(s: &str) -> Result<f64> {
    let t = s.trim();
    t.parse::<f64>()
        .map_err(|_| RsomicsError::InvalidInput(format!("non-numeric value '{t}'")))
}

pub struct Expr {
    pub samples: Vec<String>,
    pub genes: Vec<String>,
    /// row-major [gene][sample]
    pub y: Vec<Vec<f64>>,
}

pub fn read_expr(path: &Path) -> Result<Expr> {
    let mut lines = open(path)?.lines();
    let header = lines
        .next()
        .ok_or_else(|| RsomicsError::InvalidInput("empty expression matrix".into()))?
        .map_err(RsomicsError::Io)?;
    let samples: Vec<String> = header.split('\t').skip(1).map(str::to_string).collect();
    if samples.is_empty() {
        return Err(RsomicsError::InvalidInput(
            "expression matrix needs at least one sample column".into(),
        ));
    }
    let mut genes = Vec::new();
    let mut y = Vec::new();
    for line in lines {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty() {
            continue;
        }
        let mut f = line.split('\t');
        let gene = f
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput("missing gene id".into()))?;
        let row: Vec<f64> = f.map(parse_f64).collect::<Result<_>>()?;
        if row.len() != samples.len() {
            return Err(RsomicsError::InvalidInput(format!(
                "gene '{gene}' has {} values, header declares {} samples",
                row.len(),
                samples.len()
            )));
        }
        genes.push(gene.to_string());
        y.push(row);
    }
    if genes.len() < 2 {
        return Err(RsomicsError::InvalidInput(
            "need at least two genes to fit a mean-variance trend".into(),
        ));
    }
    Ok(Expr { samples, genes, y })
}

pub struct Design {
    pub coef_names: Vec<String>,
    /// row-major [sample][coef]
    pub x: Vec<Vec<f64>>,
}

/// Design TSV: col 1 = sample id, header = coefficient names, numeric
/// model-matrix entries (one row per sample, in sample order).
pub fn read_design(path: &Path) -> Result<Design> {
    let mut lines = open(path)?.lines();
    let header = lines
        .next()
        .ok_or_else(|| RsomicsError::InvalidInput("empty design matrix".into()))?
        .map_err(RsomicsError::Io)?;
    let coef_names: Vec<String> = header.split('\t').skip(1).map(str::to_string).collect();
    if coef_names.is_empty() {
        return Err(RsomicsError::InvalidInput(
            "design matrix needs at least one coefficient column".into(),
        ));
    }
    let mut x = Vec::new();
    for line in lines {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty() {
            continue;
        }
        let mut f = line.split('\t');
        let id = f
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput("missing design row id".into()))?;
        let row: Vec<f64> = f.map(parse_f64).collect::<Result<_>>()?;
        if row.len() != coef_names.len() {
            return Err(RsomicsError::InvalidInput(format!(
                "design row '{id}' has {} values, header declares {} coefficients",
                row.len(),
                coef_names.len()
            )));
        }
        x.push(row);
    }
    if x.is_empty() {
        return Err(RsomicsError::InvalidInput("design has no rows".into()));
    }
    Ok(Design { coef_names, x })
}
