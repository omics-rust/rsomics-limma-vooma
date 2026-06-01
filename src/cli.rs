use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_limma_vooma::{read_design, read_expr, vooma, write_trend, write_weights};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-limma-vooma", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// Log-expression matrix TSV: header = sample ids, col 1 = gene ids.
    pub expr: PathBuf,
    /// Design matrix TSV: header = coefficient names, col 1 = sample ids.
    #[arg(long)]
    design: PathBuf,
    /// Precision-weights matrix destination; "-" is stdout.
    #[arg(short = 'o', long, default_value = "-")]
    output: String,
    /// Also write the fitted mean-variance trend (AveLogExpr, sqrtSD) here.
    #[arg(long)]
    trend: Option<PathBuf>,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let expr = read_expr(&self.expr)?;
        let design = read_design(&self.design)?;
        let v = vooma(&expr, &design)?;

        let mut out: Box<dyn std::io::Write> = if self.output == "-" {
            Box::new(std::io::stdout().lock())
        } else {
            Box::new(std::fs::File::create(&self.output).map_err(RsomicsError::Io)?)
        };
        write_weights(&v, &mut out)?;
        drop(out);

        if let Some(tpath) = &self.trend {
            write_trend(&v, tpath)?;
        }

        if !self.common.quiet {
            eprintln!(
                "{} genes x {} samples vooma-weighted",
                v.genes.len(),
                v.samples.len()
            );
        }
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "vooma mean-variance precision weights for a log-expression matrix.",
    origin: Some(Origin {
        upstream: "limma vooma",
        upstream_license: "GPL (>=2)",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("11343/38150"),
    }),
    usage_lines: &["<expr.tsv> --design <design.tsv> [-o weights.tsv] [--trend trend.tsv]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: None,
                long: "design",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: true,
                default: None,
                description: "Design matrix TSV (header = coefficient names, col 1 = sample ids).",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("String"),
                required: false,
                default: Some("-"),
                description: "Precision-weights matrix destination; \"-\" is stdout.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "trend",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: false,
                default: None,
                description: "Also write the fitted mean-variance trend points.",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Two-group design, weights to a file",
            command: "rsomics-limma-vooma E.tsv --design design.tsv -o weights.tsv",
        },
        Example {
            description: "Weights to stdout plus the trend curve",
            command: "rsomics-limma-vooma E.tsv --design design.tsv --trend trend.tsv > weights.tsv",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
