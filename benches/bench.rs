use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_vooma(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-limma-vooma");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let expr = manifest.join("tests/golden/expr.tsv");
    let design = manifest.join("tests/golden/design.tsv");
    c.bench_function("rsomics-limma-vooma golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args([
                    expr.to_str().unwrap(),
                    "--design",
                    design.to_str().unwrap(),
                    "-o",
                    "/dev/null",
                ])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_vooma);
criterion_main!(benches);
