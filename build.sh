pushd WFA2-lib
make all
popd

rustup default stable
cargo build -p pgr-db --release
cargo build -p pgr-bin --release

pushd pgr-tk/
maturin build --release
maturin build --release --skip-auditwheel
popd
