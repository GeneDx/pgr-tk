#pushd WFA2-lib
#make all
#popd

rustup default stable
cargo build -p pgr-db --release
cargo build -p pgr-bin --release
cargo install --path pgr-bin

pushd pgr-tk/
maturin build --release
maturin build --release --skip-auditwheel
popd
