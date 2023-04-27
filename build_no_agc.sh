#pushd WFA2-lib
#make all
#popd

rustup default stable
cargo build -p pgr-db --release --no-default-features
cargo build -p pgr-bin --release --no-default-features
cargo install --path pgr-bin

pushd pgr-tk/
maturin build --release --no-default-features
maturin build --release --skip-auditwheel --no-default-features
popd
