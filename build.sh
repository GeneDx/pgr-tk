pushd rs-wfa/
rm -rf WFA
git clone https://github.com/cschin/WFA.git --depth=1
popd

rustup default stable
cargo build -p libwfa --release
cargo build -p pgr-db --release
cargo build -p pgr-bin --release

pushd pgr-tk/
maturin build --release
maturin build --release --skip-auditwheel
popd
