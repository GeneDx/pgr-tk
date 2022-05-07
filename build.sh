pushd rs-wfa/
rm -rf WFA
git clone https://github.com/cschin/WFA.git --depth=1
popd

cargo build -p libwfa --release
cargo build -p pgr-db --release

#pushd pgr-py-lite/
#maturin build --release
#maturin build --release --skip-auditwheel
#popd
