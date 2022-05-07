
git submodule update --init
pushd agc
make -f makefile.release
popd

pushd rs-wfa/
rm -rf WFA
git clone https://github.com/cschin/WFA.git --depth=1
popd

rustup default stable
cargo build -p libwfa --release
cargo build -p pgr-db --release

pushd pgr-py-lite/
maturin build --release
maturin build --release --skip-auditwheel
popd
