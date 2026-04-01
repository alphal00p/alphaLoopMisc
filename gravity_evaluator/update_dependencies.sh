#!/usr/bin/env bash
#
# First update gammaloop python API
cd /Users/vjhirsch/Documents/Work/gammaloop_main
maturin build -m crates/gammaloop-api/Cargo.toml --features=ufo_support,python_api --profile=release
python3 -m pip install /Users/vjhirsch/Documents/Work/gammaloop_main/target/wheels/gammaloop-0.3.3-cp312-abi3-macosx_11_0_arm64.whl --upgrade
# Then update symbolica community dependency
cd /Users/vjhirsch/HEP_programs/TMP/sbc
maturin build --release
python3 -m pip install /Users/vjhirsch/HEP_programs/TMP/sbc/target/wheels/symbolica-1.4.0-cp37-abi3-macosx_11_0_arm64.whl --upgrade
python3 -m pip install prettytable --upgrade
python3 -m pip install progressbar2 --upgrade
