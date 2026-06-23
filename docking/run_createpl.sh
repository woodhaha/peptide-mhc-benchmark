#!/bin/bash
# Run HDOCKlite createpl for all 6 docking results
# Run from WSL: bash run_createpl.sh

SRC="/mnt/d/Researching/Peptide epitope/docking"
WORK=~/hdock_work
DEST="${SRC}/hdock_models"

rm -rf "${WORK}"
mkdir -p "${WORK}"
mkdir -p "${DEST}"

# Copy binaries and data
cp "${SRC}/HDOCKlite-v1.1/hdock" "${WORK}/"
cp "${SRC}/HDOCKlite-v1.1/createpl" "${WORK}/"
cp "${SRC}/receptor_1DUZ.pdb" "${WORK}/"

for name in KRAS_G12V_neo KRAS_G12V_wt p53_R248W_neo p53_R248W_wt KRAS_G12V_enh KRAS_G12V_enh_wt; do
    cp "${SRC}/${name}.out" "${WORK}/"
done

cd "${WORK}"
chmod +x hdock createpl

echo "=== HDOCKlite createpl: building top 10 models ==="
for name in KRAS_G12V_neo KRAS_G12V_wt p53_R248W_neo p53_R248W_wt KRAS_G12V_enh KRAS_G12V_enh_wt; do
    echo -n "${name}: "
    ./createpl "${name}.out" "${name}_top10.pdb" -nmax 10 -complex -models > /dev/null 2>&1
    if [ -f "${name}_top10.pdb" ]; then
        models=$(grep -c "^MODEL" "${name}_top10.pdb" 2>/dev/null)
        size=$(wc -c < "${name}_top10.pdb")
        echo "OK (${models} models, ${size} bytes)"
        cp "${name}_top10.pdb" "${DEST}/"
    else
        echo "FAILED"
    fi
done

echo ""
echo "=== Done ==="
ls -la "${DEST}/"
