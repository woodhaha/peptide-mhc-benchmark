#!/bin/bash
set -e
SRC="/mnt/d/Researching/Peptide epitope/docking"
WORK="$HOME/hdock_work"
DEST="$SRC/hdock_models"

rm -rf "$WORK"
mkdir -p "$WORK"
mkdir -p "$DEST"

cd "$WORK"
cp "$SRC/HDOCKlite-v1.1/createpl" .
cp "$SRC/receptor_1DUZ.pdb" .
for name in KRAS_G12V_neo KRAS_G12V_wt p53_R248W_neo p53_R248W_wt KRAS_G12V_enh KRAS_G12V_enh_wt; do
    cp "$SRC/${name}.out" .
done
chmod +x createpl

echo "=== Building top 10 docking models ==="
for name in KRAS_G12V_neo KRAS_G12V_wt p53_R248W_neo p53_R248W_wt KRAS_G12V_enh KRAS_G12V_enh_wt; do
    echo -n "$name: "
    ./createpl "${name}.out" "${name}_top10.pdb" -nmax 10 -complex -models 2>/dev/null
    if [ -f "${name}_top10.pdb" ]; then
        models=$(grep -c '^MODEL' "${name}_top10.pdb" 2>/dev/null || echo 0)
        size=$(wc -c < "${name}_top10.pdb")
        echo "OK (${models} models, ${size} bytes)"
        cp "${name}_top10.pdb" "$DEST/"
    else
        echo "FAILED"
    fi
done

echo ""
echo "=== Models in $DEST ==="
ls -la "$DEST/"
