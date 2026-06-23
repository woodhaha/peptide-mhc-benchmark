#!/bin/bash
SRC="/mnt/d/Researching/Peptide epitope/docking"
WORK="$HOME/hdock_final2"
DEST="$SRC/hdock_models"
rm -rf "$WORK" && mkdir -p "$WORK" && mkdir -p "$DEST"
cd "$WORK"
cp "$SRC/HDOCKlite-v1.1/createpl" .
cp "$SRC/receptor_1DUZ.pdb" .
for name in KRAS_G12V_neo KRAS_G12V_wt p53_R248W_neo p53_R248W_wt KRAS_G12V_enh KRAS_G12V_enh_wt; do
    cp "$SRC/${name}.out" .
    cp "$SRC/peptide_${name}.pdb" . 2>/dev/null || true
done
chmod +x createpl
echo "=== Complex models (receptor + ligand) ==="
for name in KRAS_G12V_neo KRAS_G12V_wt p53_R248W_neo p53_R248W_wt KRAS_G12V_enh KRAS_G12V_enh_wt; do
    echo -n "$name: "
    ./createpl "${name}.out" "${name}_complex_top10.pdb" -nmax 10 -complex 2>/dev/null
    if [ -f "${name}_complex_top10.pdb" ]; then
        s=$(wc -c < "${name}_complex_top10.pdb")
        echo "OK (${s} bytes)"
        cp "${name}_complex_top10.pdb" "$DEST/"
    else
        echo "FAILED"
    fi
done
echo "Done"
ls -la "$DEST/"*complex*