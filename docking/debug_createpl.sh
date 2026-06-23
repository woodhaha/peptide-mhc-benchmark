#!/bin/bash
SRC="/mnt/d/Researching/Peptide epitope/docking"
WORK=~/hdock_work
DEST="${SRC}/hdock_models"

rm -rf "${WORK}"
mkdir -p "${WORK}"
mkdir -p "${DEST}"

cp "${SRC}/HDOCKlite-v1.1/hdock" "${WORK}/"
cp "${SRC}/HDOCKlite-v1.1/createpl" "${WORK}/"

cd "${WORK}"
chmod +x hdock createpl
pwd
ls -la

echo "---"
echo "createpl exists: $(test -f ./createpl && echo YES || echo NO)"
echo "createpl executable: $(test -x ./createpl && echo YES || echo NO)"
echo "ldd: $(ldd ./createpl 2>&1 | head -3)"

echo "---"
cp "${SRC}/KRAS_G12V_neo.out" "${WORK}/"
echo "out file: $(ls -la KRAS_G12V_neo.out)"
echo "Running: ./createpl KRAS_G12V_neo.out test.pdb -nmax 5 -complex -models"
./createpl KRAS_G12V_neo.out test.pdb -nmax 5 -complex -models
echo "exit: $?"
ls -la test*
