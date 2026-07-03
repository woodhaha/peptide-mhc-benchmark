#!/usr/bin/env python3
"""
Submit all 6 neoepitope-MHC docking jobs to HDOCK.
Usage: python submit_hdock.py
"""
import requests, re, time, json, sys

URL = 'http://hdock.phys.hust.edu.cn/'

JOBS = [
    # (job_name, receptor_file, ligand_file, description)
    ('KRAS_G12V_neo', 'receptor_1DUZ.pdb', 'peptide_KRAS_G12V_neo.pdb',
     'K-E63 salt bridge: K at P2 compensated by GLU63'),
    ('KRAS_G12V_wt',  'receptor_1DUZ.pdb', 'peptide_KRAS_G12V_wt.pdb',
     'WT control: K at P2, G at P9 — no compensation'),
    ('p53_R248W_neo', 'receptor_1DUZ.pdb', 'peptide_p53_R248W_neo.pdb',
     'P7 Trp burial: R->W at P7 improves MHC complementarity'),
    ('p53_R248W_wt',  'receptor_1DUZ.pdb', 'peptide_p53_R248W_wt.pdb',
     'WT control: R at P7 — charged, solvent-exposed'),
    ('KRAS_G12V_enh', 'receptor_1DUZ.pdb', 'peptide_KRAS_G12V_enh.pdb',
     'Canonical anchor control: V at P2, V at P9'),
    ('KRAS_G12V_enh_wt','receptor_1DUZ.pdb','peptide_KRAS_G12V_enh_wt.pdb',
     'WT control for canonical anchor comparison'),
]

print('=' * 60)
print('  HDOCK SUBMISSION — 6 neoepitope docking jobs')
print('=' * 60)
print(f'  Server: {URL}')
print(f'  Receptor: receptor_1DUZ.pdb (HLA-A*02:01, 1.80A)')
print()

results = {}

for jobname, receptor_f, ligand_f, desc in JOBS:
    print(f'[{len(results)+1}/{len(JOBS)}] {jobname}')
    print(f'  {desc}')

    try:
        with open(receptor_f, 'rb') as fr, open(ligand_f, 'rb') as fl:
            files = {
                'pdbfile1': (receptor_f, fr, 'application/octet-stream'),
                'pdbfile2': (ligand_f, fl, 'application/octet-stream'),
            }
            r = requests.post(URL, files=files, data={'jobname': jobname, 'hdock': '1'}, timeout=120)

        # Extract job ID from response
        job_id = None
        text = r.text
        # Try multiple patterns
        patterns = [
            r'hdock_[a-zA-Z0-9]+',
            r'[a-f0-9]{16,32}',
            r'job[_-]?id[\"\\s:=]+([a-zA-Z0-9_]+)',
            r'Your job ID[:\s]+([a-zA-Z0-9_]+)',
            r'<b>([A-Z0-9]{6,20})</b>',
            r'H[Dd]ock[_-]?([a-zA-Z0-9]+)',
        ]
        for pat in patterns:
            matches = re.findall(pat, text, re.IGNORECASE)
            if matches:
                job_id = matches[0] if isinstance(matches[0], str) else matches[0][0]
                break

        # Save response snippet for debugging
        html_preview = text[:1000].replace('\n', ' ')[:300]

        results[jobname] = {
            'status': r.status_code,
            'job_id': job_id or 'check manually',
            'desc': desc,
            'html_snippet': html_preview,
        }
        status = 'OK' if r.status_code == 200 else f'HTTP {r.status_code}'
        print(f'  -> {status} | Job ID: {results[jobname]["job_id"]}')
        if not job_id:
            print(f'  -> Response preview: {html_preview[:200]}')

    except Exception as e:
        results[jobname] = {'status': 'error', 'error': str(e)[:100], 'desc': desc}
        print(f'  -> ERROR: {e}')

    time.sleep(3)  # Polite delay

# Summary
print(f'\n{"=" * 60}')
print('  SUBMISSION SUMMARY')
print('=' * 60)

result_urls = {}
for name, r in results.items():
    marker = '✅' if r['status'] == 200 else '❌'
    jid = r.get('job_id', 'N/A')
    print(f'{marker} {name:<22s}  ID: {jid}')
    if jid and jid != 'check manually':
        result_urls[name] = f'http://hdock.phys.hust.edu.cn/results/{jid}'

# Save
with open('hdock_submission_results.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f'\nSaved to hdock_submission_results.json')
print(f'\nCheck results at: http://hdock.phys.hust.edu.cn/')
if result_urls:
    print('Result URLs (when ready):')
    for name, url in result_urls.items():
        print(f'  {name}: {url}')
