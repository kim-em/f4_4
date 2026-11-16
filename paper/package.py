#!/usr/bin/env python3
"""Packages the generated TeX as a self-contained upload and checks that it compiles alone."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import tempfile
import zipfile

PAPER = Path(__file__).resolve().parent
JOB = 'f4-fsymbols'
out = PAPER / 'verso' / '_out' / 'submission'
shutil.rmtree(out, ignore_errors=True)
out.mkdir(parents=True)
source = (PAPER / f'{JOB}.tex').read_text()
if '\\bibliography{' in source:
    raise SystemExit('The TeX source must include its bibliography')
archive = out / 'source.zip'
with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED) as z:
    z.writestr(f'{JOB}.tex', source)

# Compile the upload in an empty directory with pdflatex only, as the venue does.
with tempfile.TemporaryDirectory() as tmp:
    with zipfile.ZipFile(archive) as z:
        z.extractall(tmp)
    for _ in range(3):
        subprocess.run(['pdflatex', '-interaction=nonstopmode', '-halt-on-error', f'{JOB}.tex'],
                       cwd=tmp, check=True, stdout=subprocess.DEVNULL)
    log = (Path(tmp) / f'{JOB}.log').read_text(errors='replace')
    for bad in ['undefined', 'Overfull', 'Rerun to get']:
        if bad in log:
            raise SystemExit(f'The upload does not compile cleanly: {bad!r} in the log')
    shutil.copy(Path(tmp) / f'{JOB}.pdf', out / f'{JOB}.pdf')

revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=PAPER, text=True).strip()
(out / 'package.json').write_text(json.dumps({
    'source_revision': revision,
    'entry_point': f'{JOB}.tex',
    'sha256': hashlib.sha256(archive.read_bytes()).hexdigest(),
}, indent=2) + '\n')
print(archive)
