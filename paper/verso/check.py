#!/usr/bin/env python3
"""Checks of the generated paper: TeX log, numbering, Verso credit, and the source's markup."""
from pathlib import Path
import json
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
JOB = 'f4-fsymbols'
OUT = ROOT / '_out' / 'pdf'
failures = []

# TeX log: no errors, overfull boxes, undefined references or citations.
log = (OUT / f'{JOB}.log').read_text(errors='replace')
for pattern in ['Overfull', 'undefined', 'Rerun to get', 'There were multiply-defined labels', '! ']:
    if pattern in log:
        failures.append(f'TeX log contains {pattern!r}')

# Every label has the number computed by the renderer.
expected = json.loads((OUT / 'labels.json').read_text())
aux = (OUT / f'{JOB}.aux').read_text()
actual = {m[1]: m[2] for m in re.finditer(r'\\newlabel\{([^}]*)\}\{\{([^}]*)\}', aux)}
for label, info in expected.items():
    if actual.get(label) != info['number']:
        failures.append(f'label {label}: renderer {info["number"]}, TeX {actual.get(label)}')

# The PDF ends with the Verso credit.
text = subprocess.check_output(['pdftotext', str(OUT / f'{JOB}.pdf'), '-'], text=True)
if not re.sub(r'\s+', ' ', text).strip().endswith('Document prepared in Verso.'):
    failures.append('the PDF does not end with the Verso credit')

# The source contains no raw TeX blocks and no layout commands outside mathematics.
source = (ROOT / 'Paper' / 'Text.lean').read_text()
prose = re.sub(r'```.*?```|%%%.*?%%%', '', source, flags=re.S)
prose = re.sub(r'\$\$?(`+).*?\1', '', prose, flags=re.S)
prose = re.sub(r'(`+).*?\1', '', prose, flags=re.S)
if m := re.search(r'\\[A-Za-z]+', prose.split('=>', 1)[1]):
    failures.append(f'TeX command outside mathematics in the source: {m[0]}')
if m := re.search(r'\\(?:vspace|hspace|vskip|hskip|newpage|clearpage|pagebreak|linebreak|smallskip|medskip|bigskip|noindent)\b', source):
    failures.append(f'layout command in the source: {m[0]}')

if failures:
    print('\n'.join('FAIL: ' + f for f in failures))
    sys.exit(1)
print(f'PASS: {len(expected)} labels numbered as computed, clean TeX log, Verso credit, no TeX outside mathematics in the source')
