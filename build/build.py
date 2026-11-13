"""Combine core.js and template.html into one offline HTML file.

Usage:
    python build/build.py

Writes pathway_layout_scorer.html in the project root.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), 'pathway_layout_scorer.html')


def read(name):
    with open(os.path.join(HERE, name), encoding='utf-8') as f:
        return f.read()


t = read('template.html')
if '__CORE__' not in t:
    sys.exit('__CORE__ placeholder not found in template.html')
t = t.replace('__CORE__', read('core.js'), 1)

with open(OUT, 'w', encoding='utf-8', newline='\n') as f:
    f.write(t)
print(f'{len(t)} chars written to {OUT}')
