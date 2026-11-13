"""Combine core.js, layout_rules.json and template.html into one offline HTML file.

Usage:
    python build/build.py

Writes pathway_layout_scorer.html in the project root.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), 'pathway_layout_scorer.html')


def read(name):
    with open(os.path.join(HERE, name), encoding='utf-8') as f:
        return f.read()


t = read('template.html')
for placeholder in ('__CORE__', '__RULES__'):
    if placeholder not in t:
        sys.exit(f'{placeholder} placeholder not found in template.html')

rules = json.loads(read('layout_rules.json'))  # fails loudly on invalid JSON
rules_js = json.dumps(rules, ensure_ascii=False).replace('</', '<\\/')

t = t.replace('__RULES__', rules_js, 1)
t = t.replace('__CORE__', read('core.js'), 1)

with open(OUT, 'w', encoding='utf-8', newline='\n') as f:
    f.write(t)
print(f'{len(t)} chars written to {OUT} ({len(rules["rules"])} rules)')
