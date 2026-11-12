import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
t = open(os.path.join(HERE, 'template.html'), encoding='utf-8').read()
t = t.replace('__CORE__', open(os.path.join(HERE, 'core.js'), encoding='utf-8').read()) \
     .replace('__EXAMPLES__', json.dumps(ex).replace('</', '<\\/'))
open(os.path.join(HERE, 'pathway_layout_scorer.html'), 'w', encoding='utf-8').write(t)
print(len(t), "chars written,", len(ex), "examples embedded")
