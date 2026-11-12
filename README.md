# Pathway Layout Scorer
 
A browser tool that scores the layout of metabolic pathway maps against the layout rules defined in the BioHackathon 2026 metabolic map layout project. Drop one or more map files on the page to get a 0–100 score per rule, an overall score, and a map preview that marks crossings, cofactor hubs, overlaps and the detected cycle.
 
Everything runs in the browser. Files are not uploaded anywhere.
 
## Supported formats
 
The format is detected from the file content, so the extension does not matter.
 
- **SBGN-ML** (0.2 and 0.3)
- **SBML with the Layout package** (Level 3, or Level 2 with the layout annotation)
- **SBML from CellDesigner** (lines are read as straight, so bends are not counted)
- **Escher JSON**
## Rules
 
| Rule | Type | Measured as |
|---|---|---|
| Cycles look like cycles | required | circle fit × direction consistency of the largest cycle (or one you name) |
| Cofactors drawn as clones | required | cofactors drawn once but linked to more than one reaction |
| No crossings | required | line crossings per line |
| Consistent participant placement | required | share of reactions with substrates and products on opposite sides (≥ 120°) |
| Aspect ratio | required | width ÷ height |
| Avoid long lines | recommended | longest line ÷ median line |
| Avoid bends | recommended | bend points per line |
| Density | recommended | node area ÷ bounding-box area |
| No overlapping nodes | recommended | overlapping node pairs |
 
Each rule scores 100 at its pass threshold and 0 at its fail threshold. Required rules count twice in the overall score. The thresholds, weights, cofactor list and the cycle to check can all be changed under **Scoring settings** in the page. The defaults are first proposals for the group to tune.
 
Reproducibility, fixed positions, versioning and interoperable output depend on the layout tool rather than on one map, so they are not scored here.
 
The measurements are the same as those of `score_rules.py`. On the 21 SBGN-ML maps of the Metabolism Regulation gallery, both give identical values, so their CSV outputs can be combined.
 
## Files
 
- `core.js`: parsers for the four formats, the metrics and the rule scoring
- `template.html`: the web page
- `build.py`: combines both (plus optional example maps) into one HTML file
## Build
 
```
python build.py "path/to/folder/with/example/maps"
```
 
This writes `pathway_layout_scorer.html`, a single file that works offline in any browser. The folder is optional: the script looks there for `F021-ureacycle.sbgn` and `F022-bcaa.sbgn` (or else any two `.sbgn` files) to show as examples. On Windows, keep the quotes around the path.
 
To change a rule, edit `core.js` (the cofactor list and default thresholds are near the top, the rules are in `RULES`) and run `build.py` again.
 
## Output
 
**Copy CSV** (or **Download CSV**) exports one row per map, with the same metric columns as `score_rules.py` plus `score_overall` and one `score_<rule>` column per rule.
