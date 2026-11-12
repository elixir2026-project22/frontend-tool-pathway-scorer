/* Layout scoring core: a port of score_rules.py (Biohackathon 2026, metabolic map layout rules).
   Works in the browser (DOMParser) and in node (pass a parsed XML Document).

   Rules (kept in step with sbml_pathway_layout.py / svg_layout_metrics.py):
     hard  Reproducible rendering            (from the Python tools; needs opts.reference = a second rendering)
           Labels at fixed positions per type (from the Python tools; needs label geometry in the file)
           Cycles look like cycles            (Python definition: drawn as a closed loop in cycle order)
           Cofactors drawn as clones
           No crossings
           No line through a node             (from the Python tools)
           Consistent participant placement
           Aspect ratio
     soft  Avoid long lines, Avoid bends, Density, No overlapping nodes */
(function (root) {
  "use strict";

  const CURRENCY_DEFAULT = [
    "atp", "adp", "amp", "gtp", "gdp", "gmp", "utp", "udp", "ump", "ctp", "cdp",
    "nad", "nad+", "nadh", "nadp", "nadp+", "nadph", "fad", "fadh2", "fmn",
    "h2o", "water", "co2", "o2", "h+", "h", "pi", "ppi", "phosphate",
    "pyrophosphate", "diphosphate", "coa", "coa-sh", "coash", "nh3", "nh4+",
    "ammonia", "hco3-", "bicarbonate", "h2o2", "e-", "camp", "sam", "sah",
    "acceptor", "reduced acceptor", "glutathione", "nh4",
  ];
  const SUB = { "₀": "0", "₁": "1", "₂": "2", "₃": "3", "₄": "4", "₅": "5", "₆": "6", "₇": "7", "₈": "8", "₉": "9", "⁺": "+", "⁻": "-" };
  const PROCESS = new Set(["process", "omitted process", "uncertain process", "association", "dissociation"]);
  const METABOLITE = new Set(["simple chemical", "simple chemical multimer", "nucleic acid feature", "unspecified entity"]);
  const INPUT = new Set(["consumption"]);
  const OUTPUT = new Set(["production"]);
  const MODIFIER = new Set(["catalysis", "modulation", "stimulation", "inhibition", "necessary stimulation", "absolute inhibition"]);
  const SKIP_AREA = new Set(["compartment", "submap"]);

  function normName(s) {
    s = (s || "").replace(/[₀-₉⁺⁻]/g, (c) => SUB[c]).trim().toLowerCase();
    return s.replace(/\s+/g, " ");
  }

  // ------------------------------------------------------------- DOM helpers
  const ln = (el) => el.localName || String(el.nodeName).split(":").pop();
  function kids(el) {
    const out = [];
    for (let c = el.firstChild; c; c = c.nextSibling) if (c.nodeType === 1) out.push(c);
    return out;
  }
  function attr(el, name) {
    if (!el || !el.attributes) return null;
    for (let i = 0; i < el.attributes.length; i++) {
      const a = el.attributes[i];
      const n = a.localName || String(a.name).split(":").pop();
      if (n === name) return a.value;
    }
    return null;
  }
  function allEls(doc) { return Array.from(doc.getElementsByTagName("*")); }
  function descendants(el, out = []) {
    out.push(el);
    for (const k of kids(el)) descendants(k, out);
    return out;
  }
  const kid = (el, name) => kids(el).find((c) => ln(c) === name) || null;

  // ------------------------------------------------------------- model
  class Model {
    constructor(name, format) {
      this.name = name; this.format = format;
      this.glyphs = new Map(); this.arcs = [];
    }
    addGlyph(g) {
      g.cx = g.x + g.w / 2; g.cy = g.y + g.h / 2;
      g.keys = g.keys || [];
      this.glyphs.set(g.id, g);
    }
    nodes(classes, exclude) {
      const out = [];
      for (const g of this.glyphs.values())
        if ((!classes || classes.has(g.cls)) && !(exclude && exclude.has(g.cls))) out.push(g);
      return out;
    }
  }

  // ------------------------------------------------------------- SBGN-ML
  function parseSBGN(doc, name) {
    const mp = allEls(doc).find((e) => ln(e) === "map");
    if (!mp) throw new Error("not an SBGN-ML file (no <map> element)");
    const m = new Model(name, "SBGN-ML");
    const owner = new Map();
    for (const g of kids(mp)) {
      if (ln(g) !== "glyph") continue;
      const gid = attr(g, "id");
      const bb = kid(g, "bbox");
      if (!bb) continue;
      const f = (k) => parseFloat(attr(bb, k) || 0);
      const lab = kid(g, "label");
      const lbb = lab ? kid(lab, "bbox") : null;
      m.addGlyph({
        id: gid, cls: attr(g, "class"), x: f("x"), y: f("y"), w: f("w"), h: f("h"),
        label: lab ? attr(lab, "text") || "" : "",
        clone: kids(g).some((c) => ln(c) === "clone"),
        labelBox: lbb ? { x: parseFloat(attr(lbb, "x") || 0), y: parseFloat(attr(lbb, "y") || 0), w: parseFloat(attr(lbb, "w") || 0), h: parseFloat(attr(lbb, "h") || 0) } : null,
      });
      for (const d of descendants(g)) {
        const t = ln(d);
        if ((t === "glyph" || t === "port") && attr(d, "id")) owner.set(attr(d, "id"), gid);
      }
    }
    for (const a of kids(mp)) {
      if (ln(a) !== "arc") continue;
      const pts = [];
      for (const c of kids(a)) {
        const t = ln(c);
        if (t === "start" || t === "next" || t === "end") pts.push([parseFloat(attr(c, "x")), parseFloat(attr(c, "y"))]);
      }
      const src = owner.get(attr(a, "source")), tgt = owner.get(attr(a, "target"));
      if (src === undefined || tgt === undefined || pts.length < 2) continue;
      m.arcs.push({ id: attr(a, "id"), cls: attr(a, "class"), src, tgt, pts, bends: pts.length - 2 });
    }
    return m;
  }

  // ------------------------------------------------------------- SBML Layout
  const MACRO_SBO = new Set(["SBO:0000252", "SBO:0000245", "SBO:0000246", "SBO:0000253", "SBO:0000297", "SBO:0000243", "SBO:0000354", "SBO:0000278", "SBO:0000250", "SBO:0000251"]);
  function bezier(p0, c1, c2, p1, n) {
    const out = [];
    for (let i = 1; i < n; i++) {
      const t = i / n, u = 1 - t;
      out.push([
        u * u * u * p0[0] + 3 * u * u * t * c1[0] + 3 * u * t * t * c2[0] + t * t * t * p1[0],
        u * u * u * p0[1] + 3 * u * u * t * c1[1] + 3 * u * t * t * c2[1] + t * t * t * p1[1],
      ]);
    }
    return out;
  }
  const pt = (el) => [parseFloat(attr(el, "x") || 0), parseFloat(attr(el, "y") || 0)];
  function sbmlCurve(curve) {
    const lst = curve && kid(curve, "listOfCurveSegments");
    const segs = lst ? kids(lst).filter((s) => ln(s) === "curveSegment") : [];
    if (!segs.length) return null;
    const pts = [];
    segs.forEach((s) => {
      const st = kid(s, "start"), en = kid(s, "end");
      if (!st || !en) return;
      const a = pt(st), b = pt(en);
      const last = pts[pts.length - 1];
      if (!last || last[0] !== a[0] || last[1] !== a[1]) pts.push(a);
      const b1 = kid(s, "basePoint1"), b2 = kid(s, "basePoint2");
      if (/cubic/i.test(attr(s, "type") || "") && (b1 || b2)) pts.push(...bezier(a, b1 ? pt(b1) : a, b2 ? pt(b2) : b, b, 6));
      pts.push(b);
    });
    return { pts, bends: Math.max(0, segs.length - 1) };
  }
  function parseSBML(doc, name) {
    const els = allEls(doc);
    const layout = els.find((e) => ln(e) === "layout");
    if (!layout) throw new Error("this SBML file has no drawing: it needs the Layout extension or CellDesigner coordinates");
    const species = new Map(), comps = new Map();
    for (const e of els) {
      if (ln(e) === "species" && attr(e, "id")) species.set(attr(e, "id"), e);
      if (ln(e) === "compartment" && attr(e, "id")) comps.set(attr(e, "id"), e);
    }
    const m = new Model(name, "SBML");
    const lay = descendants(layout);
    const box = (el) => {
      const bb = kid(el, "boundingBox");
      if (!bb) return null;
      const p = kid(bb, "position"), d = kid(bb, "dimensions");
      return { x: parseFloat(attr(p, "x") || 0), y: parseFloat(attr(p, "y") || 0), w: parseFloat(attr(d, "width") || 0), h: parseFloat(attr(d, "height") || 0) };
    };
    for (const cg of lay.filter((e) => ln(e) === "compartmentGlyph")) {
      const b = box(cg); if (!b) continue;
      const c = comps.get(attr(cg, "compartment"));
      m.addGlyph({ id: attr(cg, "id"), cls: "compartment", ...b, label: c ? attr(c, "name") || attr(c, "id") : "", clone: false });
    }
    const perSpecies = new Map();
    for (const sg of lay.filter((e) => ln(e) === "speciesGlyph")) {
      const b = box(sg); if (!b) continue;
      const sid = attr(sg, "species"), sp = species.get(sid);
      const sbo = sp ? attr(sp, "sboTerm") || "" : "";
      const typ = sp ? attr(sp, "speciesType") || "" : "";
      const macro = MACRO_SBO.has(sbo) || /protein|complex|gene|rna|receptor|channel/i.test(typ);
      m.addGlyph({ id: attr(sg, "id"), cls: macro ? "macromolecule" : "simple chemical", ...b,
        label: sp ? attr(sp, "name") || sid : sid || "", clone: false, species: sid });
      perSpecies.set(sid, (perSpecies.get(sid) || 0) + 1);
    }
    for (const g of m.glyphs.values()) if (g.species && perSpecies.get(g.species) > 1) g.clone = true;
    // complex members are separate glyphs drawn inside the complex box: mark them as nested
    const sp = [...m.glyphs.values()].filter((g) => g.species);
    for (const g of sp) {
      const host = sp.find((c) => c !== g && c.w * c.h > g.w * g.h && g.x >= c.x && g.y >= c.y && g.x + g.w <= c.x + c.w && g.y + g.h <= c.y + c.h);
      if (host) g.nestedIn = host.id;
    }
    const ROLE = { substrate: "consumption", sidesubstrate: "consumption", product: "production", sideproduct: "production",
      modifier: "modulation", activator: "stimulation", inhibitor: "inhibition" };
    for (const rg of lay.filter((e) => ln(e) === "reactionGlyph")) {
      const rid = attr(rg, "id");
      const srgs = descendants(rg).filter((e) => ln(e) === "speciesReferenceGlyph");
      let c = null;
      const b = box(rg);
      const cv = sbmlCurve(kid(rg, "curve"));
      if (b && (b.w || b.h)) c = [b.x + b.w / 2, b.y + b.h / 2];
      else if (cv && cv.pts.length) c = [(cv.pts[0][0] + cv.pts[cv.pts.length - 1][0]) / 2, (cv.pts[0][1] + cv.pts[cv.pts.length - 1][1]) / 2];
      else {
        const st = srgs.map((s) => sbmlCurve(kid(s, "curve"))).filter(Boolean).map((x) => x.pts[0]);
        if (st.length) c = [st.reduce((s, p) => s + p[0], 0) / st.length, st.reduce((s, p) => s + p[1], 0) / st.length];
      }
      if (!c) continue;
      const r = attr(rg, "reaction");
      m.addGlyph({ id: rid, cls: "process", x: c[0] - 5, y: c[1] - 5, w: 10, h: 10, label: r || "", clone: false });
      for (const s of srgs) {
        const gid = attr(s, "speciesGlyph");
        const cls = ROLE[attr(s, "role") || ""];
        if (!cls || !m.glyphs.has(gid)) continue;
        const g = m.glyphs.get(gid);
        const cu = sbmlCurve(kid(s, "curve"));
        const geo = cu && cu.pts.length >= 2 ? cu : { pts: [c, [g.cx, g.cy]], bends: 0 };
        const out = cls === "production";
        m.arcs.push({ id: attr(s, "id"), cls, src: out ? rid : gid, tgt: out ? gid : rid, pts: geo.pts, bends: geo.bends,
          side: /^side/.test(attr(s, "role") || "") });
      }
    }
    for (const tg of lay.filter((e) => ln(e) === "textGlyph")) {      // label geometry: text glyph boxes
      const g = m.glyphs.get(attr(tg, "graphicalObject")); const b = box(tg);
      if (g && b && !g.labelBox) g.labelBox = b;
    }
    return m;
  }


  // ------------------------------------------------------------- CellDesigner (SBML Level 2 + CellDesigner annotations)
  const CD_WORDS = { minus: "-", plus: "+", slash: "/", space: " ", underscore: "_", alpha: "α", beta: "β", gamma: "γ",
    delta: "δ", epsilon: "ε", kappa: "κ", lambda: "λ", mu: "μ", omega: "ω", sub: "", endsub: "", super: "", endsuper: "",
    br: " ", colon: ":", comma: ",", dot: ".", apos: "'", minus2: "-" };
  const cdName = (s) => (s || "").replace(/_([a-z]+\d?)_/gi, (m, w) => (w.toLowerCase() in CD_WORDS ? CD_WORDS[w.toLowerCase()] : m));
  const CD_METABOLITE = new Set(["SIMPLE_MOLECULE", "ION", "DRUG", "UNKNOWN"]);
  function parseCellDesigner(doc, name) {
    const els = allEls(doc);
    const info = new Map(); // species id -> {name, cls}
    for (const e of els) {
      if (ln(e) !== "species" || !attr(e, "id")) continue;
      const clsEl = descendants(e).find((d) => ln(d) === "class");
      const cls = clsEl ? (clsEl.textContent || "").trim() : "";
      if (!info.has(attr(e, "id")) || cls) info.set(attr(e, "id"), { name: cdName(attr(e, "name") || attr(e, "id")), cls });
    }
    const comps = new Map();
    for (const e of els) if (ln(e) === "compartment" && attr(e, "id")) comps.set(attr(e, "id"), cdName(attr(e, "name") || attr(e, "id")));
    const m = new Model(name, "CellDesigner");
    const bounds = (el) => {
      const b = kid(el, "bounds");
      return b ? { x: parseFloat(attr(b, "x") || 0), y: parseFloat(attr(b, "y") || 0), w: parseFloat(attr(b, "w") || 0), h: parseFloat(attr(b, "h") || 0) } : null;
    };
    for (const e of els.filter((x) => ln(x) === "compartmentAlias")) {
      const b = bounds(e); if (!b || !(b.w || b.h)) continue;
      m.addGlyph({ id: attr(e, "id"), cls: "compartment", ...b, label: comps.get(attr(e, "compartment")) || "", clone: false });
    }
    const perSpecies = new Map();
    for (const e of els.filter((x) => ln(x) === "speciesAlias" || ln(x) === "complexSpeciesAlias")) {
      const b = bounds(e); if (!b) continue;
      const sid = attr(e, "species"), sp = info.get(sid) || { name: sid, cls: "" };
      const cls = ln(e) === "complexSpeciesAlias" ? "complex" : CD_METABOLITE.has(sp.cls) || !sp.cls ? "simple chemical" : sp.cls === "DEGRADED" ? "source and sink" : "macromolecule";
      m.addGlyph({ id: attr(e, "id"), cls, ...b, label: sp.name, clone: false, species: sid, nestedIn: attr(e, "complexSpeciesAlias") || undefined });
      if (!attr(e, "complexSpeciesAlias")) perSpecies.set(sid, (perSpecies.get(sid) || 0) + 1);
    }
    for (const g of m.glyphs.values()) if (g.species && !g.nestedIn && perSpecies.get(g.species) > 1) g.clone = true;
    const MOD = { CATALYSIS: "catalysis", UNKNOWN_CATALYSIS: "catalysis", INHIBITION: "inhibition", UNKNOWN_INHIBITION: "inhibition",
      PHYSICAL_STIMULATION: "stimulation", TRIGGER: "necessary stimulation", MODULATION: "modulation" };
    for (const r of els.filter((x) => ln(x) === "reaction")) {
      const rid = attr(r, "id");
      const d = descendants(r);
      const aliasesOf = (tag) => d.filter((x) => ln(x) === tag).map((x) => attr(x, "alias")).filter((a) => m.glyphs.has(a));
      const bre = aliasesOf("baseReactant"), bpr = aliasesOf("baseProduct");
      if (!bre.length && !bpr.length) continue;
      const cen = (ids) => [mean(ids.map((i) => m.glyphs.get(i).cx)), mean(ids.map((i) => m.glyphs.get(i).cy))];
      const ends = [bre.length ? cen(bre) : null, bpr.length ? cen(bpr) : null].filter(Boolean);
      const c = ends.length === 2 ? [(ends[0][0] + ends[1][0]) / 2, (ends[0][1] + ends[1][1]) / 2] : ends[0];
      m.addGlyph({ id: rid, cls: "process", x: c[0] - 5, y: c[1] - 5, w: 10, h: 10, label: rid, clone: false });
      const line = (gid) => [[m.glyphs.get(gid).cx, m.glyphs.get(gid).cy], c];
      let k = 0;
      for (const a of bre) m.arcs.push({ id: rid + ":" + k++, cls: "consumption", src: a, tgt: rid, pts: line(a), bends: 0 });
      for (const a of aliasesOf("reactantLink")) m.arcs.push({ id: rid + ":" + k++, cls: "consumption", src: a, tgt: rid, pts: line(a), bends: 0, side: true });
      for (const a of bpr) m.arcs.push({ id: rid + ":" + k++, cls: "production", src: rid, tgt: a, pts: line(a).reverse(), bends: 0 });
      for (const a of aliasesOf("productLink")) m.arcs.push({ id: rid + ":" + k++, cls: "production", src: rid, tgt: a, pts: line(a).reverse(), bends: 0, side: true });
      for (const mod of d.filter((x) => ln(x) === "modification")) {
        const cls = MOD[attr(mod, "type")] || "modulation";
        for (const a of (attr(mod, "aliases") || "").split(",").map((x) => x.trim()).filter((x) => m.glyphs.has(x)))
          m.arcs.push({ id: rid + ":" + k++, cls, src: a, tgt: rid, pts: line(a), bends: 0 });
      }
    }
    m.note = "CellDesigner lines are read as straight lines from each node to its reaction, so bend points are not counted.";
    return m;
  }

  // ------------------------------------------------------------- Escher JSON
  function parseEscher(data, name) {
    const body = Array.isArray(data) ? data[1] : data;
    if (!body || !body.nodes || !body.reactions) throw new Error("not an Escher map (needs nodes and reactions)");
    const m = new Model(name, "Escher");
    const nodes = body.nodes;
    const perId = new Map();
    for (const [id, n] of Object.entries(nodes)) {
      if (n.node_type !== "metabolite") continue;
      const s = n.node_is_primary ? 30 : 20;
      const base = (n.bigg_id || "").replace(/_[a-z0-9]{1,2}$/i, "");
      m.addGlyph({ id, cls: "simple chemical", x: n.x - s / 2, y: n.y - s / 2, w: s, h: s,
        label: n.name || n.bigg_id || id, keys: [normName(base), normName(n.bigg_id || "")], clone: false, bigg: n.bigg_id,
        primary: n.node_is_primary !== false,
        labelBox: n.label_x !== undefined ? { x: n.label_x, y: n.label_y, w: 0, h: 0 } : null });
      perId.set(n.bigg_id, (perId.get(n.bigg_id) || 0) + 1);
    }
    for (const g of m.glyphs.values()) if (g.bigg && perId.get(g.bigg) > 1) g.clone = true;
    for (const [rid, r] of Object.entries(body.reactions)) {
      const segs = Object.entries(r.segments || {});
      const markers = new Set();
      for (const [, s] of segs) for (const k of [s.from_node_id, s.to_node_id]) {
        const n = nodes[k]; if (n && n.node_type !== "metabolite") markers.add(String(k));
      }
      const mk = [...markers];
      const mid = mk.find((k) => nodes[k].node_type === "midmarker") || mk[0];
      if (mid === undefined) continue;
      const mn = nodes[mid];
      const pid = "rxn:" + mid;
      m.addGlyph({ id: pid, cls: "process", x: mn.x - 5, y: mn.y - 5, w: 10, h: 10, label: r.bigg_id || r.name || rid, clone: false,
        labelBox: r.label_x !== undefined ? { x: r.label_x, y: r.label_y, w: 0, h: 0 } : null });
      const coef = new Map((r.metabolites || []).map((x) => [x.bigg_id, x.coefficient]));
      for (const [sid, s] of segs) {
        const a = nodes[s.from_node_id], b = nodes[s.to_node_id];
        if (!a || !b) continue;
        const aMet = a.node_type === "metabolite", bMet = b.node_type === "metabolite";
        if (aMet === bMet) continue;
        const metId = String(aMet ? s.from_node_id : s.to_node_id);
        const met = aMet ? a : b;
        const c = coef.get(met.bigg_id);
        if (c === undefined || c === 0) continue;
        let p0 = [a.x, a.y], p1 = [b.x, b.y];
        const c1 = s.b1 ? [s.b1.x, s.b1.y] : null, c2 = s.b2 ? [s.b2.x, s.b2.y] : null;
        let pts = [p0, ...(c1 || c2 ? bezier(p0, c1 || p0, c2 || p1, p1, 6) : []), p1];
        if (!aMet) pts = pts.reverse(); // always metabolite -> marker
        const cls = c < 0 ? "consumption" : "production";
        m.arcs.push({ id: rid + ":" + sid, cls, src: cls === "consumption" ? metId : pid, tgt: cls === "consumption" ? pid : metId, pts, bends: 0,
          side: met.node_is_primary === false });
      }
    }
    return m;
  }

  function parseText(text, name, DOMParserImpl) {
    const t = text.trimStart();
    if (t[0] === "[" || t[0] === "{") return parseEscher(JSON.parse(text), name);
    const P = DOMParserImpl || root.DOMParser;
    const doc = new P().parseFromString(text, "application/xml");
    const err = doc.getElementsByTagName("parsererror");
    if (err && err.length) throw new Error("not a layout file: expected SBGN-ML, SBML with Layout, or Escher JSON");
    const r = doc.documentElement;
    const tag = ln(r);
    if (tag === "sbgn") return parseSBGN(doc, name);
    if (tag === "sbml") {
      const els = allEls(doc);
      if (els.some((e) => ln(e) === "layout")) return parseSBML(doc, name);
      if (els.some((e) => ln(e) === "speciesAlias")) return parseCellDesigner(doc, name);
      throw new Error("this SBML file has no drawing: it needs the Layout extension or CellDesigner coordinates");
    }
    if (allEls(doc).some((e) => ln(e) === "map")) return parseSBGN(doc, name);
    throw new Error("not a layout file: expected SBGN-ML, SBML or Escher JSON, found <" + tag + ">");
  }

  // ------------------------------------------------------------- numeric helpers
  const mean = (a) => a.reduce((s, x) => s + x, 0) / a.length;
  function median(a) {
    const s = [...a].sort((x, y) => x - y), n = s.length;
    return n % 2 ? s[(n - 1) / 2] : (s[n / 2 - 1] + s[n / 2]) / 2;
  }
  const std = (a) => { const mu = mean(a); return Math.sqrt(mean(a.map((x) => (x - mu) ** 2))); };
  const dist = (p, q) => Math.hypot(p[0] - q[0], p[1] - q[1]);
  const arcLen = (pts) => { let s = 0; for (let i = 0; i < pts.length - 1; i++) s += dist(pts[i], pts[i + 1]); return s; };
  const r3 = (x) => (x === null || Number.isNaN(x) ? NaN : Math.round(x * 1000) / 1000);
  const r4 = (x) => (Number.isNaN(x) ? NaN : Math.round(x * 10000) / 10000);
  function roundHalfEven(x) {
    const f = Math.floor(x), d = x - f;
    if (d > 0.5) return f + 1;
    if (d < 0.5) return f;
    return f % 2 === 0 ? f : f + 1;
  }
  const pmod = (a, n) => ((a % n) + n) % n;
  function angleBetween(v1, v2) {
    const n1 = Math.hypot(v1[0], v1[1]), n2 = Math.hypot(v2[0], v2[1]);
    if (n1 === 0 || n2 === 0) return NaN;
    const c = Math.min(1, Math.max(-1, (v1[0] * v2[0] + v1[1] * v2[1]) / (n1 * n2)));
    return (Math.acos(c) * 180) / Math.PI;
  }
  function fitCircle(xy) {
    // least squares for [2x 2y 1] p = x^2 + y^2 via normal equations
    const M = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], v = [0, 0, 0];
    for (const [x, y] of xy) {
      const row = [2 * x, 2 * y, 1], b = x * x + y * y;
      for (let i = 0; i < 3; i++) { v[i] += row[i] * b; for (let j = 0; j < 3; j++) M[i][j] += row[i] * row[j]; }
    }
    const A = M.map((r, i) => [...r, v[i]]);
    for (let i = 0; i < 3; i++) {
      let p = i;
      for (let k = i + 1; k < 3; k++) if (Math.abs(A[k][i]) > Math.abs(A[p][i])) p = k;
      [A[i], A[p]] = [A[p], A[i]];
      if (Math.abs(A[i][i]) < 1e-12) return null;
      for (let k = 0; k < 3; k++) if (k !== i) {
        const f = A[k][i] / A[i][i];
        for (let j = i; j < 4; j++) A[k][j] -= f * A[i][j];
      }
    }
    const cx = A[0][3] / A[0][0], cy = A[1][3] / A[1][1], c = A[2][3] / A[2][2];
    return { cx, cy, r: Math.sqrt(Math.max(c + cx * cx + cy * cy, 1e-9)) };
  }
  function orient(a, b, c) { return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]); }
  function segX(p1, p2, p3, p4) {
    const d1 = orient(p3, p4, p1), d2 = orient(p3, p4, p2), d3 = orient(p1, p2, p3), d4 = orient(p1, p2, p4);
    if (!(d1 * d2 < 0 && d3 * d4 < 0)) return null;
    const t = d1 / (d1 - d2);
    return [p1[0] + t * (p2[0] - p1[0]), p1[1] + t * (p2[1] - p1[1])];
  }

  // ------------------------------------------------------------- metrics
  function isCurrency(g, cur) {
    if (cur.has(normName(g.label))) return true;
    return g.keys.some((k) => k && cur.has(k));
  }

  // Every crossing between lines that share no end node, counted on all segments (no size cap). A uniform grid keeps
  // this fast on large maps: only segments whose bounding boxes share a cell are compared, and each crossing is counted
  // once, in the cell that contains the crossing point - so the result equals the all-pairs count.
  function crossings(m) {
    const segs = [];
    m.arcs.forEach((a, i) => { for (let j = 0; j < a.pts.length - 1; j++) segs.push([i, a.pts[j], a.pts[j + 1]]); });
    const ends = m.arcs.map((a) => [a.src, a.tgt]);
    const share = (x, y) => x[0] === y[0] || x[0] === y[1] || x[1] === y[0] || x[1] === y[1];
    if (segs.length < 2) return { n: 0, pts: [] };
    let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity, tot = 0;
    for (const [, p, q] of segs) {
      x0 = Math.min(x0, p[0], q[0]); y0 = Math.min(y0, p[1], q[1]); x1 = Math.max(x1, p[0], q[0]); y1 = Math.max(y1, p[1], q[1]);
      tot += Math.abs(p[0] - q[0]) + Math.abs(p[1] - q[1]);
    }
    // cell size: about the mean segment extent, but never more than ~4 cells per segment on average
    const span = Math.max(x1 - x0, y1 - y0, 1e-9);
    const cell = Math.max(tot / segs.length, span / Math.ceil(Math.sqrt(4 * segs.length)), 1e-9);
    const nx = Math.floor((x1 - x0) / cell) + 1;
    const key = (cx, cy) => cy * nx + cx;
    const ci = (v, o) => Math.floor((v - o) / cell);
    const grid = new Map();
    segs.forEach(([, p, q], k) => {
      const ax = ci(Math.min(p[0], q[0]), x0), bx = ci(Math.max(p[0], q[0]), x0);
      const ay = ci(Math.min(p[1], q[1]), y0), by = ci(Math.max(p[1], q[1]), y0);
      for (let cx = ax; cx <= bx; cx++) for (let cy = ay; cy <= by; cy++) {
        const kk = key(cx, cy); if (!grid.has(kk)) grid.set(kk, []); grid.get(kk).push(k);
      }
    });
    let n = 0; const pts = [];
    for (const [kk, list] of grid) {
      for (let u = 0; u < list.length; u++) {
        const [ai, p1, p2] = segs[list[u]];
        for (let v = u + 1; v < list.length; v++) {
          const [aj, p3, p4] = segs[list[v]];
          if (ai === aj || share(ends[ai], ends[aj])) continue;
          const x = segX(p1, p2, p3, p4);
          if (x && key(ci(x[0], x0), ci(x[1], y0)) === kk) { n++; pts.push(x); }   // count it in its own cell only
        }
      }
    }
    return { n, pts };
  }

  function overlaps(m) {
    const gs = m.nodes(null, SKIP_AREA);
    // nested glyphs (complex members) never count as overlapping their own complex or each other
    const top = (g) => { let seen = 0; while (g.nestedIn && m.glyphs.has(g.nestedIn) && seen++ < 20) g = m.glyphs.get(g.nestedIn); return g.id; };
    const related = (a, b) => top(a) === top(b);
    let n = 0; const ids = new Set();
    for (let i = 0; i < gs.length; i++) {
      const a = gs[i];
      for (let j = i + 1; j < gs.length; j++) {
        const b = gs[j];
        if ((a.nestedIn || b.nestedIn) && related(a, b)) continue;
        if (a.x < b.x + b.w && b.x < a.x + a.w && a.y < b.y + b.h && b.y < a.y + a.h) { n++; ids.add(a.id); ids.add(b.id); }
      }
    }
    return { n, ids };
  }

  function participants(m) {
    const by = new Map();
    const get = (k) => { if (!by.has(k)) by.set(k, { in: [], out: [], mod: [] }); return by.get(k); };
    for (const a of m.arcs) {
      if (INPUT.has(a.cls) && m.glyphs.has(a.tgt)) get(a.tgt).in.push(a.src);
      else if (OUTPUT.has(a.cls)) get(a.src).out.push(a.tgt);
      else if (MODIFIER.has(a.cls)) get(a.tgt).mod.push(a.src);
    }
    const opp = [], perp = [], dirs = [];
    const centroid = (ids) => [mean(ids.map((s) => m.glyphs.get(s).cx)), mean(ids.map((s) => m.glyphs.get(s).cy))];
    for (const [pid, d] of by) {
      const p = m.glyphs.get(pid);
      if (!p || !PROCESS.has(p.cls) || !d.in.length || !d.out.length) continue;
      const ci = centroid(d.in), co = centroid(d.out);
      const vin = [ci[0] - p.cx, ci[1] - p.cy], vout = [co[0] - p.cx, co[1] - p.cy];
      const ang = angleBetween(vin, vout);
      if (!Number.isNaN(ang)) opp.push(ang >= 120);
      const flow = [vout[0] - vin[0], vout[1] - vin[1]];
      if (Math.hypot(flow[0], flow[1]) > 0) {
        dirs.push(pmod(roundHalfEven(Math.atan2(flow[1], flow[0]) / (Math.PI / 2)), 4));
        for (const s of d.mod) {
          const g = m.glyphs.get(s);
          const a2 = angleBetween([g.cx - p.cx, g.cy - p.cy], flow);
          if (!Number.isNaN(a2)) perp.push(a2 >= 60 && a2 <= 120);
        }
      }
    }
    const bc = [0, 0, 0, 0]; dirs.forEach((d) => bc[d]++);
    const mx = Math.max(...bc), am = bc.indexOf(mx);
    return {
      reactions_scored: opp.length,
      inout_opposite_share: opp.length ? mean(opp.map(Number)) : NaN,
      modifier_perpendicular_share: perp.length ? mean(perp.map(Number)) : NaN,
      flow_main_dir: dirs.length ? ["right", "down", "left", "up"][am] : "",
      flow_main_dir_share: dirs.length ? mx / dirs.length : NaN,
    };
  }

  // networkx-compatible undirected graph (insertion-ordered adjacency)
  function metaboliteGraph(m, cur) {
    const ins = new Map(), outs = new Map();
    const push = (mp, k, v) => { if (!mp.has(k)) mp.set(k, []); mp.get(k).push(v); };
    for (const a of m.arcs) {
      if (INPUT.has(a.cls)) push(ins, a.tgt, a.src);
      else if (OUTPUT.has(a.cls)) push(outs, a.src, a.tgt);
    }
    const G = new Map();
    const addNode = (n) => { if (!G.has(n)) G.set(n, new Map()); };
    const pids = [...new Set([...ins.keys(), ...outs.keys()])].sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
    for (const pid of pids)
      for (const s of ins.get(pid) || [])
        for (const t of outs.get(pid) || []) {
          if (isCurrency(m.glyphs.get(s), cur) || isCurrency(m.glyphs.get(t), cur)) continue;
          addNode(s); addNode(t);
          G.get(s).set(t, true); G.get(t).set(s, true);
        }
    return G;
  }
  function cycleBasis(G) {
    const gnodes = new Map([...G.keys()].map((k) => [k, true]));
    const cycles = [];
    while (gnodes.size) {
      const keys = [...gnodes.keys()];
      const root = keys[keys.length - 1];
      gnodes.delete(root);
      const stack = [root], pred = new Map([[root, root]]), used = new Map([[root, new Set()]]);
      while (stack.length) {
        const z = stack.pop(), zused = used.get(z);
        for (const nbr of G.get(z).keys()) {
          if (!used.has(nbr)) { pred.set(nbr, z); stack.push(nbr); used.set(nbr, new Set([z])); }
          else if (nbr === z) cycles.push([z]);
          else if (!zused.has(nbr)) {
            const pn = used.get(nbr); const cyc = [nbr, z]; let p = pred.get(z);
            while (!pn.has(p)) { cyc.push(p); p = pred.get(p); }
            cyc.push(p); cycles.push(cyc); pn.add(z);
          }
        }
      }
      for (const n of pred.keys()) gnodes.delete(n);
    }
    return cycles;
  }
  // Simple cycles of an undirected graph (each once), length 3..maxLen, at most cap of them.
  function simpleCycles(G, maxLen, cap) {
    const order = [...G.keys()], idx = new Map(order.map((n, i) => [n, i])), out = [];
    for (const s0 of order) {
      const i0 = idx.get(s0), path = [s0], onPath = new Set([s0]);
      const dfs = (v) => {
        if (out.length >= cap) return;
        for (const w of G.get(v).keys()) {
          if (idx.get(w) < i0) continue;
          if (w === s0 && path.length >= 3) { if (idx.get(path[1]) < idx.get(path[path.length - 1])) out.push([...path]); continue; }
          if (onPath.has(w) || path.length >= maxLen) continue;
          path.push(w); onPath.add(w); dfs(w); path.pop(); onPath.delete(w);
        }
      };
      dfs(s0);
      if (out.length >= cap) break;
    }
    return out;
  }
  // Rule "cycles look like cycles", as defined in sbml_pathway_layout.py / svg_layout_metrics.py: a cycle passes when its
  // metabolites and reaction nodes, in cycle order, form a closed loop - a polygon that does not cross itself, winds once
  // around its centre in one direction and does not fold inwards (circle, stadium, rounded rectangle... all qualify).
  function loopShape(m, cyc, step) {
    const c = (id) => [m.glyphs.get(id).cx, m.glyphs.get(id).cy];
    const g = [mean(cyc.map((n) => c(n)[0])), mean(cyc.map((n) => c(n)[1]))];
    const pts = [];
    for (let k = 0; k < cyc.length; k++) {
      const a = cyc[k], b = cyc[(k + 1) % cyc.length];
      pts.push(c(a));
      const rs = step.get(a < b ? a + "\u0000" + b : b + "\u0000" + a) || [];
      const want = (dist(c(a), g) + dist(c(b), g)) / 2;      // the step's reaction node ON the loop
      let best = null;
      for (const r of rs) { const d = Math.abs(dist(c(r), g) - want); if (!best || d < best[0]) best = [d, r]; }
      if (best) pts.push(c(best[1]));
    }
    const cen = [mean(pts.map((p) => p[0])), mean(pts.map((p) => p[1]))];
    const th = []; let prev = null, off = 0;
    for (const p of [...pts, pts[0]]) {
      let t = Math.atan2(p[1] - cen[1], p[0] - cen[0]);
      if (prev !== null) { while (t + off - prev > Math.PI) off -= 2 * Math.PI; while (t + off - prev < -Math.PI) off += 2 * Math.PI; }
      th.push(t + off); prev = t + off;
    }
    const d = th.slice(1).map((t, k) => t - th[k]);
    const orderOk = (d.every((x) => x > -1e-9) || d.every((x) => x < 1e-9)) && Math.abs(Math.abs(th[th.length - 1] - th[0]) - 2 * Math.PI) < 0.5;
    let simple = true;
    const E = pts.map((p, k) => [p, pts[(k + 1) % pts.length]]);
    for (let x = 0; x < E.length && simple; x++)
      for (let y = x + 2; y < E.length; y++) {
        if (x === 0 && y === E.length - 1) continue;
        if (segX(E[x][0], E[x][1], E[y][0], E[y][1])) { simple = false; break; }
      }
    const turns = pts.map((b, k) => {
      const a = pts[(k - 1 + pts.length) % pts.length], cc = pts[(k + 1) % pts.length];
      const cr = (b[0] - a[0]) * (cc[1] - b[1]) - (b[1] - a[1]) * (cc[0] - b[0]);
      return cr / ((dist(a, b) * dist(b, cc)) || 1);
    });
    const sign = Math.sign(turns.reduce((s, t) => s + t, 0)) || 1;
    const folds = turns.filter((t) => t * sign < -0.35).length;
    const fc = fitCircle(pts);
    const circ = fc ? Math.max(0, 1 - std(pts.map(([x, y]) => Math.hypot(x - fc.cx, y - fc.cy))) / fc.r) : NaN;
    return { ok: orderOk && simple && folds === 0, orderOk, simple, folds, circ, pts, circle: fc,
      dir: (th[th.length - 1] - th[0]) > 0 ? "ccw" : "cw" };   // SVG/SBML y axis points down
  }
  function cycleMetrics(m, cur, cycleNames) {
    // main steps only: substrate -> product of the same reaction, no secondary lines, no cofactors
    const ins = new Map(), outs = new Map();
    const push = (mp, k, v) => { if (!mp.has(k)) mp.set(k, []); mp.get(k).push(v); };
    for (const a of m.arcs) {
      if (a.side) continue;
      if (INPUT.has(a.cls)) push(ins, a.tgt, a.src); else if (OUTPUT.has(a.cls)) push(outs, a.src, a.tgt);
    }
    const want = cycleNames && cycleNames.length ? new Set(cycleNames.map(normName)) : null;
    const okNode = (n) => m.glyphs.has(n) && !isCurrency(m.glyphs.get(n), cur) && (!want || want.has(normName(m.glyphs.get(n).label)));
    const G = new Map(), step = new Map();
    const addE = (s, t, p) => {
      if (!G.has(s)) G.set(s, new Map()); if (!G.has(t)) G.set(t, new Map());
      G.get(s).set(t, true); G.get(t).set(s, true);
      const k = s < t ? s + "\u0000" + t : t + "\u0000" + s; if (!step.has(k)) step.set(k, []); step.get(k).push(p);
    };
    const pids = [...new Set([...ins.keys(), ...outs.keys()])].sort();
    for (const p of pids) for (const s of ins.get(p) || []) for (const t of outs.get(p) || [])
      if (s !== t && okNode(s) && okNode(t)) addE(s, t, p);
    const cand = G.size ? simpleCycles(G, 16, 20000) : [];
    const shp = new Map(cand.map((c) => [c, loopShape(m, c, step)]));
    const main = []; let pool = cand.slice(); const used = new Set();
    while (pool.length) {                                         // main cycles: the longest drawn as a closed loop (if at
      const mx = Math.max(...pool.map((c) => c.length));          // least half as long as the longest loop), else the longest
      const rnd = pool.filter((c) => shp.get(c).ok && c.length >= 0.5 * mx);
      const pick = rnd.length ? rnd.reduce((b, c) => (c.length > b.length || (c.length === b.length && shp.get(c).circ > shp.get(b).circ) ? c : b))
                              : pool.filter((c) => c.length === mx).reduce((b, c) => (shp.get(c).circ > shp.get(b).circ ? c : b));
      main.push(pick); pick.forEach((n) => used.add(n));
      pool = pool.filter((c) => !c.some((n) => used.has(n)));
    }
    const empty = { cycle_len: 0, cycles_main: 0, cycles_closed_loop: 0, cycle_closed_loop_share: NaN, cycle_circularity: NaN,
      cycle_order_consistency: NaN, cycle_direction: "", cycle_members: "", cycle_alternative_routes: 0, _cycle: null };
    if (!main.length) return empty;
    const first = main[0], s0 = shp.get(first), nOk = main.filter((c) => shp.get(c).ok).length;
    return {
      cycle_len: first.length, cycles_main: main.length, cycles_closed_loop: nOk, cycle_closed_loop_share: r3(nOk / main.length),
      cycle_circularity: r3(s0.circ), cycle_order_consistency: s0.orderOk ? 1 : 0, cycle_direction: s0.dir,
      cycle_members: first.map((n) => m.glyphs.get(n).label).join("; "), cycle_alternative_routes: cand.length - main.length,
      _cycle: { ids: first, circle: s0.circle, loop: s0.pts, ok: s0.ok },
    };
  }

  // ------------------------------------------------------------- rules that come from the Python tools
  // Liang-Barsky: does segment p-q pass through box b (grown by m)?
  function segHitsBox(p, q, b, m) {
    const x0 = b.x - m, y0 = b.y - m, x1 = b.x + b.w + m, y1 = b.y + b.h + m, dx = q[0] - p[0], dy = q[1] - p[1];
    let t0 = 0, t1 = 1;
    for (const [pp, qq] of [[-dx, p[0] - x0], [dx, x1 - p[0]], [-dy, p[1] - y0], [dy, y1 - p[1]]]) {
      if (Math.abs(pp) < 1e-12) { if (qq < 0) return false; }
      else if (pp < 0) t0 = Math.max(t0, qq / pp); else t1 = Math.min(t1, qq / pp);
    }
    return t0 <= t1;
  }
  function linesThroughNodes(m) {
    const gs = m.nodes(null, SKIP_AREA).filter((g) => !g.nestedIn);   // every node, reactions included (as in the Python tools)
    let n = 0; const pts = [];
    for (const a of m.arcs) {
      for (const g of gs) {
        if (g.id === a.src || g.id === a.tgt) continue;
        const shr = { x: g.x + 2, y: g.y + 2, w: g.w - 4, h: g.h - 4 };
        if (shr.w <= 0 || shr.h <= 0) continue;
        for (let k = 0; k < a.pts.length - 1; k++) if (segHitsBox(a.pts[k], a.pts[k + 1], shr, 0)) { n++; pts.push([g.cx, g.cy]); break; }
      }
    }
    return { n, pts };
  }
  // labels at fixed positions per node type: spread of the label's offset from its node, per type (class + clone)
  function labelPositions(m) {
    const by = new Map(); let withGeo = 0;
    for (const g of m.nodes(null, SKIP_AREA)) {
      const lb = g.labelBox; if (!lb) continue;
      withGeo++;
      const k = g.cls + (g.clone ? " (clone)" : "") + (g.primary === false ? " (secondary)" : "");
      const ox = (lb.x + lb.w / 2 - g.cx) / Math.max(g.w, 1), oy = (lb.y + lb.h / 2 - g.cy) / Math.max(g.h, 1);
      if (!by.has(k)) by.set(k, []); by.get(k).push([ox, oy]);
    }
    let worst = 0;
    for (const v of by.values()) if (v.length >= 2) worst = Math.max(worst, std(v.map((p) => p[0])), std(v.map((p) => p[1])));
    const boxes = m.nodes(null, SKIP_AREA).filter((g) => g.labelBox && g.labelBox.w > 0 && g.labelBox.h > 0).map((g) => g.labelBox);
    let ov = 0;
    for (let i = 0; i < boxes.length; i++) for (let j = i + 1; j < boxes.length; j++) {
      const a = boxes[i], b = boxes[j];
      const ix = Math.min(a.x + a.w, b.x + b.w) - Math.max(a.x, b.x), iy = Math.min(a.y + a.h, b.y + b.h) - Math.max(a.y, b.y);
      if (ix > 0 && iy > 0 && ix * iy > 0.15 * Math.min(a.w * a.h, b.w * b.h)) ov++;
    }
    return { labels_with_geometry: withGeo, label_types: by.size, label_offset_spread: withGeo ? r3(worst) : NaN, label_overlaps: ov };
  }
  // reproducible: the same drawing as a reference rendering (same nodes at the same places, same lines)
  function sameLayout(m, ref) {
    if (!ref) return null;
    const near = (a, b) => Math.abs(a - b) < 1e-3;
    let diff = 0;
    for (const g of m.glyphs.values()) {
      const h = ref.glyphs.get(g.id);
      if (!h || !near(g.x, h.x) || !near(g.y, h.y) || !near(g.w, h.w) || !near(g.h, h.h)) diff++;
    }
    diff += Math.abs(ref.glyphs.size - m.glyphs.size);
    const sig = (mm) => mm.arcs.map((a) => a.src + ">" + a.tgt + ":" + a.pts.map((p) => p.map((v) => v.toFixed(2)).join(",")).join(";")).sort().join("|");
    return { nodes_differing: diff, lines_identical: sig(m) === sig(ref) };
  }

  function scoreModel(m, opts = {}) {
    const cur = new Set((opts.currency || CURRENCY_DEFAULT).map(normName));
    const nodes = m.nodes(null, SKIP_AREA);
    const allg = [...m.glyphs.values()];
    const row = { map: m.name, format: m.format, glyphs: allg.length, arcs: m.arcs.length,
      processes: m.nodes(PROCESS).length, metabolites: m.nodes(METABOLITE).length };
    const extra = { crossPts: [], hubIds: new Set(), overlapIds: new Set(), cycle: null, bbox: null };
    if (!allg.length || !m.arcs.length) return { row, extra };
    const x0 = Math.min(...allg.map((g) => g.x)), y0 = Math.min(...allg.map((g) => g.y));
    const x1 = Math.max(...allg.map((g) => g.x + g.w)), y1 = Math.max(...allg.map((g) => g.y + g.h));
    const W = x1 - x0, H = y1 - y0;
    extra.bbox = { x0, y0, x1, y1 };
    const nodeSize = nodes.length ? median(nodes.map((g) => Math.max(g.w, g.h))) : 1;

    const cm = cycleMetrics(m, cur, opts.cycle);
    extra.cycle = cm._cycle; delete cm._cycle;
    Object.assign(row, cm);

    const curG = m.nodes(METABOLITE).filter((g) => isCurrency(g, cur));
    const deg = new Map();
    for (const a of m.arcs) for (const e of [a.src, a.tgt]) {
      if (!deg.has(e)) deg.set(e, new Set());
      deg.get(e).add(e === a.src ? a.tgt : a.src);
    }
    const lens = m.arcs.map((a) => arcLen(a.pts));
    const med = median(lens);
    const curIds = new Set(curG.map((g) => g.id));
    const curArcs = m.arcs.filter((a) => curIds.has(a.src) || curIds.has(a.tgt)).map((a) => arcLen(a.pts));
    const labels = m.nodes(METABOLITE).map((g) => normName(g.label));
    const cnt = new Map(); labels.forEach((l) => cnt.set(l, (cnt.get(l) || 0) + 1));
    const hubs = curG.filter((g) => !g.clone && (deg.get(g.id) || new Set()).size > 1);
    hubs.forEach((g) => extra.hubIds.add(g.id));
    Object.assign(row, {
      currency_glyphs: curG.length,
      currency_hub_violations: hubs.length,
      currency_long_arc_share: curArcs.length ? r3(mean(curArcs.map((l) => (l > 2 * med ? 1 : 0)))) : NaN,
      clone_markers: allg.filter((g) => g.clone).length,
      duplicated_metabolite_labels: [...cnt].filter(([l, c]) => l && c > 1).length,
    });

    const cr = crossings(m);
    extra.crossPts = cr.pts;
    row.edge_crossings = cr.n;
    row.crossings_per_arc = Number.isNaN(cr.n) ? NaN : r3(cr.n / m.arcs.length);

    const pp = participants(m);
    for (const [k, v] of Object.entries(pp)) row[k] = typeof v === "number" && !Number.isInteger(v) ? r3(v) : v;

    row.aspect_ratio = H ? r3(W / H) : NaN;
    const rel = lens.map((l) => l / nodeSize);
    const rmed = median(rel);
    Object.assign(row, {
      edge_len_median: r3(rmed), edge_len_max_rel: r3(Math.max(...rel) / rmed), edge_len_cv: r3(std(rel) / mean(rel)),
    });
    const bends = m.arcs.map((a) => a.bends);
    row.bends_per_arc = r3(mean(bends));
    row.bent_arc_share = r3(mean(bends.map((b) => (b > 0 ? 1 : 0))));
    const area = nodes.reduce((s, g) => s + (g.nestedIn ? 0 : g.w * g.h), 0);
    row.density = W * H ? r4(area / (W * H)) : NaN;
    const ov = overlaps(m);
    row.node_overlaps = ov.n; extra.overlapIds = ov.ids;
    const lt = linesThroughNodes(m); row.lines_through_nodes = lt.n; extra.throughPts = lt.pts;
    Object.assign(row, labelPositions(m));
    const same = sameLayout(m, opts.reference);
    if (same) { row.reference_nodes_differing = same.nodes_differing; row.reference_lines_identical = same.lines_identical; }
    return { row, extra };
  }

  // ------------------------------------------------------------- rules -> scores
  const DEFAULT_THRESHOLDS = {
    reproducible: { pass: 0, fail: 1 },
    labels: { pass: 0.02, fail: 0.2 },
    cycle: { pass: 1, fail: 0.5 },
    through: { pass: 0, fail: 5 },
    clones: { pass: 0, fail: 3 },
    crossings: { pass: 0.02, fail: 0.3 },
    placement: { pass: 0.8, fail: 0.5 },
    aspect: { lo: 0.5, hi: 2, flo: 0.25, fhi: 4 },
    longlines: { pass: 6, fail: 15 },
    bends: { pass: 0.5, fail: 2 },
    density: { lo: 0.05, hi: 0.3, flo: 0.02, fhi: 0.5 },
    overlaps: { pass: 0, fail: 5 },
    weights: { hard: 2, soft: 1 },
  };
  const lowBetter = (v, p, f) => (v <= p ? 1 : v >= f ? 0 : (f - v) / (f - p));
  const highBetter = (v, p, f) => (v >= p ? 1 : v <= f ? 0 : (v - f) / (p - f));
  function band(v, t, log) {
    const T = log ? Math.log : (x) => x;
    if (v >= t.lo && v <= t.hi) return 1;
    if (v < t.lo) return v <= t.flo ? 0 : (T(v) - T(t.flo)) / (T(t.lo) - T(t.flo));
    return v >= t.fhi ? 0 : (T(t.fhi) - T(v)) / (T(t.fhi) - T(t.hi));
  }
  const pct = (x) => Math.round(x * 100) + "%";
  const fmt = (x, d = 2) => (Number.isNaN(x) || x === undefined ? "–" : Number(x).toFixed(d).replace(/\.?0+$/, "") || "0");

  const RULES = [
    { id: "reproducible", kind: "hard", name: "Reproducible rendering",
      eval(r, t) {
        if (r.reference_nodes_differing === undefined) return { na: true, value: "–", detail: "Give a second rendering of the same map (opts.reference) to check it." };
        const v = r.reference_nodes_differing + (r.reference_lines_identical ? 0 : 1);
        return { v, s: lowBetter(v, t.pass, t.fail), value: v ? "differs" : "identical", unit: "",
          detail: `${r.reference_nodes_differing} node(s) placed differently; lines ${r.reference_lines_identical ? "identical" : "differ"} from the reference rendering.` };
      } },
    { id: "labels", kind: "hard", name: "Labels at fixed positions per type",
      eval(r, t) {
        if (!r.labels_with_geometry) return { na: true, value: "–", detail: "This file has no label positions (SBML text glyphs, Escher label_x/label_y or SBGN label boxes)." };
        const v = r.label_offset_spread + (r.label_overlaps ? t.fail : 0);
        return { v, s: lowBetter(v, t.pass, t.fail), value: fmt(r.label_offset_spread, 3), unit: "spread",
          detail: `${r.labels_with_geometry} labels in ${r.label_types} node types: the label's offset from its node varies by at most ${fmt(r.label_offset_spread * 100, 1)}% of the node size within a type; ${r.label_overlaps} overlapping labels. Pass ≤ ${t.pass * 100}%.` };
      } },
    { id: "cycle", kind: "hard", name: "Cycles look like cycles",
      eval(r, t) {
        if (!r.cycles_main) return { na: true, value: "no cycle", detail: "No cycle of 3 or more metabolites (main steps) found. Name one under Settings to check it." };
        const v = r.cycle_closed_loop_share;
        return { v, s: highBetter(v, t.pass, t.fail), value: `${r.cycles_closed_loop}/${r.cycles_main}`, unit: "closed loops",
          detail: `${r.cycles_main} main cycle(s), ${r.cycles_closed_loop} drawn as a closed loop in cycle order (no self-crossing, winds once, no inward fold; circle, stadium or rounded rectangle). Longest: ${r.cycle_len} metabolites, ${r.cycle_direction === "cw" ? "clockwise" : "counter-clockwise"}; ${r.cycle_alternative_routes} alternative route(s) not judged.` };
      } },
    { id: "clones", kind: "hard", name: "Cofactors drawn as clones",
      eval(r, t) {
        const v = r.currency_hub_violations;
        if (r.currency_glyphs === undefined) return { na: true, value: "–", detail: "" };
        return { v, s: lowBetter(v, t.pass, t.fail), value: String(v), unit: v === 1 ? "hub" : "hubs",
          detail: `${r.currency_glyphs} cofactor glyphs; ${v} of them shared by more than one reaction without a clone marker. ${fmt(r.currency_long_arc_share * 100, 0)}% of cofactor lines are over twice the median length.` };
      } },
    { id: "crossings", kind: "hard", name: "No crossings",
      eval(r, t) {
        if (r.edge_crossings === undefined || Number.isNaN(r.edge_crossings)) return { na: true, value: "–", detail: "No lines to check." };
        const v = r.crossings_per_arc;
        return { v, s: lowBetter(v, t.pass, t.fail), value: String(r.edge_crossings), unit: r.edge_crossings === 1 ? "crossing" : "crossings",
          detail: `${fmt(v, 3)} per line across ${r.arcs} lines. Pass ≤ ${t.pass} per line.` };
      } },
    { id: "placement", kind: "hard", name: "Consistent participant placement",
      eval(r, t) {
        if (!r.reactions_scored) return { na: true, value: "–", detail: "No reaction has both a substrate and a product." };
        const v = r.inout_opposite_share;
        const mod = Number.isNaN(r.modifier_perpendicular_share) ? "no modifiers" : `modifiers beside the reaction in ${pct(r.modifier_perpendicular_share)}`;
        return { v, s: highBetter(v, t.pass, t.fail), value: pct(v), unit: "opposite",
          detail: `Substrates and products on opposite sides (≥120°) in ${pct(v)} of ${r.reactions_scored} reactions; ${mod}; ${pct(r.flow_main_dir_share)} run ${r.flow_main_dir}.` };
      } },
    { id: "through", kind: "hard", name: "No line through a node",
      eval(r, t) {
        if (r.lines_through_nodes === undefined) return { na: true, value: "–", detail: "" };
        const v = r.lines_through_nodes;
        return { v, s: lowBetter(v, t.pass, t.fail), value: String(v), unit: v === 1 ? "line" : "lines",
          detail: "Lines that pass through a node they do not connect, reaction nodes included (exact segment/box test, node box shrunk by 2 px)." };
      } },
    { id: "aspect", kind: "hard", name: "Aspect ratio",
      eval(r, t) {
        if (Number.isNaN(r.aspect_ratio) || r.aspect_ratio === undefined) return { na: true, value: "–", detail: "" };
        const v = r.aspect_ratio;
        return { v, s: band(v, t, true), value: fmt(v), unit: "width ÷ height",
          detail: `Target ${t.lo}–${t.hi}. ${v > 1 ? "Wider than tall." : v < 1 ? "Taller than wide." : "Square."}` };
      } },
    { id: "longlines", kind: "soft", name: "Avoid long lines",
      eval(r, t) {
        if (r.edge_len_max_rel === undefined) return { na: true, value: "–", detail: "" };
        const v = r.edge_len_max_rel;
        return { v, s: lowBetter(v, t.pass, t.fail), value: fmt(v, 1) + "×", unit: "longest vs median",
          detail: `Median line is ${fmt(r.edge_len_median, 1)} node widths; spread (CV) ${fmt(r.edge_len_cv)}. Pass ≤ ${t.pass}×.` };
      } },
    { id: "bends", kind: "soft", name: "Avoid bends",
      eval(r, t) {
        if (r.bends_per_arc === undefined) return { na: true, value: "–", detail: "" };
        const v = r.bends_per_arc;
        return { v, s: lowBetter(v, t.pass, t.fail), value: fmt(v), unit: "bends per line",
          detail: `${pct(r.bent_arc_share)} of lines have at least one bend point. Curves in SBML/Escher count as smooth.` };
      } },
    { id: "density", kind: "soft", name: "Density",
      eval(r, t) {
        if (r.density === undefined || Number.isNaN(r.density)) return { na: true, value: "–", detail: "" };
        const v = r.density;
        return { v, s: band(v, t, false), value: pct(v), unit: "canvas covered",
          detail: `Node area ÷ bounding-box area. Target ${pct(t.lo)}–${pct(t.hi)}.` };
      } },
    { id: "overlaps", kind: "soft", name: "No overlapping nodes",
      eval(r, t) {
        if (r.node_overlaps === undefined) return { na: true, value: "–", detail: "" };
        const v = r.node_overlaps;
        return { v, s: lowBetter(v, t.pass, t.fail), value: String(v), unit: v === 1 ? "overlap" : "overlaps",
          detail: `Pairs of node boxes that overlap (compartments excluded).` };
      } },
  ];

  function evaluate(row, th) {
    const t = th || DEFAULT_THRESHOLDS;
    const out = RULES.map((rule) => {
      const e = rule.eval(row, t[rule.id]);
      if (e.na) return { ...rule, ...e, score: null, status: "na" };
      const score = Math.round(e.s * 100);
      return { ...rule, ...e, score, status: e.s >= 0.999 ? "pass" : e.s <= 0.001 ? "fail" : "check" };
    });
    let ws = 0, wt = 0;
    for (const r of out) if (r.score !== null) { const w = t.weights[r.kind]; ws += w * r.score; wt += w; }
    return { rules: out, overall: wt ? Math.round(ws / wt) : null,
      hardPass: out.filter((r) => r.kind === "hard" && r.status === "pass").length,
      hardTotal: out.filter((r) => r.kind === "hard" && r.status !== "na").length };
  }

  const COLUMNS = ["map", "format", "glyphs", "arcs", "processes", "metabolites", "cycle_len", "cycles_main", "cycles_closed_loop",
    "cycle_closed_loop_share", "cycle_alternative_routes", "cycle_circularity", "cycle_order_consistency", "cycle_direction", "cycle_members", "currency_glyphs", "currency_hub_violations",
    "currency_long_arc_share", "clone_markers", "duplicated_metabolite_labels", "edge_crossings", "crossings_per_arc",
    "reactions_scored", "inout_opposite_share", "modifier_perpendicular_share", "flow_main_dir", "flow_main_dir_share",
    "aspect_ratio", "edge_len_median", "edge_len_max_rel", "edge_len_cv", "bends_per_arc", "bent_arc_share",
    "density", "node_overlaps", "lines_through_nodes", "labels_with_geometry", "label_types", "label_offset_spread", "label_overlaps",
    "reference_nodes_differing", "reference_lines_identical"];

  const api = { parseText, parseSBGN, parseSBML, parseCellDesigner, parseEscher, scoreModel, evaluate, RULES, COLUMNS, loopShape, simpleCycles,
    CURRENCY_DEFAULT, DEFAULT_THRESHOLDS, normName, PROCESS, METABOLITE, MODIFIER, crossings };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.LayoutScore = api;
})(typeof window !== "undefined" ? window : globalThis);
