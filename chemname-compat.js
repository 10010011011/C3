/*
 * Chem3D name compatibility layer — isolated, opt-in add-on.
 *
 * Load after Chem3DCore is defined:
 *   <script src="./chemname-compat.js"></script>
 * Resolve through:
 *   Chem3DNameCompat.resolve(input, Chem3DCore)
 *
 * This module does not edit or monkey-patch the existing core. The UI can
 * switch to it in one place once tests are run. Coverage is curated, not a
 * claim of complete IUPAC parsing.
 */
(function (root) {
  'use strict';

  var entries = [
    // Phenols and polyhydric phenols.
    ['phenol', 'Oc1ccccc1'],
    ['naphthol', 'Oc1cccc2ccccc12'],
    ['1-naphthol', 'Oc1cccc2ccccc12'],
    ['2-naphthol', 'Oc1ccc2ccccc2c1'],
    ['catechol', 'Oc1ccccc1O'],
    ['catchol', 'Oc1ccccc1O'],
    ['resorcinol', 'Oc1cccc(O)c1'],
    ['hydroquinone', 'Oc1ccc(O)cc1'],
    ['quinol', 'Oc1ccc(O)cc1'],
    ['pyrogallol', 'Oc1cccc(O)c1O'],
    ['phloroglucinol', 'Oc1cc(O)cc(O)c1'],
    ['cresol', 'Cc1ccc(O)cc1'],
    ['o-cresol', 'Cc1ccccc1O'],
    ['ortho-cresol', 'Cc1ccccc1O'],
    ['m-cresol', 'Cc1cccc(O)c1'],
    ['meta-cresol', 'Cc1cccc(O)c1'],
    ['p-cresol', 'Cc1ccc(O)cc1'],
    ['para-cresol', 'Cc1ccc(O)cc1'],

    // Aromatic hydrocarbons, fused and non-benzenoid systems.
    ['azulene', 'c1ccc2cccc2cc1'],
    ['naphthalene', 'c1ccc2ccccc2c1'],
    ['anthracene', 'c1ccc2cc3ccccc3cc2c1'],
    ['phenanthrene', 'c1ccc2c(c1)ccc1ccccc12'],
    ['biphenyl', 'c1ccccc1-c1ccccc1'],
    ['styrene', 'C=Cc1ccccc1'],

    // Aromatic amines and amides.
    ['aniline', 'Nc1ccccc1'],
    ['acetanilide', 'CC(=O)Nc1ccccc1'],
    ['benzamide', 'NC(=O)c1ccccc1'],
    ['benzylamine', 'NCc1ccccc1'],

    // Aromatic acids and positional isomers.
    ['benzoic acid', 'O=C(O)c1ccccc1'],
    ['benzoicacid', 'O=C(O)c1ccccc1'],
    ['isophthalic acid', 'O=C(O)c1cccc(C(=O)O)c1'],
    ['phthalic acid', 'O=C(O)c1ccccc1C(=O)O'],
    ['terephthalic acid', 'O=C(O)c1ccc(C(=O)O)cc1'],
    ['cinnamic acid', 'O=C(O)/C=C/c1ccccc1'],
    ['trans-cinnamic acid', 'O=C(O)/C=C/c1ccccc1'],
    ['phenylacetic acid', 'O=C(O)Cc1ccccc1'],

    // Anhydrides, esters and carbonyl-rich structures.
    ['benzoic anhydride', 'O=C(OC(=O)c1ccccc1)c1ccccc1'],
    ['benzoicacid anhydride', 'O=C(OC(=O)c1ccccc1)c1ccccc1'],
    ['acetic anhydride', 'CC(=O)OC(=O)C'],
    ['succinic anhydride', 'O=C1OC(=O)CCC1'],
    ['methyl benzoate', 'COC(=O)c1ccccc1'],
    ['ethyl benzoate', 'CCOC(=O)c1ccccc1'],
    ['methyl cinnamate', 'COC(=O)/C=C/c1ccccc1'],

    // Heteroaromatic retained names.
    ['pyridine', 'n1ccccc1'],
    ['pyrrole', '[nH]1cccc1'],
    ['furan', 'o1cccc1'],
    ['thiophene', 's1cccc1'],
    ['imidazole', 'c1ncc[nH]1'],
    ['indole', 'c1ccc2[nH]ccc2c1'],
    ['quinoline', 'c1ccc2ncccc2c1'],
    ['isoquinoline', 'c1ccc2cnccc2c1'],

    // Carbonyl compounds often exposing suffix/prefix conflicts.
    ['formaldehyde', 'C=O'],
    ['acetaldehyde', 'CC=O'],
    ['acetone', 'CC(=O)C'],
    ['acetophenone', 'CC(=O)c1ccccc1'],
    ['benzophenone', 'O=C(c1ccccc1)c1ccccc1'],
    ['benzaldehyde', 'O=Cc1ccccc1'],
    ['benzil', 'O=C(c1ccccc1)C(=O)c1ccccc1']
  ];

  function cleanName(value) {
    return String(value == null ? '' : value)
      .toLowerCase()
      .replace(/[‐‑‒–—]/g, '-')
      .replace(/[\u00a0\t\r\n]+/g, ' ')
      .replace(/\s+/g, ' ')
      .replace(/\s*([,()])\s*/g, '$1')
      .replace(/\s*-\s*/g, '-')
      .replace(/\bbenzonic\b/g, 'benzoic')
      .replace(/\bbenzon(ic|o)\b/g, 'benzoic')
      .replace(/\bcatchol\b/g, 'catechol')
      .trim();
  }

  var synonymMap = Object.create(null);
  entries.forEach(function (entry) {
    synonymMap[cleanName(entry[0])] = entry[1];
  });

  function canonicalSpelling(name) {
    var s = cleanName(name);
    // Accept typography variations such as "3-Ethanoyl oxy -5- Methoxy carbonyl".
    s = s.replace(/\b(\d+)\s*[- ]\s*(?=(?:ethanoyl|hydroxy|methoxy|oxo|methyl|ethyl|propyl|butyl|chloro|bromo|fluoro|iodo)\b)/g, '$1-');
    s = s.replace(/\bethanoyl\s+oxy\b/g, 'acetyloxy');
    s = s.replace(/\bmethoxy\s+carbonyl\b/g, 'methoxycarbonyl');
    s = s.replace(/\bbenzon(ic|o)\s+anhydride\b/g, 'benzoic anhydride');
    return s;
  }

  function toGraph(smiles, core, label) {
    if (!core || typeof core.parseSmiles !== 'function') {
      return { err: 'Chem3DCore.parseSmiles is unavailable; load this module after the core.' };
    }
    var graph = core.parseSmiles(smiles);
    if (!graph || !graph.atoms) {
      return { err: graph && graph.err ? graph.err : 'Could not parse compatibility SMILES for ' + label };
    }
    graph.name = label;
    graph.note = 'Name compatibility layer: ' + label;
    return graph;
  }

  function resolve(raw, core) {
    core = core || root.Chem3DCore;
    var original = String(raw == null ? '' : raw);
    var cleaned = canonicalSpelling(original);
    if (!cleaned) return { err: 'Type a compound name, formula, or SMILES.' };

    var known = synonymMap[cleaned] || synonymMap[cleaned.replace(/\s/g, '')];
    if (known) return toGraph(known, core, cleaned);

    // Delegate systematic names and formulae to the existing implementation.
    if (core && typeof core.fromName === 'function') {
      var named = core.fromName(cleaned);
      if (named && named.src) {
        var graph = /^smiles:/i.test(named.src)
          ? toGraph(named.src.replace(/^smiles:/i, ''), core, named.name || cleaned)
          : (core.resolve ? core.resolve(named.src) : null);
        if (graph && graph.atoms) {
          graph.name = named.name || cleaned;
          graph.note = named.note || 'Interpreted by Chem3D core';
          return graph;
        }
      }
    }
    if (core && typeof core.resolve === 'function') {
      var fallback = core.resolve(original);
      return fallback && fallback.atoms ? fallback :
        (fallback || { err: 'Could not resolve compound name.' });
    }
    return { err: 'Chem3DCore is unavailable.' };
  }

  root.Chem3DNameCompat = {
    resolve: resolve,
    cleanName: cleanName,
    canonicalSpelling: canonicalSpelling,
    knownNames: function () { return Object.keys(synonymMap).sort(); },
    version: '1.0.1'
  };
})(typeof window !== 'undefined' ? window : globalThis);
