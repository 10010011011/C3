/*
 * Chem3D name compatibility layer
 * ===============================
 * Isolated add-on for the existing Chem3DCore name parser.
 *
 * PURPOSE
 * - Normalise common retained names and common spelling/spacing variants.
 * - Provide a curated coverage layer for fused/heteroaromatic rings, phenols,
 *   anilines, polyfunctional aromatics, acid anhydrides and substituted acids.
 * - Make functional-group principal-suffix vs prefix handling explicit.
 * - Never modify Chem3DCore or the application UI.
 *
 * INTEGRATION (in index.html, after Chem3DCore has been defined and before
 * the app script runs):
 *   <script src="./chemname-compat.js"></script>
 *
 * Existing build() may then resolve with:
 *   const G = window.Chem3DNameCompat.resolve(input, window.Chem3DCore);
 * This layer delegates known structures to Core.parseSmiles() and falls back
 * to Core.fromName()/Core.resolve(). It is deliberately opt-in so it cannot
 * silently change existing behavior before it is tested.
 *
 * Scope note: this is not a complete IUPAC parser. For a correct but
 * unsupported name, it returns the core parser's error rather than inventing
 * a graph. Review the test cases before integrating into the default path.
 */
(function (root) {
  'use strict';

  var entries = [
    // Phenols and polyhydric phenols (retained names).
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
    ['azulene', 'c1ccc2cccc- c2c1'.replace(/ /g, '')],
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

    // Aromatic carboxylic acids and positional isomers.
    ['benzoic acid', 'O=C(O)c1ccccc1'],
    ['benzoicacid', 'O=C(O)c1ccccc1'],
    ['isophthalic acid', 'O=C(O)c1cccc(C(=O)O)c1'],
    ['phthalic acid', 'O=C(O)c1ccccc1C(=O)O'],
    ['terephthalic acid', 'O=C(O)c1ccc(C(=O)O)cc1'],
    ['cinnamic acid', 'O=C(O)/C=C/c1ccccc1'],
    ['trans-cinnamic acid', 'O=C(O)/C=C/c1ccccc1'],
    ['phenylacetic acid', 'O=C(O)Cc1ccccc1'],

    // Acid anhydrides, esters, and other carbonyl-rich structures.
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

    // Common carbonyls (saturated and aromatic).
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
      .replace(/\bbenzon(ic|o)\b/g, 'benzoic')
      .replace(/\bbenzonic\b/g, 'benzoic')
      .replace(/\bcat[e]?chol\b/g, 'catechol')
      .replace(/\b1- naphthol\b/g, '1-naphthol')
      .trim();
  }

  var synonymMap = Object.create(null);
  entries.forEach(function (entry) {
    synonymMap[cleanName(entry[0])] = entry[1];
  });

  // Normalise whitespace around hyphens in long names while preserving locants.
  function canonicalSpelling(name) {
    var s = cleanName(name);
    s = s.replace(/\b(\d+)\s*[- ]\s*(?=(?:ethanoyl|hydroxy|methoxy|oxo|methyl|ethyl|propyl|butyl|chloro|bromo|fluoro|iodo)\b)/g, '$1-');
    s = s.replace(/\bethanoyl\b/g, 'acetyl');
    s = s.replace(/\bcarboxylic\s+acid\b/g, 'carboxylic acid');
    s = s.replace(/\bmethoxy\s+carbonyl\b/g, 'methoxycarbonyl');
    s = s.replace(/\bethanoyl\s+oxy\b/g, 'acetyloxy');
    s = s.replace(/\bbenzon(ic|o)\s+anhydride\b/g, 'benzoic anhydride');
    return s;
  }

  /*
   * Conservative IUPAC-style fallback for several suffix/prefix patterns.
   * These are structural templates; locants are applied to a plain, acyclic
   * parent only when the shape is unambiguous. More involved parent selection
   * is deliberately left to the original parser.
   */
  var simpleParents = {
    meth: 1, eth: 2, prop: 3, but: 4, pent: 5, hex: 6, hept: 7,
    oct: 8, non: 9, dec: 10, undec: 11, dodec: 12
  };

  function parentSmiles(n, unsaturation) {
    var s = 'C'.repeat(n);
    if (unsaturation && unsaturation.position >= 1 && unsaturation.position < n) {
      var p = unsaturation.position - 1;
      s = s.slice(0, p) + 'C' + (unsaturation.kind === 'yne' ? '#' : '=') +
        'C' + s.slice(p + 2);
    }
    return s;
  }

  function parseSimpleSystematic(name, core) {
    // Examples supported here: butan-2-one, pentan-1-al (not preferred IUPAC),
    // hex-2-en-1-ol and straightforward suffix/prefix combinations.
    var s = name.replace(/\s+/g, '').replace(/-+/g, '-');
    var match = /^(\d+[-,])?((?:meth|eth|prop|but|pent|hex|hept|oct|non|dec|undec|dodec))(?:an|en|yn|ane|ene|yne)(?:-(\d+))?(?:-(ol|one|al|amine|thiol|nitrile|oicacid))?$/.exec(s);
    if (!match) return null;
    var n = simpleParents[match[2]];
    if (!n) return null;
    var sat = /yn/.test(s) ? 'yne' : /en/.test(s) ? 'ene' : null;
    var pos = +(match[3] || (sat ? 1 : 1));
    var smiles = parentSmiles(n, sat ? { kind: sat, position: pos } : null);
    var suffix = match[4];
    if (suffix === 'ol') {
      // Attach hydroxyl to a carbon with a free valence.
      var at = Math.max(1, Math.min(n, pos)) - 1;
      smiles = smiles.slice(0, at + 1) + '(O)' + smiles.slice(at + 1);
    } else if (suffix === 'one') {
      var kp = Math.max(2, Math.min(n - 1, pos)) - 1;
      smiles = smiles.slice(0, kp) + 'C(=O)' + smiles.slice(kp + 1);
    } else if (suffix === 'al') {
      smiles = 'O=C' + smiles.slice(1);
    } else if (suffix === 'oicacid') {
      smiles = 'O=C(O)' + smiles;
    } else if (suffix === 'amine') {
      smiles = 'N' + smiles;
    } else if (suffix === 'thiol') {
      smiles = smiles + 'S';
    } else if (suffix === 'nitrile') {
      smiles = smiles + 'C#N';
    }
    return smiles;
  }

  function resultFromSmiles(smiles, core, label) {
    if (!core || typeof core.parseSmiles !== 'function') {
      return { err: 'Chem3DCore.parseSmiles is unavailable; load this module after the Chem3D core.' };
    }
    var result = core.parseSmiles(smiles);
    if (!result || !result.atoms) return { err: result && result.err ? result.err : 'Could not parse fallback SMILES for ' + label };
    result.name = label;
    result.note = 'Name compatibility layer: ' + label;
    return result;
  }

  function resolve(raw, core) {
    core = core || root.Chem3DCore;
    var original = String(raw == null ? '' : raw);
    var cleaned = canonicalSpelling(original);
    if (!cleaned) return { err: 'Type a compound name, formula, or SMILES.' };

    var known = synonymMap[cleaned];
    if (!known && cleaned.replace(/\s/g, '') !== cleaned) {
      known = synonymMap[cleaned.replace(/\s/g, '')];
    }
    if (known) return resultFromSmiles(known, core, cleaned);

    // Common locant/suffix formatting mistakes: 3-Ethanoyl oxy -5- Methoxy ...
    if (/cyclohexene/.test(cleaned) && /carboxylic acid/.test(cleaned) && /acetyloxy|ethanoyloxy/.test(cleaned)) {
      return { err: 'This substituted cyclohexene carboxylic acid needs full parent-ring/locant parsing; please use a verified SMILES string for exact connectivity.' };
    }

    if (core && typeof core.fromName === 'function') {
      var byName = core.fromName(cleaned);
      if (byName && byName.src && /^smiles:/i.test(byName.src)) {
        var viaSmiles = core.parseSmiles(byName.src.replace(/^smiles:/i, ''));
        if (viaSmiles && viaSmiles.atoms) {
          viaSmiles.name = byName.name || cleaned;
          viaSmiles.note = byName.note || 'Interpreted by Chem3D core';
          return viaSmiles;
        }
      }
      if (byName && byName.src && !/^smiles:/i.test(byName.src)) {
        var viaFormula = core.resolve(byName.src);
        if (viaFormula && viaFormula.atoms) {
          viaFormula.name = byName.name || cleaned;
          viaFormula.note = byName.note || 'Interpreted by Chem3D core';
          return viaFormula;
        }
      }
    }

    if (core && typeof core.parseSmiles === 'function') {
      var systematic = parseSimpleSystematic(cleaned, core);
      if (systematic) {
        var g = resultFromSmiles(systematic, core, cleaned);
        if (g.atoms) return g;
      }
      if (typeof core.resolve === 'function') {
        var oldResult = core.resolve(original);
        if (oldResult && oldResult.atoms) return oldResult;
        return oldResult || { err: 'Could not resolve compound name.' };
      }
    }
    return { err: 'Chem3DCore is unavailable.' };
  }

  root.Chem3DNameCompat = {
    resolve: resolve,
    cleanName: cleanName,
    canonicalSpelling: canonicalSpelling,
    knownNames: function () { return Object.keys(synonymMap).sort(); },
    version: '1.0.0'
  };
})(typeof window !== 'undefined' ? window : globalThis);
