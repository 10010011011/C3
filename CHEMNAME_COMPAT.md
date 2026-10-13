# Chem3D name compatibility layer

New file: [chemname-compat.js](./chemname-compat.js)

This module is deliberately isolated. It does not patch or replace the existing parser and does not change the UI until explicitly integrated.

## What it covers

- Retained phenol names and isomers: naphthol, catechol (including the typo `catchol`), resorcinol, hydroquinone/quinol, pyrogallol and cresols.
- Fused/non-benzenoid and heteroaromatic compounds: naphthalene, anthracene, phenanthrene, azulene, pyridine, pyrrole, furan, thiophene, indole and quinoline.
- Carbonyl and amide examples: acetanilide, acetone, acetophenone, benzaldehyde, benzophenone and benzil.
- Carboxylic acids, positional isomers, esters and anhydrides: cinnamic acid, phthalic acid, isophthalic acid, terephthalic acid, methyl benzoate and benzoic anhydride.
- Name cleanup for case, repeated whitespace, hyphen spacing, and common `benzonic acid` / `catchol` spelling variants.

## Integration

After the inline `Chem3DCore` definition in `index.html`, load this file before the application script:

```html
<script src="./chemname-compat.js"></script>
```

In the application's single `build()` resolver call, use:

```js
const G = window.Chem3DNameCompat
  ? window.Chem3DNameCompat.resolve(smi, Core)
  : smart(smi);
```

Keep this opt-in until testing the compounds below. Existing parsing remains the fallback for formulas, SMILES, and other names.

## Important limitations

This is a curated compatibility table, not a universal chemical nomenclature engine. Isomers such as 1- and 2-naphthol require distinct structures; confirm the connectivity against a chemistry reference before scientific use. The complex substituted cyclohexene carboxylic-acid example needs full locant and parent selection logic and is not safely reconstructed from an unverified name string. That example should be given as a verified SMILES or added only after its exact intended IUPAC name and structure are confirmed.

The compatibility resolver's API returns the same graph shape expected by the existing renderer: `{ atoms, bonds, name, note }`, or `{ err }`.
