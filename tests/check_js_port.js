// Checks that the browser model (docs/model.js) reproduces R/model.R.
// Reference values were produced by simulate_measles() in R.
// Run from the project root:   node tests/check_js_port.js

const { simulate } = require("../docs/model.js");

const cases = [
  { p: { R0: 15, birth_rate: 35, importation: 1, mcv1_cov: 0.65, mcv2_cov: 0.40,
         mcv1_age: 9, sia_interval: 4, sia_reach: 0.40 },
    incidence: 750.1559534295, peakReff: 1.1744139223 },
  { p: { R0: 15, birth_rate: 35, importation: 1, mcv1_cov: 0.92, mcv2_cov: 0.85,
         mcv1_age: 9, sia_interval: 3, sia_reach: 0.60 },
    incidence: 0.7955953764, peakReff: 0.3461103496 },
  { p: { R0: 13.3, birth_rate: 22, importation: 3.7, mcv1_cov: 0.81, mcv2_cov: 0.55,
         mcv1_age: 10.5, sia_interval: 2.7, sia_reach: 0.33 },
    incidence: 9.7531150109, peakReff: 0.6573603626 },
];

let ok = true;
for (const c of cases) {
  const r = simulate(c.p);
  const d1 = Math.abs(r.incidence - c.incidence) / c.incidence;
  const d2 = Math.abs(r.peakReff - c.peakReff) / c.peakReff;
  const pass = d1 < 1e-8 && d2 < 1e-8;
  ok = ok && pass;
  console.log(`${pass ? "PASS" : "FAIL"}  incidence ${r.incidence.toFixed(6)} vs ${c.incidence}` +
              `   peakReff ${r.peakReff.toFixed(6)} vs ${c.peakReff}`);
}
process.exit(ok ? 0 : 1);
