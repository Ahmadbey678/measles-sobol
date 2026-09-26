// model.js
// JavaScript port of simulate_measles() in R/model.R, used by the web page to
// run the model live in the browser. Kept line-for-line equivalent to the R
// version; tests/check_js_port.js checks the two agree.

(function (root) {
  const MCV2_EFFICACY = 0.95;
  const SIA_EFFICACY = 0.92;
  const STEPS_PER_YEAR = 26;

  // 77% at 9 months to 92% at 12 months (Uzicanin & Zimmerman 2011).
  function mcv1Efficacy(ageMonths) {
    return 0.77 + (ageMonths - 9) / 3 * (0.92 - 0.77);
  }

  // p: {R0, birth_rate, mcv1_cov, mcv2_cov, mcv1_age, sia_interval, sia_reach, importation}
  function simulate(p, opts) {
    const o = Object.assign({ years: 30, burnIn: 10, N: 1e6 }, opts || {});
    const N = o.N;
    const nSteps = o.years * STEPS_PER_YEAR;
    const mu = p.birth_rate / 1000 / STEPS_PER_YEAR;
    const e1 = mcv1Efficacy(p.mcv1_age);
    const v = p.mcv1_cov * e1 +
              p.mcv1_cov * p.mcv2_cov * (1 - e1) * MCV2_EFFICACY;
    const imp = p.importation * N / 1e6;
    const siaEvery = Math.round(p.sia_interval * STEPS_PER_YEAR);
    const siaHit = p.sia_reach * SIA_EFFICACY;

    let S = N / p.R0;
    let I = N * 1e-4;
    let R = N - S - I;

    const recordFrom = o.burnIn * STEPS_PER_YEAR;
    let cases = 0;
    let peakReff = 0;
    let stepsAbove = 0;
    const traj = { year: [], reff: [], casesPer100k: [] };

    for (let t = 1; t <= nSteps; t++) {
      const births = mu * N;
      const newInf = S * (1 - Math.exp(-p.R0 * (I + imp) / N));

      const Snext = S - newInf + births * (1 - v) - mu * S;
      const Rnext = R + I + births * v - mu * R;
      I = newInf;
      S = Snext;
      R = Rnext;

      const moved = (t % siaEvery === 0) ? S * siaHit : 0;
      S -= moved;
      R += moved;

      if (t > recordFrom) {
        cases += newInf;
        const reffBefore = p.R0 * (S + moved) / N;
        peakReff = Math.max(peakReff, reffBefore);
        if (p.R0 * S / N > 1) stepsAbove++;
        traj.year.push((t - recordFrom) / STEPS_PER_YEAR);
        traj.reff.push(p.R0 * S / N);
        traj.casesPer100k.push(newInf / N * 1e5);
      }
    }

    const recordedYears = o.years - o.burnIn;
    return {
      incidence: cases / recordedYears / N * 1e5,
      peakReff: peakReff,
      fracAbove: stepsAbove / (recordedYears * STEPS_PER_YEAR),
      routineProtected: v,
      traj: traj
    };
  }

  const api = { simulate, mcv1Efficacy };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.MeaslesModel = api;
})(this);
