"""
Independent cross-check of the R analysis, written in Python with SALib.

Same model, same input ranges, but a different sampling scheme
(Saltelli / Sobol' sequence instead of plain random draws) and a different
estimator implementation. If the ranking of inputs matches the R results,
the conclusions do not depend on one package's implementation.

Run from the project root:   python python/crosscheck_salib.py
"""
import numpy as np
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol as sobol_analyze

# Same ranges as R/model.R
PARAMS = [
    ("R0",           12,   18),
    ("birth_rate",   15,   40),
    ("mcv1_cov",     0.60, 0.95),
    ("mcv2_cov",     0.30, 0.90),
    ("mcv1_age",     9,    12),
    ("sia_interval", 2,    5),
    ("sia_reach",    0.20, 0.80),
    ("importation",  0.1,  5.0),
]
problem = {
    "num_vars": len(PARAMS),
    "names":  [p[0] for p in PARAMS],
    "bounds": [[p[1], p[2]] for p in PARAMS],
}

MCV2_EFFICACY = 0.95
SIA_EFFICACY = 0.92


def simulate(X, years=30, burn_in=10, N=1e6):
    """Vectorised port of simulate_measles() in R/model.R."""
    R0, br, c1, c2, age, sia_int, reach, imp_rate = X.T
    spy = 26
    mu = br / 1000 / spy
    e1 = 0.77 + (age - 9) / 3 * (0.92 - 0.77)
    v = c1 * e1 + c1 * c2 * (1 - e1) * MCV2_EFFICACY
    imp = imp_rate * N / 1e6
    sia_every = np.round(sia_int * spy).astype(int)
    sia_hit = reach * SIA_EFFICACY

    S = N / R0
    I = np.full(len(X), N * 1e-4)
    R = N - S - I
    cases = np.zeros(len(X))
    peak = np.zeros(len(X))

    for t in range(1, years * spy + 1):
        births = mu * N
        new_inf = S * (1 - np.exp(-R0 * (I + imp) / N))
        S_next = S - new_inf + births * (1 - v) - mu * S
        R_next = R + I + births * v - mu * R
        I, S, R = new_inf, S_next, R_next

        moved = np.where(t % sia_every == 0, S * sia_hit, 0.0)
        S = S - moved
        R = R + moved

        if t > burn_in * spy:
            cases += new_inf
            peak = np.maximum(peak, R0 * (S + moved) / N)

    incidence = cases / (years - burn_in) / N * 1e5
    return np.log10(incidence), peak


if __name__ == "__main__":
    X = sobol_sample.sample(problem, 8192, calc_second_order=False, seed=2027)
    print(f"Model runs: {len(X):,}")
    y_inc, y_peak = simulate(X)

    for title, y in [("Average burden (log10 cases/100k/yr)", y_inc),
                     ("Worst-case outbreak risk (peak Reff)", y_peak)]:
        Si = sobol_analyze.analyze(problem, y, calc_second_order=False,
                                   seed=2027, print_to_console=False)
        order = np.argsort(-Si["ST"])
        print(f"\n== {title} ==")
        print(f"{'input':<14}{'first_order':>12}{'total':>8}")
        for i in order:
            print(f"{problem['names'][i]:<14}{Si['S1'][i]:>12.3f}{Si['ST'][i]:>8.3f}")
