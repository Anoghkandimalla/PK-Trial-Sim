import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Typical rate constants (per hour): plasma->tumor, tumor->plasma,
# plasma->heart, heart->plasma, elimination
TYPICAL_K = np.array([0.3, 0.1, 0.05, 0.2, 0.15])
VARIABILITY = 0.3          # log-normal spread between patients
DOSE = 10.0                # starting plasma concentration
TUMOR_TARGET = 2.0         # needed in tumor to be effective
HEART_LIMIT = 0.5          # above this in heart = toxic
HOURS = 48
N_PATIENTS = 1000


def model(t, y, kpt, ktp, kph, khp, kel):
    Cp, Ct, Ch = y
    dCp = -(kpt + kph + kel) * Cp + ktp * Ct + khp * Ch
    dCt = kpt * Cp - ktp * Ct
    dCh = kph * Cp - khp * Ch
    return [dCp, dCt, dCh]


def simulate(k):
    t_eval = np.linspace(0, HOURS, 500)
    return solve_ivp(model, [0, HOURS], [DOSE, 0, 0], args=tuple(k), t_eval=t_eval)


def run_trial(seed=42):
    rng = np.random.default_rng(seed)
    successes = 0
    for _ in range(N_PATIENTS):
        k = rng.lognormal(np.log(TYPICAL_K), VARIABILITY)
        sol = simulate(k)
        if sol.y[1].max() >= TUMOR_TARGET and sol.y[2].max() < HEART_LIMIT:
            successes += 1
    return successes / N_PATIENTS


def plot_typical():
    sol = simulate(TYPICAL_K)
    plt.plot(sol.t, sol.y[0], label="Plasma")
    plt.plot(sol.t, sol.y[1], label="Tumor")
    plt.plot(sol.t, sol.y[2], label="Heart")
    plt.axhline(TUMOR_TARGET, ls="--", c="gray", label="Tumor target")
    plt.axhline(HEART_LIMIT, ls=":", c="red", label="Heart limit")
    plt.xlabel("Hours")
    plt.ylabel("Concentration")
    plt.title("Drug concentration, typical patient")
    plt.legend()
    plt.savefig("concentration_curves.png", dpi=150)


if __name__ == "__main__":
    p = run_trial()
    print(f"Predicted probability of success: {p:.1%}")
    plot_typical()
    print("Saved concentration_curves.png")
