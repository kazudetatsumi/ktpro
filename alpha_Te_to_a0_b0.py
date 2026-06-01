import numpy as np
import matplotlib.pyplot as plt
"""
This script obtains a set of sample parameters Te and alpha of single Bragg
edge transmission parameters of my run_rits.py so that the corresponding rits
parameters a0 and b0 have similar histograms to the experimentally analyzed
results shown by Dr. Oikawa and transmission values are physically valid.
The first code was generated AI studio Gemini 3.1 Pro preview, then I modified
so that script composed of functions and added some English description.
Kazuyoshi TATSUMI 2026/05/24
To Fight back Against the vultures by publishing a paper of 3D denoise on
Bragg edge imaging.
"""


def get_constants():
    # starting and ending tof points [mirosec]
    ti, te = 23000, 26020
    # exchanged constants
    C, k = te - ti, te / (te - ti)

    # desired ranges of a0 and b0 of rits background parameters:
    # Transmission=T(tof)=exp(-a0 - b0*0.0001*tof[microsec])
    # Te = T(te), Tiext=T(ti), alpha=Tiext/Te
    a0_range = (0.0, 0.05)
    b0_range = (0.0, 0.1)

    # set parameters for weight function of the second dice
    # (set them to show the sample distribution close to
    # experimentally analyzed distribution of a0 and b0)
    mu_a0, sigma_a0 = 0.02, 0.008              # a0: a narrow hill around 0.0
    #n_b0 = 3.0                                # b0: 3乗で右肩上がりにする
    mu_b0, sigma_b0 = mu_a0*2.5, sigma_a0*2.5  # b0: a 2.5 times expanded
                                               # hill than a0's.
    return ti, te, C, k, a0_range, b0_range, mu_a0, sigma_a0, mu_b0, sigma_b0


def get_weights(a0, b0, mu_a0, sigma_a0, mu_b0, sigma_b0):
    # weight of a0 (Gaussian)
    wa0 = np.exp(-(a0 - mu_a0)**2 / (2 * sigma_a0**2))
    # weight of b0  (positively skewed function:
    # normalizing [0. 0.1] to [0. 1.] and taking nth power)
    # wb0 = (b0 / b0_range[1])**n_b0
    # wight of b0 (Gaussian)
    wb0 = np.exp(-(b0 - mu_b0)**2 / (2 * (sigma_b0)**2))
    return wa0 * wb0


# execute Monte Carlo sampling
def execute_montecalro_sampling(
    a0_range, b0_range, k, C, mu_a0, sigma_a0, mu_b0, sigma_b0,
    n_target=5000,    # gather 5000 points
    rg=None,
):
    samples_alphaTe = []  # (alpha, Te)
    samples_a0b0 = []     # (a0, b0)
    # outer flame for sampling from the trapezoid region
    alpha_min, alpha_max = 1.0, np.exp(b0_range[1] * C / 10000)
    Te_min, Te_max = np.exp(-k * np.log(alpha_max) - a0_range[1]), 1.0
    while len(samples_alphaTe) < n_target:
        # 1. 1st dice: uniform sampling within a rectangle enveloping
        #    the trapezoid region
        if rg is None:
            try_alpha = np.random.uniform(alpha_min, alpha_max)
            try_Te = np.random.uniform(Te_min, Te_max)
        else:
            try_alpha = rg.uniform(alpha_min, alpha_max)
            try_Te = rg.uniform(Te_min, Te_max)
        # checking whether the sample is inside of the trapezoid
        a0 = -k * np.log(try_alpha) - np.log(try_Te)
        b0 = 10000 * np.log(try_alpha) / C
        if (a0_range[0] <= a0 <= a0_range[1]) and\
                (b0_range[0] <= b0 <= b0_range[1]):
            # 2. 2nd dice: weighting the sampling of a0 and b0.
            p = get_weights(a0, b0, mu_a0, sigma_a0, mu_b0, sigma_b0)
            if rg is None:
                random_value = np.random.rand()
            else:
                random_value = rg.random()
                if random_value < p:
                    samples_alphaTe.append([try_alpha, try_Te])
                    samples_a0b0.append([a0, b0])
    return np.array(samples_alphaTe), np.array(samples_a0b0)


# visualization of the results
def plot_results(samples_alphaTe, samples_a0b0):
    fig, ax = plt.subplots(1, 2, figsize=(14, 5))
    # Left figure: distribution on alpha-Te region
    ax[0].scatter(samples_alphaTe[:, 0], samples_alphaTe[:, 1], s=1, alpha=0.5,
                  color='green')
    ax[0].set_title("Sampled Points in ($\\alpha, T_e$) Space")
    ax[0].set_xlabel("$\\alpha$")
    ax[0].set_ylabel("$T_e$")
    # Right figure: histogram of a0 and b0 like to experimental distributions.
    ax2 = ax[1].twinx()
    ax[1].hist(samples_a0b0[:, 0], bins=40, color='blue', alpha=0.5,
               label='a0 (Left axis)')
    ax2.hist(samples_a0b0[:, 1], bins=40, color='red', alpha=0.5,
             label='b0 (Right axis)')
    ax[1].set_title("Resulting Histograms of $a_0$ and $b_0$")
    ax[1].set_xlabel("Value")
    ax[1].legend(loc='upper left')
    ax2.legend(loc='upper right')
    plt.tight_layout()
    plt.show()


def run_alpha_Te_to_a0_b0(n_target=5000, plot=False, rg=None):
    ti, te, C, k, a0_range, b0_range, mu_a0, sigma_a0, mu_b0, sigma_b0 =\
        get_constants()
    samples_alphaTe, samples_a0b0 = execute_montecalro_sampling(
            a0_range, b0_range, k, C, mu_a0, sigma_a0, mu_b0, sigma_b0,
            n_target=n_target, rg=rg)
    if plot:
        plot_results(samples_alphaTe, samples_a0b0)
    return samples_alphaTe, samples_a0b0
