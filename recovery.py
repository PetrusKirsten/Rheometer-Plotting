import re
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf

from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from statsmodels.stats.multicomp import pairwise_tukeyhsd

from rheology.plotting import Recovery
from rheology.plotting import powerLaw
from rheology.plotting import arraySplit
from rheology.plotting import exportFit

def extract_cacl2(sample_name):
    match = re.search(r'CL[-_](\d+)', sample_name)
    return int(match.group(1)) if match else 0

def extract_polymer(sample_name):
    if 'kCar' in sample_name:
        return 'kappa'
    elif 'iCar' in sample_name:
        return 'iota'
    else:
        return 'unknown'

def build_dataframe(data_list):
    records = []
    for entry in data_list:
        cacl2 = extract_cacl2(entry['Sample'])
        polymer = extract_polymer(entry['Sample'])
        formulation = f"{polymer} - CaCl2 {cacl2} mM"
        if np.isfinite(entry["k'"]):
            record = {
                'Polymer': polymer,
                'CaCl2_mM': cacl2,
                'Formulation': formulation,
                "G0_prime": entry["k'"]
            }
            records.append(record)
    df = pd.DataFrame(records)
    df = df.dropna()
    return df

def run_two_way_anova(df):
    if df.shape[0] < 2:
        raise ValueError("Not enough data to perform ANOVA. Please ensure the DataFrame has sufficient rows.")
    if df['Polymer'].nunique() > 1:
        formula = 'G0_prime ~ C(Polymer) + C(CaCl2_mM) + C(Polymer):C(CaCl2_mM)'
    else:
        formula = 'G0_prime ~ C(CaCl2_mM)'
    model = smf.ols(formula, data=df).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)
    return anova_table

def run_posthoc_tukey(df):
    df['Polymer_CaCl2'] = df['Polymer'] + ' - ' + df['CaCl2_mM'].astype(str) + ' mM'
    result = pairwise_tukeyhsd(endog=df['G0_prime'], groups=df['Polymer_CaCl2'], alpha=0.05)
    print(' \n', result.summary())
    return result

def plot_posthoc_tukey(result):
    plt.figure(figsize=(10, 6))
    result.plot_simultaneous(comparison_name=None)
    plt.title("Tukey HSD Post Hoc Test (Full Formulation Names)")
    plt.xlabel("Mean Difference in G'₀ with 95% CI")
    plt.ylabel("Formulation Comparisons")
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()

def plot_interaction(df):
    if df.empty:
        print("Empty DataFrame! Nothing to plot.")
        return
    dodge = df['Polymer'].nunique() > 1
    plt.figure(figsize=(10, 6))
    sns.pointplot(data=df, x='CaCl2_mM', y='G0_prime', hue='Polymer',
                  dodge=dodge, markers=['o', 's'], capsize=0.1, err_kws={'linewidth': 1}, palette='tab10')
    plt.title("Interaction between CaCl2 concentration and Polymer for G'₀")
    plt.ylabel("G'₀")
    plt.xlabel("CaCl₂ concentration (mM)")
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.show()


# def blends(folderPath):
#     filePath = [
#         # St
#         folderPath + "/10St/10_0WSt-viscRec_1.xlsx",
#         folderPath + "/10St/10_0WSt-viscRec_2.xlsx",
#
#         # St + kCar
#         folderPath + "/10St_kC/10_0WSt_kCar-viscoelasticRecovery-Flow_2a.xlsx",
#         folderPath + "/10St_kC/10_0WSt_kCar-viscoelasticRecovery-Flow_3a.xlsx",
#         # folderPath + "/10St_kC/10_0WSt_kCar-viscoelasticRecovery-Flow_4a.xlsx",
#
#         # St + iCar
#         folderPath + "/10St_iC/10_0WSt_iCar-viscoRecoveryandFlow_2.xlsx",
#         # folderPath + "10St_iC/10_0WSt_iCar-viscoRecoveryandFlow_1.xlsx",
#         folderPath + "/10St_iC/10_0WSt_iCar-viscoRecoveryandFlow_3.xlsx",
#         folderPath + "/10St_iC/10_0WSt_iCar-viscoRecoveryandFlow_4.xlsx",
#
#         # St/CL_7
#         folderPath + "/10St_CL_7/10_0St_CL-recovery-1.xlsx",
#         # folderPath + "/10St_CL_7/10_0St_CL-recovery-2.xlsx",
#         folderPath + "/10St_CL_7/10_0St_CL-recovery-3.xlsx",
#
#         # St + kCar/CL_7
#         folderPath + "/10St_kC_CL_7/10_0St_kC_CL-recovery-1.xlsx",
#         folderPath + "/10St_kC_CL_7/10_0St_kC_CL-recovery-3_off.xlsx",
#         folderPath + "/10St_kC_CL_7/10_0St_kC_CL-recovery-4.xlsx",
#
#         # St + iCar/CL_7
#         folderPath + "/10St_iC_CL_7/10_0St_iC_CL-recovery-1.xlsx",
#         folderPath + "/10St_iC_CL_7/10_0St_iC_CL-recovery-2.xlsx",
#         folderPath + "/10St_iC_CL_7/10_0St_iC_CL-recovery-3.xlsx",
#
#         # St/CL_14
#         folderPath + "/10St_CL_14/St_CL_14-viscoelasticRecovery-1.xlsx",
#         folderPath + "/10St_CL_14/St_CL_14-viscoelasticRecovery-2.xlsx",
#         folderPath + "/10St_CL_14/St_CL_14-viscoelasticRecovery-3.xlsx",
#         folderPath + "/10St_CL_14/St_CL_14-viscoelasticRecovery-4.xlsx",
#
#         # St + iCar/CL_14
#         folderPath + "/10St_iC_CL_14/0St_iC_CL_14-viscoelasticRecovery-1.xlsx",
#         folderPath + "/10St_iC_CL_14/0St_iC_CL_14-viscoelasticRecovery-2.xlsx",
#         folderPath + "/10St_iC_CL_14/0St_iC_CL_14-viscoelasticRecovery-3.xlsx",
#         folderPath + "/10St_iC_CL_14/0St_iC_CL_14-viscoelasticRecovery-4.xlsx",
#
#         # St/CL_28
#         folderPath + "/10St_CL_28/St_CL_28-viscoelasticRecovery-1.xlsx",
#         folderPath + "/10St_CL_28/St_CL_28-viscoelasticRecovery-2.xlsx",
#         folderPath + "/10St_CL_28/St_CL_28-viscoelasticRecovery-3.xlsx",
#         folderPath + "/10St_CL_28/St_CL_28-viscoelasticRecovery-4.xlsx",
#     ]
#
#     keySamples = {
#         # No CL
#         'St': [], 'St + kCar': [], 'St + iCar': [],
#         # CL 7
#         'St/CL_7': [], 'St + kCar/CL_7': [], 'St + iCar/CL_7': [],
#         # CL 14
#         'St/CL_14': [],  # 'St + kCar/CL_14': [], TBD
#         'St + iCar/CL_14': [],
#         # CL 28
#         'St/CL_28': [],
#         # 'St + kCar/CL_28': [], TBD
#         # 'St + iCar/CL_28': [] TBD
#     }
#     nSamples = [
#         # No CL
#         2, 2, 3,
#         # CL 7
#         2, 3, 3,
#         # CL 14
#         4, 4,
#         # CL 28
#         4,
#     ]
#     cSamples = [
#         # No CL
#         'silver', 'hotpink', 'lightskyblue',
#         # CL 7
#         'grey', 'mediumvioletred', 'royalblue',
#         # CL 14
#         'dimgrey', 'r',
#         # CL 28
#         'k',
#     ]
#
#     starches = Recovery(
#         filePath, 'blends',
#         keySamples, nSamples, cSamples)
#
#     starches.plotGmodulus(
#         f"Elastic modulus $G'$ (Pa)", (1 * 10 ** 0, 1 * 10 ** 5),
#         f'Frequency (Hz)', (.075, 100),
#         show=False)
#
#     starches.plotBars(
#         [0.5, 1800, 500, 30],
#         [None, 4, None, None],
#         show=False)
#
#     starches.plotHeatMap(show=False)

def kappa(folderPath):
    filePath = [
        # kC
        folderPath + "/kC/kC-viscoelasticRecovery-1.xlsx",
        folderPath + "/kC/kC-viscoelasticRecovery-2.xlsx",
        folderPath + "/kC/kC-viscoelasticRecovery-3.xlsx",

        # kC CL 7
        folderPath + "/kC_CL_7/kC_CL-viscoelasticRecovery-1.xlsx",
        folderPath + "/kC_CL_7/kC_CL-viscoelasticRecovery-2.xlsx",
        folderPath + "/kC_CL_7/kC_CL-viscoelasticRecovery-3.xlsx",
        folderPath + "/kC_CL_7/kC_CL-viscoelasticRecovery-4.xlsx",

        # kC CL 14
        folderPath + "/kC_CL_14/kC_CL_14-viscoelasticRecovery-1.xlsx",
        folderPath + "/kC_CL_14/kC_CL_14-viscoelasticRecovery-2.xlsx",

        # kC CL 21
        folderPath + "/kC_CL_21/kC_CL_21-viscoelasticRecovery-1.xlsx",
        # folderPath + "/kC_CL_21/kC_CL_21-viscoelasticRecovery-2.xlsx",
        folderPath + "/kC_CL_21/kC_CL_21-viscoelasticRecovery-3.xlsx",

        # kC CL 28
        folderPath + "/kC_CL_28/kC_CL_28-viscoelasticRecovery-1.xlsx",
        folderPath + "/kC_CL_28/kC_CL_28-viscoelasticRecovery-2.xlsx",

        # kC CL 42
        folderPath + "/kC_CL_42/kC_CL_42-viscoelasticRecovery-1.xlsx",
        folderPath + "/kC_CL_42/kC_CL_42-viscoelasticRecovery-2.xlsx",
        # folderPath + "/kC_CL_42/kC_CL_42-viscoelasticRecovery-3.xlsx",
        # folderPath + "/kC_CL_42/kC_CL_42-viscoelasticRecovery-4.xlsx",
    ]

    kappa_keySamples = {
        'kCar': [],
        'kCar/CL-7': [],
        'kCar/CL-14': [],
        'kCar/CL-21': [],
        'kCar/CL-28': [],
        'kCar/CL-42': []
    }
    kappa_nSamples = [3, 4, 2, 3, 2, 4]
    kappa_cSamples = ['lightsteelblue', '#A773FF', '#892F99', '#AB247B', '#E64B83', '#FF0831']

    kappas = Recovery(
        filePath, 'kappas',
        kappa_keySamples, kappa_nSamples, kappa_cSamples)

    kappas.plotGmodulus(
        f"Elastic modulus $G'$ (Pa)", (1 * 10 ** (-2), 1 * 10 ** 5),
        f'Frequency (Hz)', (.075, 100),
        show=False, save=False)

    kappas.plotBars(
        [0, None, None, None],
        significance_list=[
            ['A', 'B', 'B', 'B', 'B', 'B',
             'A*', 'B*', 'C', 'C', 'B*', 'B*',],

            ['X', 'A', 'B', 'C', 'D', 'E',
             'X', 'A*', 'B*', 'C*', 'D*', 'D*',],

            ['X', 'A', 'B', 'C', 'D', 'E',
             'X', 'A*', 'B*', 'C*', 'D*', 'D*', ],

            ['X', 'A', 'B', 'B', 'B', 'B',
             'X', 'A*', 'A*', 'B*', 'C*', 'C*',],
        ],
        show=False, save=True)


def iota(folderPath):
    filePath = [
        # iC CL 7
        folderPath + "/iC_CL_7/iC_CL_7-viscoelasticRecovery-1.xlsx",
        folderPath + "/iC_CL_7/iC_CL_7-viscoelasticRecovery-2.xlsx",

        # iC CL 14
        folderPath + "/iC_CL_14/iC_CL_14-viscoelasticRecovery-1.xlsx",
        folderPath + "/iC_CL_14/iC_CL_14-viscoelasticRecovery-2.xlsx",

        # iC CL 21
        folderPath + "/iC_CL_21/iC_CL_21-viscoelasticRecovery-1.xlsx",
        folderPath + "/iC_CL_21/iC_CL_21-viscoelasticRecovery-2.xlsx",

        # iC CL 28
        folderPath + "/iC_CL_28/iC_CL_28-viscoelasticRecovery-1.xlsx",
        folderPath + "/iC_CL_28/iC_CL_28-viscoelasticRecovery-2.xlsx",

        # iC CL 42
        folderPath + "/iC_CL_42/iC_CL_42-viscoelasticRecovery-1.xlsx",
        folderPath + "/iC_CL_42/iC_CL_42-viscoelasticRecovery-2.xlsx",
        folderPath + "/iC_CL_42/iC_CL_42-viscoelasticRecovery-3.xlsx",
    ]

    keySamples = {
        'iCar/CL_7': [],
        'iCar/CL_14': [],
        'iCar/CL_21': [],
        'iCar/CL_28': [],
        'iCar/CL_42': []
    }
    nSamples = [2, 2, 2, 2, 3, ]
    cSamples = ['greenyellow', '#80ed99', '#57cc99', '#38a3a5', '#22577a']

    iotas = Recovery(
        filePath, 'iota',
        keySamples, nSamples, cSamples)

    iotas.plotGmodulus(
        f"Elastic modulus $G'$ (Pa)", (1 * 10 ** (-2), 3 * 10 ** 3),
        f'Frequency (Hz)', (.075, 100),
        show=False)

    iotas.plotBars(
        [None, None, None, None],
        significance_list=[
            ['A', 'B', 'B', 'B', 'B',
             'A', 'A', 'B*', 'B', 'A',],

            [' ', 'A', 'B', 'A', 'A',
             ' ', 'A*', 'B*', 'B*', 'A*',],

            ['A', 'B', 'C', 'B', 'B',
             'A*', 'B*', 'C*', 'C*', 'C*',],

            [' ', 'A', 'A', 'A', 'A',
             'A*', 'B*', 'C*', 'C*', 'B*',],
        ],
        show=False, save=True)

    # iotas.plotHeatMap(show=False)


def starch(folderPath):
    filePath = [
        # St
        folderPath + "/10St/10_0WSt-viscRec_1.xlsx",
        folderPath + "/10St/10_0WSt-viscRec_2.xlsx",

        # St/CL_7
        folderPath + "/10St_CL_7/10_0St_CL-recovery-1.xlsx",
        # folderPath + "/10St_CL_7/10_0St_CL-recovery-2.xlsx",
        folderPath + "/10St_CL_7/10_0St_CL-recovery-3.xlsx",

        # St/CL_14
        folderPath + "/10St_CL_14/St_CL_14-viscoelasticRecovery-1.xlsx",
        folderPath + "/10St_CL_14/St_CL_14-viscoelasticRecovery-2.xlsx",
        folderPath + "/10St_CL_14/St_CL_14-viscoelasticRecovery-3.xlsx",
        folderPath + "/10St_CL_14/St_CL_14-viscoelasticRecovery-4.xlsx",

        # St/CL_28
        folderPath + "/10St_CL_28/St_CL_28-viscoelasticRecovery-1.xlsx",
        folderPath + "/10St_CL_28/St_CL_28-viscoelasticRecovery-2.xlsx",
        folderPath + "/10St_CL_28/St_CL_28-viscoelasticRecovery-3.xlsx",
        folderPath + "/10St_CL_28/St_CL_28-viscoelasticRecovery-4.xlsx",
    ]

    keySamples = {
        # No CL
        'St': [],
        # CL 7
        'St/CL_7': [],
        # CL 14
        'St/CL_14': [],
        # CL 28
        'St/CL_21': [],
    }
    nSamples = [
        # No CL
        2,
        # CL 7
        2,
        # CL 14
        4,
        # CL 28
        4,
    ]
    cSamples = [
        # No CL
        '#E1C96B',
        # CL 7
        '#FFE138',
        # CL 14
        '#F1A836',
        # CL 28
        '#E36E34',
    ]

    starches = Recovery(
        filePath, 'starches',
        keySamples, nSamples, cSamples)

    starches.plotGmodulus(
        f"Elastic modulus $G'$ (Pa)", (1 * 10 ** 0, 1 * 10 ** 5),
        f'Frequency (Hz)', (.075, 100),
        show=False, save=False)

    starches.plotBars(
        [None, 4, None, None],
        significance_list=[
            ['A', 'B', 'A', 'A',
             'A', 'B*', 'A', 'A',],

            ['A', 'B', 'C', 'C',
             'A*', 'B*', 'C*', 'D*',],

            ['A', 'B', 'C', 'C',
             'A*', 'B*', 'C*', 'D*',],

            ['A', 'B', 'A', 'A',
             'A', ' ', 'A', 'A',],
        ],
        show=False, save=True)

    # starches.plotHeatMap(show=False, save=True)


def starch_kappa(folderPath):
    filePath = [
        # St + kCar
        folderPath + "/10St_kC/10_0WSt_kCar-viscoelasticRecovery-Flow_2a.xlsx",
        folderPath + "/10St_kC/10_0WSt_kCar-viscoelasticRecovery-Flow_3a.xlsx",
        # folderPath + "/10St_kC/10_0WSt_kCar-viscoelasticRecovery-Flow_4a.xlsx",

        # St + kCar/CL_7
        folderPath + "/10St_kC_CL_7/10_0St_kC_CL-recovery-1.xlsx",
        folderPath + "/10St_kC_CL_7/10_0St_kC_CL-recovery-3_off.xlsx",
        folderPath + "/10St_kC_CL_7/10_0St_kC_CL-recovery-4.xlsx",

        # St + kCar/CL_14
        folderPath + "/10St_kC_CL_14/St_kC_CL_14-viscoelasticRecovery-1.xlsx",
        folderPath + "/10St_kC_CL_14/St_kC_CL_14-viscoelasticRecovery-2.xlsx",
        folderPath + "/10St_kC_CL_14/St_kC_CL_14-viscoelasticRecovery-3.xlsx",
        folderPath + "/10St_kC_CL_14/St_kC_CL_14-viscoelasticRecovery-4.xlsx",

        # St + kCar/CL_28
        folderPath + "/10St_kC_CL_28/St_kC_CL_28-viscoelasticRecovery-1.xlsx",
        folderPath + "/10St_kC_CL_28/St_kC_CL_28-viscoelasticRecovery-2.xlsx",
        folderPath + "/10St_kC_CL_28/St_kC_CL_28-viscoelasticRecovery-3.xlsx",

    ]

    keySamples = {
        # No CL
        'St + kCar': [],
        # CL 7
        'St + kCar/CL_7': [],
        # CL 14
        'St + kCar/CL_14': [],
        # CL 28
        'St + kCar/CL_21': [],
    }
    nSamples = [
        # No CL
        2,
        # CL 7
        3,
        # CL 14
        4,
        # CL 28
        3,
    ]
    cSamples = [
        # No CL
        '#F780A7',
        # CL 7
        '#CC69B5',
        # CL 14
        '#A251C3',
        # CL 28
        '#773AD1',
    ]

    starches = Recovery(
        filePath, 'starches-kappas',
        keySamples, nSamples, cSamples)

    starches.plotGmodulus(
        f"Elastic modulus $G'$ (Pa)", (1 * 10 ** 0, 1 * 10 ** 5),
        f'Frequency (Hz)', (.075, 100),
        show=False, save=False)

    starches.plotBars(
        [None, None, None, None],
        significance_list=[
            ['A', 'B', 'C', 'C',
             'A*', 'B*', 'C*', 'C*',],

            ['A', 'B', 'C', 'D',
             'A*', 'B*', 'C*', 'D*',],

            ['A', 'B', 'C', 'B',
             'A*', 'B*', 'C*', 'D*', ],

            ['A', 'B', 'C', 'C',
             'A*', 'B*', 'B*', 'A*',],
        ],
        show=False, save=True)

    # starches.plotHeatMap(show=False, save=True)


def starch_iota(folderPath):
    filePath = [
        # St + iCar
        folderPath + "/10St_iC/10_0WSt_iCar-viscoRecoveryandFlow_2.xlsx",
        # folderPath + "10St_iC/10_0WSt_iCar-viscoRecoveryandFlow_1.xlsx",
        folderPath + "/10St_iC/10_0WSt_iCar-viscoRecoveryandFlow_3.xlsx",
        folderPath + "/10St_iC/10_0WSt_iCar-viscoRecoveryandFlow_4.xlsx",

        # St + iCar/CL_7
        folderPath + "/10St_iC_CL_7/10_0St_iC_CL-recovery-1.xlsx",
        folderPath + "/10St_iC_CL_7/10_0St_iC_CL-recovery-2.xlsx",
        folderPath + "/10St_iC_CL_7/10_0St_iC_CL-recovery-3.xlsx",

        # St + iCar/CL_14
        folderPath + "/10St_iC_CL_14/0St_iC_CL_14-viscoelasticRecovery-1.xlsx",
        folderPath + "/10St_iC_CL_14/0St_iC_CL_14-viscoelasticRecovery-2.xlsx",
        folderPath + "/10St_iC_CL_14/0St_iC_CL_14-viscoelasticRecovery-3.xlsx",
        folderPath + "/10St_iC_CL_14/0St_iC_CL_14-viscoelasticRecovery-4.xlsx",

        # St + iCar/CL_28
        folderPath + "/10St_iC_CL_28/St_iC_CL_28-viscoelasticRecovery-1.xlsx",
        folderPath + "/10St_iC_CL_28/St_iC_CL_28-viscoelasticRecovery-2.xlsx",
        folderPath + "/10St_iC_CL_28/St_iC_CL_28-viscoelasticRecovery-3.xlsx",

    ]

    keySamples = {
        # No CL
        'St + iCar': [],
        # CL 7
        'St + iCar/CL_7': [],
        # CL 14
        'St + iCar/CL_14': [],
        # CL 28
        'St + iCar/CL_21': []
    }
    nSamples = [
        # No CL
        3,
        # CL 7
        3,
        # CL 14
        4,
        # CL 28
        3,
    ]
    cSamples = [
        # No CL
        'lightskyblue',
        # CL 7
        '#62BDC1',
        # CL 14
        '#31A887',
        # CL 28
        '#08653A',
    ]

    starches = Recovery(
        filePath, 'starches-iotas',
        keySamples, nSamples, cSamples)

    starches.plotGmodulus(
        f"Elastic modulus $G'$ (Pa)", (1 * 10 ** 0, 1 * 10 ** 5),
        f'Frequency (Hz)', (.075, 100),
        show=False, save=False)

    starches.plotBars(
        [0.18, 3000, 250, 18],
        significance_list=[
            ['A', 'B', 'C', 'D',
             'A*', 'B*', 'B*', 'C*',],

            ['A', 'B', 'C', 'D',
             'A*', 'B*', 'A*', 'C*',],

            ['A', 'B', 'C', 'B',
             'A*', 'A*', 'B*', 'C*', ],

            ['A', 'A', 'B', 'C',
             'A*', 'B*', 'B*', 'C*',],
        ],
        show=False, save=True)

    # starches.plotHeatMap(show=False, save=True)


if __name__ == '__main__':
    # path = "C:/Users/petrus.kirsten/PycharmProjects/Rheometer-Plotting/data/by sample"  # CEBB
    # path = "C:/Users/Petrus Kirsten/Documents/GitHub/RheometerPlots/data/by sample"   # Personal
    path = "D:/Documents/GitHub/Rheometer-Plotting/data/by sample"   # New Personal

    # kappa(path)
    # iota(path)
    # blends(path)
    # starch(path)
    # starch_kappa(path)
    starch_iota(path)

    # df = build_dataframe(statisticalData)
    # anova_table = run_two_way_anova(df)
    # print('\n', '\n', anova_table, '\n')
    #
    # if anova_table['PR(>F)'][0] < 0.05:
    #     print("Significative differences found in ANOVA.")
    # else:
    #     print("No significant differences found in ANOVA.")
    #
    # tukey_result = run_posthoc_tukey(df)
    # plot_posthoc_tukey(tukey_result)
    # # plot_interaction(df)

    plt.show()
