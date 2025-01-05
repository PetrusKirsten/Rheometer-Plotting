from matplotlib import pyplot as plt
from rheology.plotting import Recovery


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
        [.6, 10000, 850, 16],
        [0, None, None, None],
        show=False, save=True)

    # kappas.plotHeatMap(show=False, save=False)


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
        [3.75, 50, 9, 8],
        [None, None, None, None],
        show=False, save=True)

    iotas.plotHeatMap(show=False)


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
        '#FFE138',
        # CL 7
        '#F1A836',
        # CL 14
        '#EA8B35',
        # CL 28
        '#E36E34',
    ]

    starches = Recovery(
        filePath, 'starches',
        keySamples, nSamples, cSamples)

    starches.plotGmodulus(
        f"Elastic modulus $G'$ (Pa)", (1 * 10 ** 0, 1 * 10 ** 5),
        f'Frequency (Hz)', (.075, 100),
        show=False, save=True)

    starches.plotBars(
        [0.6, 1800, 150, 25],
        [None, 4, None, None],
        show=False, save=True)

    starches.plotHeatMap(show=False, save=True)


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
        'hotpink',
        # CL 7
        'mediumvioletred',
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
        show=False, save=True)

    starches.plotBars(
        [0.18, 12500, 850, 30],
        [None, None, None, None],
        show=False, save=True)

    starches.plotHeatMap(show=False, save=True)


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
        '#34D4CF',
        # CL 7
        '#2B7ED4',
        # CL 14
        '#1941B8',
        # CL 28
        '#00934E',
    ]

    starches = Recovery(
        filePath, 'starches-iotas',
        keySamples, nSamples, cSamples)

    starches.plotGmodulus(
        f"Elastic modulus $G'$ (Pa)", (1 * 10 ** 0, 1 * 10 ** 5),
        f'Frequency (Hz)', (.075, 100),
        show=False, save=True)

    starches.plotBars(
        [0.18, 3000, 250, 18],
        [None, 4, None, None],
        show=False, save=True)

    starches.plotHeatMap(show=False, save=True)


if __name__ == '__main__':
    # path = "C:/Users/petrus.kirsten/PycharmProjects/Rheometer-Plotting/data/by sample"  # CEBB
    path = "C:/Users/Petrus Kirsten/Documents/GitHub/RheometerPlots/data/by sample"   # Personal

    # kappa(path)
    # iota(path)
    # blends(path)
    # starch(path)
    # starch_kappa(path)
    starch_iota(path)

    # plt.show()
