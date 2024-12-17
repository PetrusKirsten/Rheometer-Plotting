from matplotlib import pyplot as plt
from rheology.plotting import DynamicCompression
from rheology.plotting import BreakageCompression


def starch(folderPath):
    filePath = [
        # St
        folderPath + "/10St/10St-compression-1.xlsx",
        folderPath + "/10St/10St-compression-2.xlsx",
        folderPath + "/10St/10St-compression-3.xlsx",
        folderPath + "/10St/10St-compression-4.xlsx",

        # St/CL_7
        # folderPath + "/10St_CL_7/10_0St_CaCl2.xlsx",
        folderPath + "/10St_CL_7/10_0St_CL-compression-1_seg.xlsx",
        folderPath + "/10St_CL_7/10_0St_CL-compression-2_seg.xlsx",
        folderPath + "/10St_CL_7/10_0St_CL-compression-3_seg.xlsx",

        # St/CL_14
        folderPath + "/10St_CL_14/St_CL_14-compression-1.xlsx",
        folderPath + "/10St_CL_14/St_CL_14-compression-2.xlsx",
        folderPath + "/10St_CL_14/St_CL_14-compression-3.xlsx",
        folderPath + "/10St_CL_14/St_CL_14-compression-4.xlsx",
        folderPath + "/10St_CL_14/St_CL_14-compression-5.xlsx",
        folderPath + "/10St_CL_14/St_CL_14-compression-7.xlsx",

        # St/CL_28
        folderPath + "/10St_CL_28/St_CL_28-compression-1.xlsx",
        folderPath + "/10St_CL_28/St_CL_28-compression-2.xlsx",
        folderPath + "/10St_CL_28/St_CL_28-compression-3.xlsx",
    ]

    keySamples = {
        # No CL
        'St': [],
        # CL 7
        'St/CL_7': [],
        # CL 14
        'St/CL_14': [],
        # CL 28
        'St/CL_28': [],
    }
    nSamples = [
        # No CL
        4,
        # CL 7
        3,
        # CL 14
        6,
        # CL 28
        3,
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

    starches = DynamicCompression(
        filePath, 'starches',
        keySamples, nSamples, cSamples)

    starches.plotGraphs(
        150,
        150,
        [180, 28],
        show=False
    )


def starch_kappa(folderPath):
    filePath = [
        # St
        folderPath + "/10St_kC/10St_kC-compression-1.xlsx",
        folderPath + "/10St_kC/10St_kC-compression-2.xlsx",
        folderPath + "/10St_kC/10St_kC-compression-3.xlsx",
        folderPath + "/10St_kC/10St_kC-compression-4.xlsx",

        # St/CL_7
        folderPath + "/10St_kC_CL_7/10_0St_kC_CL-compression-1_seg.xlsx",
        folderPath + "/10St_kC_CL_7/10_0St_kC_CL-compression-2_seg.xlsx",

        # St/CL_14
        folderPath + "/10St_kC_CL_14/St_kC_CL_14-compression-1.xlsx",
        folderPath + "/10St_kC_CL_14/St_kC_CL_14-compression-3.xlsx",
        # folderPath + "/10St_kC_CL_14/St_kC_CL_14-compression-4.xlsx",

        # St/CL_28
        folderPath + "/10St_kC_CL_28/St_kC_CL_28-compression-1.xlsx",
        folderPath + "/10St_kC_CL_28/St_kC_CL_28-compression-3.xlsx",
        folderPath + "/10St_kC_CL_28/St_kC_CL_28-compression-7.xlsx",
        folderPath + "/10St_kC_CL_28/St_kC_CL_28-compression-8.xlsx",
    ]

    keySamples = {
        # No CL
        'St + kCar': [],
        # CL 7
        'St + kCar/CL_7': [],
        # CL 14
        'St + kCar/CL_14': [],
        # CL 28
        'St + kCar/CL_28': [],
    }
    nSamples = [
        # No CL
        4,
        # CL 7
        2,
        # CL 14
        2,
        # CL 28
        4,
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

    starches_kappa = DynamicCompression(
        filePath, 'starches_kappa',
        keySamples, nSamples, cSamples)

    starches_kappa.plotGraphs(
        200,
        200,
        [300, 25],
        show=False
    )


def starch_iota(folderPath):
    filePath = [
        # St + iCar
        folderPath + "/10St_iC/10St_kC-compression-1.xlsx",
        folderPath + "/10St_iC/10St_kC-compression-2.xlsx",
        folderPath + "/10St_iC/10St_kC-compression-3.xlsx",
        folderPath + "/10St_iC/10St_kC-compression-4.xlsx",

        # St + iCar/CL_7
        folderPath + "/10St_iC_CL_7/10_0St_iC_CL-compression-1_seg.xlsx",
        folderPath + "/10St_iC_CL_7/10_0St_iC_CL-compression-2b_seg.xlsx",
        folderPath + "/10St_iC_CL_7/10_0St_iC_CL-compression-3_seg.xlsx",
        folderPath + "/10St_iC_CL_7/10_0St_iC_CL-compression-4_seg.xlsx",
        folderPath + "/10St_iC_CL_7/10_0St_iC_CL-compression-5_seg.xlsx",

        # St + iCar/CL_14
        folderPath + "/10St_iC_CL_14/0St_iC_CL_14-compression-1.xlsx",
        folderPath + "/10St_iC_CL_14/0St_iC_CL_14-compression-2.xlsx",
        folderPath + "/10St_iC_CL_14/0St_iC_CL_14-compression-3.xlsx",
        folderPath + "/10St_iC_CL_14/0St_iC_CL_14-compression-4.xlsx",
        folderPath + "/10St_iC_CL_14/0St_iC_CL_14-compression-5.xlsx",

        # St + iCar/CL_28
        folderPath + "/10St_iC_CL_28/St_iC_CL_28-compression-1.xlsx",
        folderPath + "/10St_iC_CL_28/St_iC_CL_28-compression-2.xlsx",
        folderPath + "/10St_iC_CL_28/St_iC_CL_28-compression-3.xlsx",

    ]

    keySamples = {
        # No CL
        'St + iCar': [],
        # CL 7
        'St + iCar/CL_7': [],
        # CL 14
        'St + iCar/CL_14': [],
        # CL 28
        'St + iCar/CL_28': []
    }
    nSamples = [
        # No CL
        4,
        # CL 7
        5,
        # CL 14
        5,
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
        '#00934E',
    ]

    starches_iota = DynamicCompression(
        filePath, 'starches_iota',
        keySamples, nSamples, cSamples)

    starches_iota.plotGraphs(
        100,
        100,
        [150, 20],
        show=False
    )


if __name__ == '__main__':
    path = "C:/Users/petrus.kirsten/PycharmProjects/Rheometer-Plotting/data/by sample"  # CEBB
    # path = "C:/Users/Petrus Kirsten/Documents/GitHub/RheometerPlots/data/by sample"  # Personal

    starch(path)
    starch_kappa(path)
    starch_iota(path)

    # blends(path)

    plt.show()