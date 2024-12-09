from rheology.plotting import Recovery


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
        filePath, 'kappas-teste',
        kappa_keySamples, kappa_nSamples, kappa_cSamples)

    kappas.plotGmodulus(
        f"Elastic modulus $G'$ (Pa)", (1 * 10 ** (-2), 1 * 10 ** 5),
        f'Frequency (Hz)', (.075, 100),
        show=False)

    kappas.plotBars(
        [.6, 10000, 850, 15],
        [0, None, None, None],
        show=False)

    kappas.plotHeatMap(show=True)


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
    nSamples = [2, 2, 2, 2, 3,]
    cSamples = ['greenyellow', '#80ed99', '#57cc99', '#38a3a5', '#22577a']

    iotas = Recovery(
        filePath, 'iotas-teste',
        keySamples, nSamples, cSamples)

    iotas.plotGmodulus(
        f"Elastic modulus $G'$ (Pa)", (1 * 10 ** (-2), 3 * 10 ** 3),
        f'Frequency (Hz)', (.075, 100),
        show=False)

    iotas.plotBars(
        [3.5, 45, 8.5, 7],
        [None, None, None, None],
        show=False)

    iotas.plotHeatMap(show=True)


def starch(folderPath):
    filePath = [
        # St
        folderPath + "/10St/10_0WSt-viscRec_1.xlsx",
        folderPath + "/10St/10_0WSt-viscRec_2.xlsx",

        # St + kCar
        folderPath + "/10St_kC/10_0WSt_kCar-viscoelasticRecovery-Flow_2a.xlsx",
        folderPath + "/10St_kC/10_0WSt_kCar-viscoelasticRecovery-Flow_3a.xlsx",
        # folderPath + "/10St_kC/10_0WSt_kCar-viscoelasticRecovery-Flow_4a.xlsx",

        # 0St + iCar
        folderPath + "/10St_iC/10_0WSt_iCar-viscoRecoveryandFlow_2.xlsx",
        # folderPath + "10St_iC/10_0WSt_iCar-viscoRecoveryandFlow_1.xlsx",
        folderPath + "/10St_iC/10_0WSt_iCar-viscoRecoveryandFlow_3.xlsx",
        folderPath + "/10St_iC/10_0WSt_iCar-viscoRecoveryandFlow_4.xlsx",

        # St/CL_7
        folderPath + "/10St_CL_7/10_0St_CL-recovery-1.xlsx",
        # folderPath + "/10St_CL_7/10_0St_CL-recovery-2.xlsx",
        folderPath + "/10St_CL_7/10_0St_CL-recovery-3.xlsx",

        # St + kCar/CL_7
        folderPath + "/10St_kC_CL_7/10_0St_kC_CL-recovery-1.xlsx",
        folderPath + "/10St_kC_CL_7/10_0St_kC_CL-recovery-3_off.xlsx",
        folderPath + "/10St_kC_CL_7/10_0St_kC_CL-recovery-4.xlsx",

        # St + iCar/CL_7
        folderPath + "/10St_iC_CL_7/10_0St_iC_CL-recovery-1.xlsx",
        folderPath + "/10St_iC_CL_7/10_0St_iC_CL-recovery-2.xlsx",
        folderPath + "/10St_iC_CL_7/10_0St_iC_CL-recovery-3.xlsx",
    ]

    keySamples = {
        'St': [],
        'St + kCar': [],
        'St + iCar': [],
        'St/CL_7': [],
        'St + kCar/CL_7': [],
        'St + iCar/CL_7': []
    }
    nSamples = [2, 2, 3,
                2, 3, 3,]
    cSamples = ['silver', 'hotpink', 'lightskyblue',
                'grey', 'mediumvioletred', 'royalblue']

    starches = Recovery(
        filePath, 'starches-teste',
        keySamples, nSamples, cSamples)

    starches.plotGmodulus(
        f"Elastic modulus $G'$ (Pa)", (1 * 10 ** 0, 1 * 10 ** 5),
        f'Frequency (Hz)', (.075, 100),
        show=False)

    starches.plotBars(
        [0.5, 1800, 500, 30],
        [None, 4, None, None],
        show=False)

    starches.plotHeatMap(show=True)


if __name__ == '__main__':
    path = "C:/Users/petrus.kirsten/PycharmProjects/Rheometer-Plotting/data/by sample"  # CEBB
    # path = "C:/Users/Petrus Kirsten/Documents/GitHub/RheometerPlots/data/by sample"  # Personal

    # kappa(path)
    # iota(path)
    starch(path)
