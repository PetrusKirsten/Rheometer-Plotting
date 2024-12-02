import numpy as np
import pandas as pd
from pathlib import Path

from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from matplotlib.ticker import MultipleLocator
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def fonts(folder_path, s=12, m=14):
    """Configures font properties for plots."""
    font_path = folder_path + 'HelveticaNeueThin.otf'
    helvetica_thin = FontProperties(fname=font_path)

    font_path = folder_path + 'HelveticaNeueLight.otf'
    helvetica_light = FontProperties(fname=font_path)

    font_path = folder_path + 'HelveticaNeueMedium.otf'
    helvetica_medium = FontProperties(fname=font_path)

    font_path = folder_path + 'HelveticaNeueBold.otf'
    helvetica_bold = FontProperties(fname=font_path)

    plt.rc('font', size=s)  # controls default text sizes
    plt.rc('axes', titlesize=s)  # fontsize of the axes title
    plt.rc('axes', labelsize=s)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=s)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=s)  # fontsize of the tick labels
    plt.rc('legend', fontsize=s)  # legend fontsize
    plt.rc('figure', titlesize=m)  # fontsize of the figure title


def powerLaw(omega, kPrime, nPrime):
    return kPrime * (omega ** nPrime)


def arraySplit(xArr, yArr, startValue, endValue):
    startIndex, endIndex = np.where(xArr >= startValue)[0][0], np.where(xArr <= endValue)[0][-1]

    return xArr[startIndex:endIndex], yArr[startIndex:endIndex]


def exportFit(
        sample,
        data, err,
        table
):
    keys = ("k'", "n'")
    values = (data, err)

    dictData = {'Sample': sample}
    iParams = 0
    for key, value in zip(keys, range(len(values[0]))):
        dictData[f'{key}'] = values[0][iParams]
        dictData[f'± {key}'] = values[1][iParams]
        iParams += 1

    table.append(dictData)

    return table


def getCteMean(values, tolerance=100):
    """
    :param values: to be analysed
    :param tolerance: the difference betweem two points data
    :return: the mean of the "cte" region and its indexes
    """
    diffs = np.abs(np.diff(values))  # Calcular as diferenças entre valores consecutivos

    constantRegions = diffs < tolerance  # Identificar regiões onde a diferença está abaixo do valor de tolerância

    # Encontrar os índices onde a condição é satisfeita
    iStart, iEnd = None, None
    lengthMax, lengthCurrent, currentStart = 0, 0, 0
    for i, is_constant in enumerate(constantRegions):
        if is_constant:
            if lengthCurrent == 0:
                currentStart = i
            lengthCurrent += 1
        else:
            if lengthCurrent > lengthMax:
                lengthMax = lengthCurrent
                iStart = currentStart
                iEnd = i
            lengthCurrent = 0

    if lengthCurrent > lengthMax:  # Checar se a última sequência é a maior constante
        iStart = currentStart
        iEnd = len(values) - 1

    if iStart is None or iEnd is None:  # Se nenhuma região constante foi encontrada
        return None, None, None

    mean = np.mean(values[iStart:iEnd + 1])  # Calcular a média da região constante encontrada
    stddev = np.std(values[iStart:iEnd + 1])  # Calcular a média da região constante encontrada

    mean = round(mean, -1)
    stddev = round(stddev, -1)

    return mean, stddev, iStart, iEnd


def getSamplesInfos(
        # quantity
        st_n, st_kc_n, st_ic_n,
        stCL_n, st_kcCL_n, st_icCL_n,
        kc_n, kcCL_n,
        # colors
        st_color, st_kc_color, st_ic_color,
        stCL_color, st_kcCL_color, st_icCL_color,
        kc_color, kcCL_color
):
    number_samples = [
        st_n, st_kc_n, st_ic_n,
        stCL_n, st_kcCL_n, st_icCL_n,
        kc_n, kcCL_n]

    colors_samples = [
        st_color, st_kc_color, st_ic_color,
        stCL_color, st_kcCL_color, st_icCL_color,
        kc_color, kcCL_color]

    return number_samples, colors_samples


def getSamplesData(
        dataPath,
        number_samples
):
    def getSegments(dataframe):
        freq = dataframe['f in Hz'].to_numpy()
        elastic = dataframe["G' in Pa"].to_numpy()
        loss = dataframe['G" in Pa'].to_numpy()

        # Identifying segments in the data
        seg2, seg3 = (
            dataframe.index[dataframe['SegIndex'] == seg].to_list()[0] for seg in ['2|1', '3|1'])

        # Slice segments
        segments = lambda arr: (arr[seg2:seg3])

        return {
            'freq': segments(freq),
            'storage': segments(elastic),
            'loss': segments(loss)
        }

    samples = {
        '0St': [], '0St + kCar': [], '0St + iCar': [],
        '0St/CL': [], '0St + kCar/CL': [], '0St + iCar/CL': [],
        'kCar': [], 'kCar/CL': []
    }
    sample_keys = list(samples.keys())
    sample_labels = (
            [sample_keys[0]] * number_samples[0] +
            [sample_keys[1]] * number_samples[1] +
            [sample_keys[2]] * number_samples[2] +
            [sample_keys[3]] * number_samples[3] +
            [sample_keys[4]] * number_samples[4] +
            [sample_keys[5]] * number_samples[5] +
            [sample_keys[6]] * number_samples[6] +
            [sample_keys[7]] * number_samples[7]
    )

    for sample_type, path in zip(sample_labels, dataPath):
        df = pd.read_excel(path)
        segments = getSegments(df)
        samples[sample_type].append(segments)

    dict_data = {}
    for sample_type in samples:
        dict_data[f'{sample_type}_freq'] = [s['freq'] for s in samples[sample_type]]
        dict_data[f'{sample_type}_storage'] = [s['storage'] for s in samples[sample_type]]
        dict_data[f'{sample_type}_loss'] = [s['loss'] for s in samples[sample_type]]

    return dict_data, sample_keys


def plotFreqSweeps(sampleName,
                   ax, x, yP, yD, yPerr, yDerr,
                   axTitle, yLabel, yLim, xLabel, xLim,
                   curveColor,
                   individualData=False, logScale=True,
                   startVal=0, endVal=31, tableData=None):
    def legendLabel():
        legend = ax.legend(
            loc='lower right', ncols=2,
            fancybox=False, frameon=True, framealpha=0.9, fontsize=9)
        legend.get_frame().set_facecolor('w')
        legend.get_frame().set_edgecolor('whitesmoke')

    def configPlot(idSample=0):
        dotCteMean = 'k'
        idSample = idSample + 1 if individualData else 'Mean'
        axisColor = '#303030'

        ax.set_title(axTitle, size=10, color='k')
        ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(0.75)
        ax.spines[['top', 'bottom', 'left', 'right']].set_color(axisColor)
        ax.tick_params(axis='both', which='both', colors=axisColor)

        ax.grid(True, which='major', axis='y', linestyle='-', linewidth=1, color='lightgray', alpha=0.5)
        ax.grid(True, which='minor', axis='y', linestyle='--', linewidth=.75, color='lightgray', alpha=0.5)

        ax.set_xlabel(f'{xLabel}', color=axisColor)
        ax.set_xscale('log' if logScale else 'linear')
        ax.set_xlim(xLim)

        ax.set_ylabel(f'{yLabel}', color=axisColor)
        ax.set_yscale('log' if logScale else 'linear')
        ax.set_ylim(yLim)
        # ax.errorbar(
        #     [x[indexStart_storage], x[indexEnd_storage]], [yP[indexStart_storage], yP[indexEnd_storage]], yerr=0,
        #     color=dotCteMean, alpha=0.75,
        #     fmt='.', markersize=4, mfc=dotCteMean, mec=dotCteMean, mew=1,
        #     capsize=0, lw=1, linestyle='',
        #     label=f'', zorder=4)

    configPlot()

    x_toFit, y_toFit = arraySplit(x, yP, startVal, endVal)
    params, covariance = curve_fit(powerLaw, x_toFit, y_toFit)  # p0=(y_mean[0], y_mean[-1], 100))
    errors = np.sqrt(np.diag(covariance))
    tableData = exportFit(
        f'{sampleName}',
        params, errors,
        tableData)

    ax.errorbar(
        x[:-2], yP[:-2], yPerr[:-2],
        color=curveColor, alpha=.65,
        fmt='D' if 'CL' in sampleName else 'o', markersize=4.5 if 'CL' in sampleName else 5.25,
        mfc=curveColor, mec=curveColor, mew=1,
        capsize=2, lw=1, linestyle='',
        label=f'{sampleName}', zorder=3)
    # label=f'{sampleName}_{idSample} | '
    #       + "$\overline{G'} \\approx$" + f'{meanStorage:.0f} ± {storageMeanErr:.0f} ' + '$Pa$',
    # ax.errorbar(
    #     x, yD, yDerr,
    #     color=curveColor, alpha=1,
    #     fmt=markerStyle, markersize=7, mfc='w', mec=curveColor, mew=0.75,
    #     capsize=3, lw=0.75, linestyle=':',
    #     zorder=3)

    # if axTitle == 'After breakage':
    # rectConfig = [(44, 0), xLim[-1] - 44, 10000]
    # rect = Rectangle(*rectConfig, linewidth=.75, edgecolor='#303030', facecolor='snow', alpha=1, zorder=1)
    # ax.add_patch(rect)

    if axTitle == 'Before breakage':
        legendLabel()

    return tableData


def midAxis(color, ax):
    ax[0].spines['right'].set_color(color)
    ax[1].spines['left'].set_color(color)
    ax[1].yaxis.tick_right()
    ax[1].yaxis.set_label_position('right')


def main(dataPath, fileName):
    fonts('C:/Users/petrus.kirsten/AppData/Local/Microsoft/Windows/Fonts/')
    plt.style.use('seaborn-v0_8-ticks')
    fig, axes = plt.subplots(
        figsize=(9, 7), ncols=1, nrows=1,
        gridspec_kw={'width_ratios': [1]}, facecolor='snow')
    axFreqS = axes

    fig.suptitle(f'Viscoelastic recovery by frequency sweeps assay.')
    yTitle, yLimits = f"Storage modulus $G'$ (Pa)", (10 ** (-2), 3 * 10 ** 4)
    xTitle, xLimits = f'Frequency (Hz)', (0.08, 100)

    nSamples, colorSamples = getSamplesInfos(
        2, 3, 3,
        3, 1, 3,
        3, 4,
        'lightgray', 'hotpink', 'deepskyblue',
        'darkgray', 'crimson', 'mediumblue',
        'mediumorchid', 'rebeccapurple')
    data, labels = getSamplesData(dataPath, nSamples)

    dataDict = {
        labels[0]: ([], [], []),
        labels[1]: ([], [], []),
        labels[2]: ([], [], []),
        labels[3]: ([], [], []),
        labels[4]: ([], [], []),
        labels[5]: ([], [], []),
        labels[6]: ([], [], []),
    }

    meanBefore, meanAfter = [], []
    meanBeforeErr, meanAfterErr = [], []
    dataFittingBef, dataFittingAft = [], []

    for key, (x, gP, gD) in dataDict.items():
        x.append(data[f'{key}_freq'])
        gP.append(data[f'{key}_storage'])
        gD.append(data[f'{key}_loss'])

    for key, color in zip(dataDict, colorSamples):
        gP, gD = np.mean(dataDict[key][1], axis=1)[0], np.mean(dataDict[key][2], axis=1)[0]
        gPerr, gDerr = np.std(dataDict[key][1], axis=1)[0], np.std(dataDict[key][2], axis=1)[0]

        meanStorage, storageMeanErr, fitStart, fitEnd = getCteMean(gP)
        meanBefore.append(meanStorage)
        meanBeforeErr.append(storageMeanErr)

        dataFittingBef = plotFreqSweeps(
            sampleName=key,
            ax=axFreqS, x=np.mean(dataDict[key][0], axis=1)[0],
            yP=gP, yD=gD, yPerr=gPerr, yDerr=gDerr,
            axTitle='Before breakage', yLabel=yTitle, yLim=yLimits, xLabel=xTitle, xLim=xLimits, curveColor=color,
            logScale=True, startVal=fitStart, endVal=fitEnd, tableData=dataFittingBef)

    plt.subplots_adjust(
        wspace=0.155,
        top=0.91, bottom=0.1,
        left=0.1, right=0.96)
    plt.show()

    dirSave = Path(*Path(filePath[0]).parts[:Path(filePath[0]).parts.index('data') + 1])
    fig.savefig(f'{dirSave}' + f'\\{fileName}' + '.png', facecolor='w', dpi=600)
    print(f'\n\n· Chart saved at\n{dirSave}.')


if __name__ == '__main__':
    folderPath = "C:/Users/petrus.kirsten/PycharmProjects/RheometerPlots/data"
    # folderPath = "C:/Users/Petrus Kirsten/Documents/GitHub/RheometerPlots/data"

    filePath = [

        # 0St
        folderPath + "/031024/10_0WSt/10_0WSt-viscRec_1.xlsx",
        folderPath + "/031024/10_0WSt/10_0WSt-viscRec_2.xlsx",

        # 0St + kCar
        folderPath + "/091024/10_0WSt_kCar/10_0WSt_kCar-viscoelasticRecovery-Flow_2a.xlsx",
        folderPath + "/091024/10_0WSt_kCar/10_0WSt_kCar-viscoelasticRecovery-Flow_3a.xlsx",
        folderPath + "/091024/10_0WSt_kCar/10_0WSt_kCar-viscoelasticRecovery-Flow_4a.xlsx",

        # 0St + iCar
        folderPath + "/031024/10_0WSt_iCar/10_0WSt_iCar-viscoRecoveryandFlow_2.xlsx",
        # folderPath + "10_0WSt_iCar/10_0WSt_iCar-viscoRecoveryandFlow_1.xlsx",
        folderPath + "/031024/10_0WSt_iCar/10_0WSt_iCar-viscoRecoveryandFlow_3.xlsx",
        folderPath + "/031024/10_0WSt_iCar/10_0WSt_iCar-viscoRecoveryandFlow_4.xlsx",

        # 0St/CL
        folderPath + "/171024/10_0St_CL/10_0St_CL-recovery-1.xlsx",
        folderPath + "/171024/10_0St_CL/10_0St_CL-recovery-2.xlsx",
        folderPath + "/171024/10_0St_CL/10_0St_CL-recovery-3.xlsx",

        # 0St+ kCar/CL
        folderPath + "/171024/10_0St_kC_CL/10_0St_kC_CL-recovery-1.xlsx",
        # folderPath + "/171024/10_0St_kC_CL/10_0St_CL-recovery-2.xlsx",

        # 0St + iCar/CL
        folderPath + "/171024/10_0St_iC_CL/10_0St_iC_CL-recovery-1.xlsx",
        folderPath + "/171024/10_0St_iC_CL/10_0St_iC_CL-recovery-2.xlsx",
        folderPath + "/171024/10_0St_iC_CL/10_0St_iC_CL-recovery-3.xlsx",

        # kC
        folderPath + "/231024/kC/kC-viscoelasticRecovery-1.xlsx",
        folderPath + "/231024/kC/kC-viscoelasticRecovery-2.xlsx",
        folderPath + "/231024/kC/kC-viscoelasticRecovery-3.xlsx",

        # kC/CL
        folderPath + "/231024/kC_CL/kC_CL-viscoelasticRecovery-1.xlsx",
        folderPath + "/231024/kC_CL/kC_CL-viscoelasticRecovery-2.xlsx",
        folderPath + "/231024/kC_CL/kC_CL-viscoelasticRecovery-3.xlsx",
        folderPath + "/231024/kC_CL/kC_CL-viscoelasticRecovery-4.xlsx",
    ]

    main(filePath, '0WSt_Car_CL-ViscoelasticRecovery')
