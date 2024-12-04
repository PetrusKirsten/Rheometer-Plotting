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
        st_kc_n,
        # colors
        st_kc_color
):
    number_samples = [st_kc_n]
    colors_samples = [st_kc_color]

    return number_samples, colors_samples


def getSamplesData(
        dataPath,
        number_samples
):
    def getSegments(dataframe):
        freq = dataframe['f in Hz'].to_numpy()
        elastic = dataframe["G' in Pa"].to_numpy()
        loss = dataframe['G" in Pa'].to_numpy()

        if 'off' in path:
            seg2, seg3 = (
                dataframe.index[dataframe['SegIndex'] == seg].to_list()[0] for seg in ['2|1', '3|1'])
            segments = lambda arr: (arr[seg2:seg3], arr[seg2:seg3])

        else:
            seg2, seg3, seg5, seg6 = (
                dataframe.index[dataframe['SegIndex'] == seg].to_list()[0] for seg in ['2|1', '3|1', '5|1', '5|31'])
            segments = lambda arr: (arr[seg2:seg3], arr[seg5:seg6+1])

        return {
            'freq': segments(freq),
            'storage': segments(elastic),
            'loss': segments(loss)
        }

    samples = {'0St + kCar': []}
    sample_keys = list(samples.keys())
    sample_labels = ([sample_keys[0]] * number_samples)

    for sample_type, path in zip(sample_labels, dataPath):
        df = pd.read_excel(path)
        segments = getSegments(df)
        samples[sample_type].append(segments)

    dict_data = {}
    for sample_type in samples:
        dict_data[f'{sample_type}_freq'] = [s['freq'][0] for s in samples[sample_type]]
        dict_data[f'{sample_type}_storage'] = [s['storage'][0] for s in samples[sample_type]]
        dict_data[f'{sample_type}_loss'] = [s['loss'][0] for s in samples[sample_type]]

        dict_data[f'{sample_type}_freq_broken'] = [s['freq'][-1] for s in samples[sample_type]]
        dict_data[f'{sample_type}_storage_broken'] = [s['storage'][-1] for s in samples[sample_type]]
        dict_data[f'{sample_type}_loss_broken'] = [s['loss'][-1] for s in samples[sample_type]]

    return dict_data, sample_keys


def plotFreqSweeps(sampleName,
                   ax, x, yP, yD, yPerr, yDerr,
                   axTitle, yLabel, yLim, xLabel, xLim,
                   curveColor,
                   individualData=False, logScale=True,
                   startVal=0, endVal=31, tableData=None):
    def legendLabel():
        legend = ax.legend(
            loc='best', fancybox=False, frameon=True,
            ncols=1, framealpha=0.9, fontsize=9)
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

        ax.grid(True, which='major', axis='y', linestyle='-', linewidth=.75, color='lightgrey', alpha=0.75)
        ax.grid(True, which='minor', axis='y', linestyle='-', linewidth=.50, color='lightgrey', alpha=0.5)

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

    ax.errorbar(
        x[:-2], yP[:-2], 0,
        color=curveColor, alpha=.95 - 0.35*sampleName,
        fmt='o',
        markersize=7,
        mfc=curveColor, mec='#383838', mew=.75,
        linestyle='',
        label=f'{sampleName+1}', zorder=3)

    if axTitle == 'After breakage':
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
    fig = plt.figure(figsize=(12, 6), facecolor='snow')
    gs = GridSpec(1, 2, width_ratios=[1.5, 1.5], height_ratios=[1])

    axPre, axPost = fig.add_subplot(gs[:, 0]), fig.add_subplot(gs[:, 1])

    fig.suptitle(f'Viscoelastic recovery by frequency sweeps assay.')
    yTitle, yLimits = f"Storage modulus $G'$ (Pa)", (2 * 10 ** 0, 2 * 10 ** 4)
    xTitle, xLimits = f'Frequency (Hz)', (0.08, 100)

    nSamples, colorSamples = getSamplesInfos(
        3, 'hotpink')
    nSamples, colorSamples = nSamples[0], colorSamples[0]

    data, labels = getSamplesData(dataPath, nSamples)

    listBefore, listAfter = {labels[0]: ([], [], [])}, {labels[0]: ([], [], [])}

    meanBefore, meanAfter = [], []
    meanBeforeErr, meanAfterErr = [], []
    dataFittingBef, dataFittingAft = [], []

    for key, (x, gP, gD) in listBefore.items():
        x.append(data[f'{key}_freq'])
        gP.append(data[f'{key}_storage'])
        gD.append(data[f'{key}_loss'])

    for key, (x, gP, gD) in listAfter.items():
        x.append(data[f'{key}_freq_broken'])
        gP.append(data[f'{key}_storage_broken'])
        gD.append(data[f'{key}_loss_broken'])

    gP_bef, gD_bef = listBefore[labels[0]][1][0], listBefore[labels[0]][2][0]
    gP_aft, gD_aft = listAfter[labels[0]][1][0], listAfter[labels[0]][2][0]
    for key in range(len(gP_bef)):
        gP_aft_id, gP_bef_id = gP_aft[key], gP_bef[key]
        meanStorage, storageMeanErr, fitStart, fitEnd = getCteMean(gP_bef_id)
        meanBefore.append(meanStorage)
        meanBeforeErr.append(storageMeanErr)

        dataFittingBef = plotFreqSweeps(  # Before axes
            sampleName=key,
            ax=axPre, x=listBefore[labels[0]][0][0][0],
            yP=gP_bef_id, yD=0, yPerr=0, yDerr=0,
            axTitle='Before breakage', yLabel=yTitle, yLim=yLimits, xLabel=xTitle, xLim=xLimits, curveColor=colorSamples,
            logScale=True, startVal=fitStart, endVal=fitEnd, tableData=dataFittingBef)

        meanStorage, storageMeanErr, _, _ = getCteMean(gP_aft_id)
        meanAfter.append(meanStorage)
        meanAfterErr.append(storageMeanErr)

        dataFittingAft = plotFreqSweeps(  # After axes
            sampleName=key,
            ax=axPost, x=listBefore[labels[0]][0][0][0],
            yP=gP_aft_id, yD=0, yPerr=0, yDerr=0,
            axTitle='After breakage', yLabel=yTitle, yLim=yLimits, xLabel=xTitle, xLim=xLimits, curveColor=colorSamples,
            logScale=True, startVal=fitStart, endVal=fitEnd, tableData=dataFittingAft)

    plt.subplots_adjust(
        wspace=0.164, hspace=0.24,
        top=0.91, bottom=0.1,
        left=0.075, right=0.99)
    plt.show()

    dirSave = Path(*Path(filePath[0]).parts[:Path(filePath[0]).parts.index('data') + 1])
    fig.savefig(f'{dirSave}' + f'\\{fileName}' + '.png', facecolor='w', dpi=600)
    print(f'\n\n· Chart saved at\n{dirSave}.')


if __name__ == '__main__':
    # folderPath = "C:/Users/petrus.kirsten/PycharmProjects/RheometerPlots/data"
    folderPath = "C:/Users/Petrus Kirsten/Documents/GitHub/RheometerPlots/data"

    filePath = [
        # 0St + kCar TODO: idenfity outlier
        folderPath + "/091024/10_0WSt_kCar/10_0WSt_kCar-viscoelasticRecovery-Flow_2a.xlsx",
        folderPath + "/091024/10_0WSt_kCar/10_0WSt_kCar-viscoelasticRecovery-Flow_3a.xlsx",
        folderPath + "/091024/10_0WSt_kCar/10_0WSt_kCar-viscoelasticRecovery-Flow_4a.xlsx",
    ]

    main(filePath, 'Outliers-ViscoelasticRecovery')
