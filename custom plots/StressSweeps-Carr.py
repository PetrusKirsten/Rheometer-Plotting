import numpy as np
import pandas as pd
from pathlib import Path
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator
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
        n_kc_0,
        n_ic_0,
        n_ic_14, n_ic_21, n_ic_28, n_ic_42,
        # colors
        color_kc_0,
        color_ic_0,
        color_ic_14, color_ic_21, color_ic_28, color_ic_42
):
    number_samples = [n_kc_0, n_ic_0, n_ic_14, n_ic_21, n_ic_28, n_ic_42]

    colors_samples = [color_kc_0, color_ic_0, color_ic_14, color_ic_21, color_ic_28, color_ic_42]

    return number_samples, colors_samples


def getSamplesData(
        dataPath,
        number_samples
):
    def getSegments(dataframe, str_seg3):
        stress = dataframe['τ in Pa'].to_numpy()
        elastic = dataframe["G' in Pa"].to_numpy()
        loss = dataframe['G" in Pa'].to_numpy()

        seg2, seg3 = (dataframe.index[dataframe['SegIndex'] == seg].to_list() for seg in ['2|1', str_seg3])
        segments = lambda arr: (arr[seg2[0]:seg3[0] + 1])

        return {
            'stress': segments(stress),
            'storage': segments(elastic),
            'loss': segments(loss)
        }

    samples = {
        'kCar': [],
        'iCar': [],
        'iCar/CL-14': [], 'iCar/CL-21': [], 'iCar/CL-28': [], 'iCar/CL-42': []
    }
    sample_keys = list(samples.keys())
    sample_labels = (
            [sample_keys[0]] * number_samples[0] +
            [sample_keys[1]] * number_samples[1] +
            [sample_keys[2]] * number_samples[2] +
            [sample_keys[3]] * number_samples[3] +
            [sample_keys[4]] * number_samples[4] +
            [sample_keys[5]] * number_samples[5])

    for sample_type, path in zip(sample_labels, dataPath):
        df = pd.read_excel(path)
        str_seg = '2|25' if 'CL' not in path else '2|40'
        segments = getSegments(df, str_seg)
        samples[sample_type].append(segments)

    dict_data = {}
    for sample_type in samples:
        dict_data[f'{sample_type} stress'] = [s['stress'] for s in samples[sample_type]]
        dict_data[f'{sample_type} storage'] = [s['storage'] for s in samples[sample_type]]
        dict_data[f'{sample_type} loss'] = [s['loss'] for s in samples[sample_type]]

    return dict_data, sample_keys


def plotStressSweeps(sampleName,
                     ax, x, yP, yD,
                     axTitle, yLabel, yLim, xLabel, xLim,
                     curveColor,
                     logScale=True):
    def legendLabel():
        """Applies consistent styling to legends in plots."""
        legend = ax.legend(fancybox=False, frameon=True, framealpha=0.9, fontsize=9)
        legend.get_frame().set_facecolor('w')
        legend.get_frame().set_edgecolor('whitesmoke')

    def configPlot():
        axisColor = '#303030'

        ax.set_title(axTitle, size=10, color='k')
        ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(0.75)
        ax.spines[['top', 'bottom', 'left', 'right']].set_color(axisColor)
        ax.tick_params(axis='both', which='both', colors=axisColor)

        ax.grid(True, which='major', axis='y', linestyle='-', linewidth=.75, color='lightgray', alpha=0.5)
        ax.grid(True, which='minor', axis='y', linestyle='-', linewidth=.5, color='lightgray', alpha=0.5)

        ax.set_xlabel(f'{xLabel}', color=axisColor)
        ax.set_xscale('log' if logScale else 'linear')
        ax.set_xlim(xLim)

        ax.set_ylabel(f'{yLabel}', color=axisColor)
        ax.set_yscale('log' if logScale else 'linear')
        ax.set_ylim(yLim)

    configPlot()

    ax.errorbar(
        x, yP, 0,
        color=curveColor, alpha=.65,
        fmt='o', markersize=7,
        mfc=curveColor, mec='#383838', mew=.75,
        linestyle='',
        label=f'{sampleName}', zorder=3)

    ax.errorbar(
        x, yD, 0,
        color=curveColor, alpha=.65,
        fmt='o', markersize=7,
        mfc='w', mec=curveColor, mew=1.25,
        linestyle='',
        label=f'', zorder=2)

    legendLabel()


def main(dataPath, fileName):
    fonts('C:/Users/petrus.kirsten/AppData/Local/Microsoft/Windows/Fonts/')
    plt.style.use('seaborn-v0_8-ticks')
    fig = plt.figure(figsize=(9, 7), facecolor='snow')
    gs = GridSpec(1, 1, width_ratios=[1], height_ratios=[1])
    axes = fig.add_subplot(gs[0, 0])

    fig.suptitle(f'Stress sweeps to linear viscoelastic region evaluation.')
    yTitle = f"Storage (●) and loss moduli (○) (Pa)"
    yLimits = (1 * 10 ** (-2), 1 * 10 ** 2)
    xTitle, xLimits = f'Shear stress (Pa)', (.1, 1175)

    nSamples, colorSamples = getSamplesInfos(
        1,
        1,
        1, 1, 1, 1,
        '#fb7e8f',
        'greenyellow',
        '#80ed99', '#57cc99', '#38a3a5', '#22577a')

    data, labels = getSamplesData(dataPath, nSamples)

    dictData = {
        labels[0]: ([], [], []),
        labels[1]: ([], [], []),
        labels[2]: ([], [], []),
        labels[3]: ([], [], []),
        labels[4]: ([], [], []),
        labels[5]: ([], [], [])
    }

    for key, (x, gP, gD) in dictData.items():
        x.append(data[f'{key} stress'])
        gP.append(data[f'{key} storage'])
        gD.append(data[f'{key} loss'])

    for key, color in zip(dictData, colorSamples):
        tau, gP, gD = dictData[key][0][0][0], dictData[key][1][0][0], dictData[key][2][0][0]

        plotStressSweeps(
            sampleName=key, axTitle='',
            ax=axes, x=tau, yP=gP, yD=gD,
            yLabel=yTitle, yLim=yLimits, xLabel=xTitle, xLim=xLimits,
            curveColor=color,
            logScale=True)

    plt.subplots_adjust(
        wspace=0.164, hspace=0.24,
        top=0.94, bottom=0.1,
        left=0.1, right=0.98)
    plt.show()

    dirSave = Path(*Path(filePath[0]).parts[:Path(filePath[0]).parts.index('data') + 1])
    fig.savefig(f'{dirSave}' + f'\\{fileName}' + '.png', facecolor='w', dpi=600)
    print(f'\n\n· Chart saved at\n{dirSave}.')


if __name__ == '__main__':
    folderPath = "C:/Users/petrus.kirsten/PycharmProjects/RheometerPlots/data"
    # folderPath = "C:/Users/Petrus Kirsten/Documents/GitHub/RheometerPlots/data"

    filePath = [
        # kC
        folderPath + "/231024/kC/kC-stressSweeps.xlsx",

        # iC
        folderPath + "/231024/iC/iC-stressSweeps-1.xlsx",

        # iC CL 14
        folderPath + "/311024/iC_CL_14/iC_CL_14-stressSweep-1.xlsx",

        # iC CL 21
        folderPath + "/311024/iC_CL_21/iC_CL_21-stressSweep-1.xlsx",

        # iC CL 28
        folderPath + "/311024/iC_CL_28/iC_CL_28-stressSweep-1.xlsx",

        # iC CL 42,
        folderPath + "/311024/iC_CL_42/iC_CL_42-stressSweep-1.xlsx",
    ]

    main(filePath, 'Car-StressSweeps')
