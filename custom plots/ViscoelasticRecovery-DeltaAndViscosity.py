import numpy as np
import pandas as pd
from pathlib import Path

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties


def fonts(folder_path, s=11, m=13):
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


def getCteMean(values, tolerance=50):
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

    return mean, stddev, iStart, iEnd


def getSamplesInfos(
        # quantity
        st_n, st_kc_n, st_ic_n,
        stCL_n, st_icCL_n,
        kc_n, kcCL_n,
        # colors
        st_color, st_kc_color, st_ic_color,
        stCL_color, st_icCL_color,
        kc_color, kcCL_color
):
    number_samples = [
        st_n, st_kc_n, st_ic_n,
        stCL_n, st_icCL_n,
        kc_n, kcCL_n]

    colors_samples = [
        st_color, st_kc_color, st_ic_color,
        stCL_color, st_icCL_color,
        kc_color, kcCL_color]

    return number_samples, colors_samples


def getSamplesData(dataPath, number_samples):

    def getSegments(dataframe):
        """
        Extracts freq, shear rate, shear stress, and delta segments from the dataframe.
        Returns tuples of constant and step segments.
        """
        freq = dataframe['f in Hz'].to_numpy()
        delta = dataframe['tan(δ) in -'].to_numpy()
        viscComplex = dataframe['|η*| in mPas'].to_numpy()

        # Identifying segments in the data
        seg2, seg3, seg5, seg6 = (
            dataframe.index[dataframe['SegIndex'] == seg].to_list()[0] for seg in ['2|1', '3|1', '5|1', '5|31'])

        # Slice segments
        segments = lambda arr: (arr[seg2:seg3], arr[seg5:seg6])  # Returns (constant segment, step segment)

        return {
            'freq': segments(freq),
            'delta': segments(delta),
            'viscosity': segments(viscComplex)
        }

    samples = {
        '0St': [], '0St + kCar': [], '0St + iCar': [],
        '0St/CL': [], '0St + iCar/CL': [],
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
            [sample_keys[6]] * number_samples[6]
    )

    for sample_type, path in zip(sample_labels, dataPath):
        df = pd.read_excel(path)
        segments = getSegments(df)
        samples[sample_type].append(segments)

    dict_freqSweeps = {}

    for sample_type in samples:
        dict_freqSweeps[f'{sample_type} freq'] = [s['freq'][0] for s in samples[sample_type]]
        dict_freqSweeps[f'{sample_type} delta'] = [s['delta'][0] for s in samples[sample_type]]
        dict_freqSweeps[f'{sample_type} viscosity'] = [s['viscosity'][0] for s in samples[sample_type]]

        dict_freqSweeps[f'{sample_type} freq_broken'] = [s['freq'][-1] for s in samples[sample_type]]
        dict_freqSweeps[f'{sample_type} delta_broken'] = [s['delta'][-1] for s in samples[sample_type]]
        dict_freqSweeps[f'{sample_type} viscosity_broken'] = [s['viscosity'][-1] for s in samples[sample_type]]

    return dict_freqSweeps, sample_keys


def plotFreqSweeps(sampleName, ax, axTitle,
                   x, yV, yPerr, yD, yDerr,
                   yLabel, yLim,
                   yLabel2, yLim2,
                   xLabel, xLim,
                   curveColor, logScale=True):

    def configPlot(axisColor='#383838'):
        axes_visc, axes_tan = ax[0], ax[1]

        axes_visc.set_title(axTitle, size=11, color='crimson')
        axes_visc.spines[['top', 'bottom', 'left', 'right']].set_linewidth(0.75)
        axes_visc.spines[['top', 'bottom', 'left', 'right']].set_color(axisColor)
        # ax.tick_params(axis='both', colors=axisColor)

        axes_visc.grid(True, which='both', axis='y', linestyle='-', linewidth=0.5, color='lightgrey', alpha=0.5)
        axes_tan.grid(True, which='both', axis='y', linestyle='-', linewidth=0.5, color='lightgrey', alpha=0.5)

        axes_tan.spines[['top', 'bottom', 'left', 'right']].set_linewidth(0.75)
        axes_visc.set_xlabel(f'{xLabel}', color=axisColor)
        axes_visc.set_xscale('log' if logScale else 'linear')
        axes_visc.set_xlim(xLim)

        axes_visc.set_ylabel(f'{yLabel}', color=axisColor)
        axes_visc.set_yscale('log' if logScale else 'linear')
        axes_visc.set_ylim(yLim)

        axes_tan.set_xlabel(f'{xLabel}', color=axisColor)
        axes_tan.set_xscale('log' if logScale else 'linear')
        axes_tan.set_xlim(xLim)

        axes_tan.set_ylabel(f'{yLabel2}', color=axisColor)
        axes_tan.set_yscale('log' if logScale else 'linear')
        axes_tan.set_ylim(yLim2)

        return axes_visc, axes_tan

    axVisc, axTan = configPlot()

    axVisc.errorbar(
        x, yV, yPerr,
        color=curveColor, alpha=.65,
        fmt='o', markersize=6, mfc=curveColor, mec=curveColor, mew=1.,
        capsize=3, lw=1, linestyle='',
        label=f'{sampleName}', zorder=3)

    axTan.errorbar(
        x, yD, yDerr,
        color=curveColor, alpha=.65,
        fmt='s', markersize=6, mfc=curveColor, mec=curveColor, mew=1,
        capsize=3, lw=1, linestyle='',
        zorder=3)


def midAxis(color, ax):
    # ax[0].spines['right'].set_color(color)
    # ax[1].spines['left'].set_color(color)
    ax[0, 1].yaxis.tick_right()
    ax[0, 1].yaxis.set_label_position('right')
    ax[1, 1].yaxis.tick_right()
    ax[1, 1].yaxis.set_label_position('right')

    ax[0, 0].tick_params(axis='x', which='both', direction='in')
    ax[0, 1].tick_params(axis='x', which='both', direction='in')


def legendLabel(ax):
    legend = ax[0, 1].legend(fancybox=False, frameon=True, framealpha=0.9, fontsize=9)
    legend.get_frame().set_facecolor('w')
    legend.get_frame().set_edgecolor('whitesmoke')


def main(dataPath):
    fonts('C:/Users/petrus.kirsten/AppData/Local/Microsoft/Windows/Fonts/')

    fileName = '0WSt_and_Car-ViscoelasticRecovery-DeltaAndViscosity'
    dirSave = Path(*Path(filePath[0]).parts[:Path(filePath[0]).parts.index('data') + 1])

    plt.style.use('seaborn-v0_8-ticks')
    fig, axes = plt.subplots(figsize=(18, 10), facecolor='w', ncols=2, nrows=2)
    fig.suptitle(f'Viscoelastic recovery by frequency sweeps assay.')

    yTitle, yLimits = f'Complex viscosity (mPa·s)', (10**1, 3*10**6)
    y2Title, y2Limits = f'tan(δ)', (3 * 10**(-2), 2*10**1)
    xTitle, xLimits = f'Frequency (Hz)', (0.06, 200)

    nSamples, colorSamples = getSamplesInfos(
        2, 3, 3,
        2, 3,
        3, 4,
        'lightgray', 'hotpink', 'deepskyblue',
        'darkgray', 'mediumblue',
        'mediumorchid', 'rebeccapurple')

    data, labels = getSamplesData(dataPath, nSamples)

    listBefore, listAfter = {
        labels[0]: ([], [], []),
        labels[1]: ([], [], []),
        labels[2]: ([], [], []),
        labels[3]: ([], [], []),
        labels[4]: ([], [], []),
        labels[5]: ([], [], []),
        labels[6]: ([], [], []),
    }, {
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

    for key, (x, tan_d, visc) in listBefore.items():
        x.append(data[f'{key} freq'])
        tan_d.append(data[f'{key} delta'])
        visc.append(data[f'{key} viscosity'])

    for key, (x, tan_d, visc) in listAfter.items():
        x.append(data[f'{key} freq_broken'])
        tan_d.append(data[f'{key} delta_broken'])
        visc.append(data[f'{key} viscosity_broken'])

    for k_a, k_b, c in zip(listAfter, listBefore, colorSamples):
        delta, visc = np.mean(listBefore[k_a][1], axis=1)[0], np.mean(listBefore[k_a][2], axis=1)[0]
        deltaErr, viscErr = np.std(listBefore[k_a][1], axis=1)[0], np.std(listBefore[k_a][2], axis=1)[0]

        meanStorage, storageMeanErr, _, _ = getCteMean(delta)
        meanBefore.append(meanStorage)
        meanBeforeErr.append(storageMeanErr)

        plotFreqSweeps(  # Before axes
            ax=axes[:, 0], sampleName=k_a,
            x=np.mean(listBefore[k_a][0], axis=1)[0],
            yV=visc, yPerr=viscErr,
            yD=delta, yDerr=deltaErr,
            axTitle='Before breakage',
            yLabel=yTitle, yLim=yLimits,
            yLabel2=y2Title, yLim2=y2Limits,
            xLabel=xTitle, xLim=xLimits,
            curveColor=c, logScale=True)

        delta, visc = np.mean(listAfter[k_a][1], axis=1)[0], np.mean(listAfter[k_a][2], axis=1)[0]
        deltaErr, viscErr = np.std(listAfter[k_a][1], axis=1)[0], np.std(listAfter[k_a][2], axis=1)[0]

        meanStorage, storageMeanErr, _, _ = getCteMean(delta)
        meanAfter.append(meanStorage)
        meanAfterErr.append(storageMeanErr)

        plotFreqSweeps(  # After axes
            ax=axes[:, 1], sampleName=k_a,
            x=np.mean(listAfter[k_a][0], axis=1)[0],
            yV=visc, yPerr=viscErr,
            yD=delta, yDerr=deltaErr,
            axTitle='After breakage',
            yLabel=yTitle, yLim=yLimits,
            yLabel2=y2Title, yLim2=y2Limits,
            xLabel=xTitle, xLim=xLimits,
            curveColor=c, logScale=True)

    legendLabel(axes)

    midAxis('#383838', axes)

    plt.subplots_adjust(wspace=0.0, hspace=0.0, top=0.92, bottom=0.08, left=0.05, right=0.95)
    plt.show()
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
        # folderPath + "/171024/10_0St_CL/10_0St_CL-recovery-2.xlsx",
        folderPath + "/171024/10_0St_CL/10_0St_CL-recovery-3.xlsx",

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

    main(dataPath=filePath)
