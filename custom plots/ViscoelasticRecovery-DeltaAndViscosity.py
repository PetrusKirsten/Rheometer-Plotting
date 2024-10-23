import numpy as np
import pandas as pd
from pathlib import Path

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties


def fonts(folder_path, s=10, m=12):
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


def getSamplesData(dataPath, nSt, nKc, nIc, nStCL, nKcCL, nIcCL):

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

    # Store data for each sample type
    samples = {'10%-0St': [], '10%-0St + kCar': [], '10%-0St + iCar': [],
               '10%-0St/Ca²⁺': [], '10%-0St + iCar/Ca²⁺': [], '10%-0St+kCar/Ca²⁺': []}

    # Determine sample types for each path
    sample_labels = (['10%-0St'] * nSt + ['10%-0St + kCar'] * nKc + ['10%-0St + iCar'] * nIc +
                     ['10%-0St/Ca²⁺'] * nStCL + ['10%-0St + iCar/Ca²⁺'] * nKcCL + ['10%-0St+kCar/Ca²⁺'] * nIcCL)

    # Read data and categorize based on sample type
    for sample_type, path in zip(sample_labels, dataPath):
        df = pd.read_excel(path)
        segments = getSegments(df)
        samples[sample_type].append(segments)

    # Initialize dictionaries to hold the results
    dict_freqSweeps = {}

    # Populate dictionaries with consolidated sample data
    for sample_type in samples:
        dict_freqSweeps[f'{sample_type}_freq'] = [s['freq'][0] for s in samples[sample_type]]
        dict_freqSweeps[f'{sample_type}_delta'] = [s['delta'][0] for s in samples[sample_type]]
        dict_freqSweeps[f'{sample_type}_viscosity'] = [s['viscosity'][0] for s in samples[sample_type]]

        dict_freqSweeps[f'{sample_type}_freq_broken'] = [s['freq'][-1] for s in samples[sample_type]]
        dict_freqSweeps[f'{sample_type}_delta_broken'] = [s['delta'][-1] for s in samples[sample_type]]
        dict_freqSweeps[f'{sample_type}_viscosity_broken'] = [s['viscosity'][-1] for s in samples[sample_type]]

    return dict_freqSweeps


def plotFreqSweeps(sampleName, ax, axTitle,
                   x, yV, yPerr, yD, yDerr,
                   yLabel, yLim,
                   yLabel2, yLim2,
                   xLabel, xLim,
                   curveColor, logScale=True):

    def configPlot(axisColor='#383838'):
        axes_visc, axes_tan = ax[0], ax[1]

        axes_visc.set_title(axTitle, size=9, color='crimson')
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
        color=curveColor, alpha=.8,
        fmt='o', markersize=7, mfc=curveColor, mec=curveColor, mew=1.,
        capsize=3, lw=1, linestyle='',
        label=f'{sampleName}', zorder=3)

    axTan.errorbar(
        x, yD, yDerr,
        color=curveColor, alpha=.8,
        fmt='s', markersize=7, mfc=curveColor, mec=curveColor, mew=1,
        capsize=3, lw=1, linestyle='',
        zorder=3)


def plotInsetMean(data, dataErr, keys, colors, ax):
    ax_inset = inset_axes(ax, width='30%', height='20%', loc='upper left')  # Create an inset plot

    xInset = np.arange(len(data))
    ax_inset.barh(xInset, width=data, xerr=0,
                  color=colors, edgecolor='black', alpha=0.85, linewidth=0.5)

    ax_inset.errorbar(y=xInset, x=data, xerr=dataErr, alpha=1,
                      color='#383838', linestyle='', capsize=2, linewidth=0.75)

    for i in range(len(data)):
        ax_inset.text(data[i] + dataErr[i] + 100, xInset[i],
                      f'{data[i]:.0f} ± {dataErr[i]:.0f}',
                      size=8, va='center_baseline', ha='left', color='black')

    ax_inset.text(
        0.5, -0.08, "Average G' values (Pa)",
        ha='center', va='top', fontsize=9, transform=ax_inset.transAxes)
    ax_inset.set_facecolor('snow')  # Change to your desired color
    ax_inset.tick_params(axis='both', labelsize=8, length=0)
    ax_inset.spines[['top', 'bottom', 'left', 'right']].set_linewidth(0.75)
    ax_inset.spines[['top', 'bottom', 'left', 'right']].set_color('dimgrey')

    ax_inset.set_yticks(np.arange(len(data)))
    ax_inset.set_yticklabels(keys)
    ax_inset.yaxis.tick_right()
    ax_inset.yaxis.set_label_position('right')

    ax_inset.set_xticks([])
    ax_inset.set_xlim(0, 2900)


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

    yTitle, yLimits = f'Complex viscosity (mPa·s)', (10**3, 4*10**6)
    y2Title, y2Limits = f'tan(δ)', (1 * 10**(-2), 10**1)
    xTitle, xLimits = f'Frequency (Hz)', (0.06, 200)

    (st_n, kc_n, ic_n,
     stCL_n, kcCL_n, icCL_n) = 2, 3, 2, 0, 0, 0
    nSamples = (st_n, kc_n, ic_n,
                stCL_n, kcCL_n, icCL_n)

    (st_color, kc_color, ic_color,
     stCL_color, kcCL_color, icCL_color) = ('sandybrown', 'hotpink', 'deepskyblue',
                                            'grey', 'sandybrown', 'deepskyblue')
    colors = (st_color, kc_color, ic_color,
              stCL_color, kcCL_color, icCL_color)

    data = getSamplesData(dataPath, *nSamples)

    listBefore = {
        '10%-0St': ([], [], []),  # (x_5st, gP_5st, gD_5st)
        '10%-0St + kCar': ([], [], []),  # (x_10st, gP_10st, gD_10st)
        '10%-0St + iCar': ([], [], []),  # (x_ic, gP_ic, gD_ic)
        # '10%-0St/Ca²⁺': ([], [], []),  # (x_kc, gP_kc, gD_kc)
        # '10%-0St + iCar/Ca²⁺': ([], [], []),  # (x_kc, gP_kc, gD_kc)
        # '10%-0St+kCar/Ca²⁺': ([], [], [])  # (x_kc, gP_kc, gD_kc)
    }
    listAfter = {
        '10%-0St': ([], [], []),  # (x_5st, gP_5st, gD_5st)
        '10%-0St + kCar': ([], [], []),  # (x_10st, gP_10st, gD_10st)
        '10%-0St + iCar': ([], [], []),  # (x_ic, gP_ic, gD_ic)
        # '10%-0St/Ca²⁺': ([], [], []),  # (x_kc, gP_kc, gD_kc)
        # '10%-0St + iCar/Ca²⁺': ([], [], []),  # (x_kc, gP_kc, gD_kc)
        # '10%-0St+kCar/Ca²⁺': ([], [], [])  # (x_kc, gP_kc, gD_kc)
    }

    meanBefore, meanAfter = [], []
    meanBeforeErr, meanAfterErr = [], []

    for key, (x, tan_d, visc) in listBefore.items():
        x.append(data[f'{key}_freq'])
        tan_d.append(data[f'{key}_delta'])
        visc.append(data[f'{key}_viscosity'])

    for key, (x, tan_d, visc) in listAfter.items():
        x.append(data[f'{key}_freq_broken'])
        tan_d.append(data[f'{key}_delta_broken'])
        visc.append(data[f'{key}_viscosity_broken'])

    for k_a, k_b, c in zip(listAfter, listBefore, colors):
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

    # Inset plot (bar plot showing meanStorage)
    # plotInsetMean(
    #     data=meanBefore, dataErr=meanBeforeErr,
    #     keys=listBefore.keys(), colors=colors, ax=axes[0])
    # plotInsetMean(
    #     data=meanAfter, dataErr=meanAfterErr,
    #     keys=listBefore.keys(), colors=colors, ax=axes[1])

    plt.subplots_adjust(wspace=0.0, hspace=0.0, top=0.92, bottom=0.08, left=0.05, right=0.95)
    plt.show()
    fig.savefig(f'{dirSave}' + f'\\{fileName}' + '.png', facecolor='w', dpi=600)

    print(f'\n\n· Chart saved at\n{dirSave}.')


if __name__ == '__main__':
    folderPath = "C:/Users/petrus.kirsten/PycharmProjects/RheometerPlots/data"
    # folderPath = "C:/Users/Petrus Kirsten/Documents/GitHub/RheometerPlots/data"
    filePath = [
        folderPath + "/031024/10_0WSt/10_0WSt-viscRec_1.xlsx",
        folderPath + "/031024/10_0WSt/10_0WSt-viscRec_2.xlsx",
        #
        # folderPath + "/091024/10_0WSt_kCar/10_0WSt_kCar-viscoelasticRecovery-Flow_2a.xlsx",
        folderPath + "/091024/10_0WSt_kCar/10_0WSt_kCar-viscoelasticRecovery-Flow_3a.xlsx",
        folderPath + "/091024/10_0WSt_kCar/10_0WSt_kCar-viscoelasticRecovery-Flow_4a.xlsx",
        #
        # folderPath + "10_0WSt_iCar/10_0WSt_iCar-viscoRecoveryandFlow_1.xlsx",
        folderPath + "/031024/10_0WSt_iCar/10_0WSt_iCar-viscoRecoveryandFlow_2.xlsx",
        folderPath + "/031024/10_0WSt_iCar/10_0WSt_iCar-viscoRecoveryandFlow_3.xlsx",
        folderPath + "/031024/10_0WSt_iCar/10_0WSt_iCar-viscoRecoveryandFlow_4.xlsx",
    ]

    main(dataPath=filePath)
