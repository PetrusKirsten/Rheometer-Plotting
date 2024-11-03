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
        n_kc_0, n_kc_14,
        n_ic_14, n_ic_21, n_ic_28, n_ic_42,
        # colors
        color_kc_0, color_kc_14,
        color_ic_14, color_ic_21, color_ic_28, color_ic_42
):
    number_samples = [n_kc_0, n_kc_14, n_ic_14, n_ic_21, n_ic_28, n_ic_42]

    colors_samples = [color_kc_0, color_kc_14, color_ic_14, color_ic_21, color_ic_28, color_ic_42]

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
        seg2, seg3, seg5, seg6 = (
            dataframe.index[dataframe['SegIndex'] == seg].to_list()[0] for seg in ['2|1', '3|1', '5|1', '5|31'])

        # Slice segments
        segments = lambda arr: (arr[seg2:seg3], arr[seg5:seg6])  # Returns (constant segment, step segment)

        return {
            'freq': segments(freq),
            'storage': segments(elastic),
            'loss': segments(loss)
        }

    samples = {
        'kCar': [], 'kCar/CL-14': [],
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
        """Applies consistent styling to legends in plots."""
        legend = ax.legend(fancybox=False, frameon=True, framealpha=0.9, fontsize=9)
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

        ax.grid(True, which='major', axis='y', linestyle='-', linewidth=.75, color='lightgray', alpha=0.5)
        ax.grid(True, which='minor', axis='y', linestyle='-', linewidth=.5, color='lightgray', alpha=0.5)

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
        color=curveColor, alpha=.85,
        fmt='none', mfc=curveColor,
        capsize=2.5, capthick=1, linestyle='', lw=1,
        label=f'', zorder=2)

    ax.errorbar(
        x[:-2], yP[:-2], 0,
        color=curveColor, alpha=.65,
        fmt='o', markersize=5.4,
        mfc=curveColor, mec='#383838', mew=.75,
        linestyle='',
        label=f'{sampleName}', zorder=3)

    # if axTitle == 'After breakage':
    # rectConfig = [(44, 0), xLim[-1] - 44, 10000]
    # rect = Rectangle(*rectConfig, linewidth=.75, edgecolor='#303030', facecolor='snow', alpha=1, zorder=1)
    # ax.add_patch(rect)

    if axTitle == 'After breakage':
        legendLabel()

    return tableData


def plotInset(data, dataErr, keys, colors, ax, recovery=None):
    ax_inset = inset_axes(ax, width='40%', height='25%', loc='lower right')

    xInset = np.arange(len(data))
    ax_inset.barh(xInset, width=data, xerr=0,
                  color=colors, edgecolor='black', alpha=.8, linewidth=0.5)

    ax_inset.errorbar(y=xInset, x=data, xerr=dataErr, alpha=1,
                      color='#303030', linestyle='', capsize=2, linewidth=0.75)

    for i in range(len(data)):
        ax_inset.text(data[i] + dataErr[i] + 100, xInset[i],
                      f'{data[i]:.0f} ± {dataErr[i]:.0f} '
                      f'~ {100 * data[i] / recovery[i]:.0f}%'
                      if recovery is not None
                      else f'{data[i]:.0f} ± {dataErr[i]:.0f}',
                      size=8, va='center_baseline', ha='left', color='black')

    ax_inset.text(
        0.5, 1.1, "Average G' values (Pa)",
        ha='center', va='top', fontsize=9, transform=ax_inset.transAxes)
    ax_inset.set_facecolor('snow')
    ax_inset.tick_params(axis='both', labelsize=8, length=0)
    ax_inset.spines[['top', 'bottom', 'left', 'right']].set_linewidth(0.75)
    ax_inset.spines[['top', 'bottom', 'left', 'right']].set_color('dimgrey')

    ax_inset.set_yticks(np.arange(len(data)))
    ax_inset.set_yticklabels(keys)
    # ax_inset.yaxis.tick_right()
    # ax_inset.yaxis.set_label_position('right')

    ax_inset.set_xticks([])
    ax_inset.set_xlim(0, 2300)

    return data, dataErr


def plotBars(
        title, axes, lim,
        data_before, data_after,
        colors, flt,
        scale_correction=None, a=.9, h='', z=1):
    def configPlot(ax, yTitle, yLim, xLim):
        ax.grid(which='major', axis='x', linestyle='-', linewidth=1, color='lightgray', alpha=0.5, zorder=-1)
        ax.grid(which='minor', axis='x', linestyle='--', linewidth=.75, color='lightgray', alpha=0.5, zorder=-1)

        ax.tick_params(axis='y', labelsize=10, length=4)
        ax.tick_params(
            axis='x', which='both', labelsize=9, pad=1, length=0,
            labeltop=True, top=True,
            labelbottom=False, bottom=False)

        ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(.75)
        ax.spines[['top', 'bottom', 'left', 'right']].set_color('#303030')
        ax.set_yticks([]), ax.set_ylim(xLim), ax.set_xlim(yLim), ax.set_xlabel(yTitle, labelpad=10, loc='left')

        ax.xaxis.tick_top(), ax.xaxis.set_label_position('top')
        ax.xaxis.set_major_locator(MultipleLocator(yLim[1] / 5))
        ax.xaxis.set_minor_locator(MultipleLocator(yLim[1] / 20))

    samples = [d['Sample'] for d in data_before]
    key = "k'" if 'Proportionality' in title else "n'"
    height_bef, height_bef_err = [d[f"{key}"] for d in data_before], [d[f"± {key}"] for d in data_before]
    height_aft, height_aft_err = [d[f"{key}"] for d in data_after], [d[f"± {key}"] for d in data_after]

    bin_width, space_samples, space_break = 1, 2.5, 2
    x = np.arange(space_samples * len(data_before))
    posList, labelsList = [], []

    configPlot(axes, title, (.0, lim), (x.min()-bin_width-.5, x.max()))

    for i in range(len(height_bef)):
        text_scale_correction = 1
        if scale_correction is not None:
            height_bef[i] = height_bef[i] / 10 if i == scale_correction else height_bef[i]
            height_bef_err[i] = height_bef_err[i] / 10 if i == scale_correction else height_bef_err[i]
            text_scale_correction = 10

        axes.barh(
            space_samples * x[i] - bin_width / space_break,
            width=height_bef[i], xerr=0,
            color=colors[i], edgecolor='#383838',
            height=bin_width, hatch=h, alpha=a, linewidth=.5,
            zorder=z)
        axes.errorbar(
            y=space_samples * x[i] - bin_width / space_break, x=height_bef[i],
            xerr=height_bef_err[i], color='#383838', alpha=.9,
            linewidth=1, capsize=5, capthick=1.05, zorder=3)
        axes.text(
            height_bef[i] + height_bef_err[i] + lim*.035,
            space_samples * x[i] - bin_width / space_break,
            f'{height_bef[i]*text_scale_correction:.{flt}f} ± {height_bef_err[i]*text_scale_correction:.{flt}f}',
            va='center_baseline', ha='left',
            color='#383838', fontsize=9)

    for i in range(len(height_aft)):
        text_scale_correction = 1
        if scale_correction is not None:
            height_aft[i] = height_aft[i] / 10 if i == scale_correction else height_aft[i]
            height_aft_err[i] = height_aft_err[i] / 10 if i == scale_correction else height_aft_err[i]
            text_scale_correction = 10

        axes.barh(
            space_samples * x[i] + bin_width / space_break,
            width=height_aft[i], xerr=0,
            left=0, color=colors[i], edgecolor='#383838',
            height=bin_width, hatch='////', alpha=a, linewidth=.5,
            zorder=2)
        axes.errorbar(
            y=space_samples * x[i] + bin_width / space_break,
            x=height_aft[i], xerr=height_aft_err[i],
            color='#383838', alpha=.99, linewidth=1, capsize=5, capthick=1.05,
            zorder=3)
        axes.text(
            height_aft[i] + height_aft_err[i] + lim*.034,
            space_samples * x[i] + bin_width / space_break,
            f'{height_aft[i]*text_scale_correction:.{flt}f} ± {height_aft_err[i]*text_scale_correction:.{flt}f}',
            va='center_baseline', ha='left',
            color='#383838', fontsize=9)

        if i == 5:
            posList.append(space_samples * x[i] - bin_width / space_break)
            posList.append(space_samples * x[i] + bin_width / space_break)
            labelsList.append('Before'), labelsList.append('After')
        if scale_correction is not None and i == scale_correction:
            posList.append(space_samples * x[i]), labelsList.append('10×')

    axes.set_yticks(posList)
    axes.set_yticklabels(labelsList)
    axes.invert_yaxis()


def midAxis(color, ax):
    ax[0].spines['right'].set_color(color)
    ax[1].spines['left'].set_color(color)
    ax[1].yaxis.tick_right()
    ax[1].yaxis.set_label_position('right')


def main(dataPath, fileName):
    fonts('C:/Users/petrus.kirsten/AppData/Local/Microsoft/Windows/Fonts/')
    plt.style.use('seaborn-v0_8-ticks')
    fig = plt.figure(figsize=(18, 7), facecolor='snow')
    gs = GridSpec(2, 3, width_ratios=[1.5, 1.5, 1.2], height_ratios=[1, 1])

    axPre, axPost = fig.add_subplot(gs[:, 0]), fig.add_subplot(gs[:, 1])
    axK, axN = fig.add_subplot(gs[0, 2]), fig.add_subplot(gs[1, 2])

    fig.suptitle(f'Viscoelastic recovery by frequency sweeps assay.')
    yTitle, yLimits = f"Storage modulus $G'$ (Pa)", (1 * 10 ** (-2), 1 * 10 ** 4)
    xTitle, xLimits = f'Frequency (Hz)', (0.09, 100)

    nSamples, colorSamples = getSamplesInfos(
        3, 4,
        2, 2, 2, 3,
        '#fb7e8f', '#e30057',
        '#80ed99', '#57cc99', '#38a3a5', '#22577a')

    data, labels = getSamplesData(dataPath, nSamples)

    listBefore, listAfter = {
        labels[0]: ([], [], []),
        labels[1]: ([], [], []),
        labels[2]: ([], [], []),
        labels[3]: ([], [], []),
        labels[4]: ([], [], []),
        labels[5]: ([], [], [])
    }, {
        labels[0]: ([], [], []),
        labels[1]: ([], [], []),
        labels[2]: ([], [], []),
        labels[3]: ([], [], []),
        labels[4]: ([], [], []),
        labels[5]: ([], [], [])
    }

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

    for k_a, k_b, c in zip(listAfter, listBefore, colorSamples):
        gP, gD = np.mean(listBefore[k_a][1], axis=1)[0], np.mean(listBefore[k_a][2], axis=1)[0]
        gPerr, gDerr = np.std(listBefore[k_a][1], axis=1)[0], np.std(listBefore[k_a][2], axis=1)[0]

        meanStorage, storageMeanErr, fitStart, fitEnd = getCteMean(gP)
        meanBefore.append(meanStorage)
        meanBeforeErr.append(storageMeanErr)

        dataFittingBef = plotFreqSweeps(  # Before axes
            sampleName=k_a,
            ax=axPre, x=np.mean(listBefore[k_a][0], axis=1)[0],
            yP=gP, yD=gD, yPerr=gPerr, yDerr=gDerr,
            axTitle='Before breakage', yLabel=yTitle, yLim=yLimits, xLabel=xTitle, xLim=xLimits, curveColor=c,
            logScale=True, startVal=fitStart, endVal=fitEnd, tableData=dataFittingBef)

        gP, gD = np.mean(listAfter[k_a][1], axis=1)[0], np.mean(listAfter[k_a][2], axis=1)[0]
        gPerr, gDerr = np.std(listAfter[k_a][1], axis=1)[0], np.std(listAfter[k_a][2], axis=1)[0]

        meanStorage, storageMeanErr, _, _ = getCteMean(gP)
        meanAfter.append(meanStorage)
        meanAfterErr.append(storageMeanErr)

        dataFittingAft = plotFreqSweeps(  # After axes
            sampleName=k_a,
            ax=axPost, x=np.mean(listAfter[k_a][0], axis=1)[0],
            yP=gP, yD=gD, yPerr=gPerr, yDerr=gDerr,
            axTitle='After breakage', yLabel=yTitle, yLim=yLimits, xLabel=xTitle, xLim=xLimits, curveColor=c,
            logScale=True, startVal=fitStart, endVal=fitEnd, tableData=dataFittingAft)

    # plotInsetMean(
    #     data=meanBefore, dataErr=meanBeforeErr,
    #     keys=listBefore.keys(), colors=colorSamples, ax=axPre)
    # plotInsetMean(
    #     data=meanAfter, dataErr=meanAfterErr,
    #     keys=listBefore.keys(), colors=colorSamples, ax=axPost,
    #     recovery=meanBefore)

    plotBars(
        "Proportionality coefficient $G_0'$ (Pa)", axK, 120,
        dataFittingBef, dataFittingAft, colorSamples, 1,
        scale_correction=1, z=1)

    plotBars(
        "Expoent index $n'$", axN, 2,
        dataFittingBef, dataFittingAft, colorSamples, 3,
        scale_correction=0, z=1)

    plt.subplots_adjust(
        wspace=0.164, hspace=0.24,
        top=0.91, bottom=0.1,
        left=0.045, right=0.98)
    plt.show()

    dirSave = Path(*Path(filePath[0]).parts[:Path(filePath[0]).parts.index('data') + 1])
    fig.savefig(f'{dirSave}' + f'\\{fileName}' + '.png', facecolor='w', dpi=600)
    print(f'\n\n· Chart saved at\n{dirSave}.')


if __name__ == '__main__':
    # folderPath = "C:/Users/petrus.kirsten/PycharmProjects/RheometerPlots/data"
    folderPath = "C:/Users/Petrus Kirsten/Documents/GitHub/RheometerPlots/data"

    filePath = [
        # kC
        folderPath + "/231024/kC/kC-viscoelasticRecovery-1.xlsx",
        folderPath + "/231024/kC/kC-viscoelasticRecovery-2.xlsx",
        folderPath + "/231024/kC/kC-viscoelasticRecovery-3.xlsx",

        # kC/CL
        folderPath + "/231024/kC_CL/kC_CL-viscoelasticRecovery-1.xlsx",
        folderPath + "/231024/kC_CL/kC_CL-viscoelasticRecovery-2.xlsx",
        folderPath + "/231024/kC_CL/kC_CL-viscoelasticRecovery-3.xlsx",
        folderPath + "/231024/kC_CL/kC_CL-viscoelasticRecovery-4.xlsx",

        # iC CL 14
        folderPath + "/311024/iC_CL_14/iC_CL_14-viscoelasticRecovery-1.xlsx",
        folderPath + "/311024/iC_CL_14/iC_CL_14-viscoelasticRecovery-2.xlsx",

        # iC CL 21
        folderPath + "/311024/iC_CL_21/iC_CL_21-viscoelasticRecovery-1.xlsx",
        folderPath + "/311024/iC_CL_21/iC_CL_21-viscoelasticRecovery-2.xlsx",

        # iC CL 28
        folderPath + "/311024/iC_CL_28/iC_CL_28-viscoelasticRecovery-1.xlsx",
        folderPath + "/311024/iC_CL_28/iC_CL_28-viscoelasticRecovery-2.xlsx",

        # iC CL 42
        folderPath + "/311024/iC_CL_42/iC_CL_42-viscoelasticRecovery-1.xlsx",
        folderPath + "/311024/iC_CL_42/iC_CL_42-viscoelasticRecovery-2.xlsx",
        folderPath + "/311024/iC_CL_42/iC_CL_42-viscoelasticRecovery-3.xlsx",
    ]

    main(filePath, 'Car-ViscoelasticRecovery')
