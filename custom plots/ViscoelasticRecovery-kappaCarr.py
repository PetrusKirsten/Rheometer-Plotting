import numpy as np
import pandas as pd
import seaborn as sns

from pathlib import Path
from math import ceil, sqrt
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
        n_kc_0,  n_kc_7, n_kc_14, n_kc_21, n_kc_28, n_kc_42,
        # colors
        color_kc_0, color_kc_7, color_kc_14, color_kc_21, color_kc_28, color_kc_42,
):
    number_samples = [n_kc_0,  n_kc_7, n_kc_14, n_kc_21, n_kc_28, n_kc_42]

    colors_samples = [color_kc_0, color_kc_7, color_kc_14, color_kc_21, color_kc_28, color_kc_42]

    return number_samples, colors_samples


def getSamplesData(
        dataPath,
        number_samples
):
    def getSegments(dataframe):
        freq = dataframe['f in Hz'].to_numpy()
        elastic = dataframe["G' in Pa"].to_numpy()
        loss = dataframe['G" in Pa'].to_numpy()
        delta = dataframe['tan(δ) in -'].to_numpy()

        seg2, seg3, seg5, seg6 = (
            dataframe.index[dataframe['SegIndex'] == seg].to_list()[0] for seg in ['2|1', '3|1', '5|1', '5|31'])

        # Slice segments
        segments = lambda arr: (arr[seg2:seg3], arr[seg5:seg6 + 1])  # Returns (constant segment, step segment)

        return {
            'freq': segments(freq),
            'storage': segments(elastic),
            'loss': segments(loss),
            'delta': segments(delta)
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
        dict_data[f'{sample_type}_delta'] = [s['delta'][0] for s in samples[sample_type]]

        dict_data[f'{sample_type}_freq_broken'] = [s['freq'][-1] for s in samples[sample_type]]
        dict_data[f'{sample_type}_storage_broken'] = [s['storage'][-1] for s in samples[sample_type]]
        dict_data[f'{sample_type}_loss_broken'] = [s['loss'][-1] for s in samples[sample_type]]
        dict_data[f'{sample_type}_delta_broken'] = [s['delta'][-1] for s in samples[sample_type]]

    return dict_data, sample_keys


def getRecoveryByFreq(storageList, data, frequencies):
    def round_to_nearest_125(value):
        return round(value / 0.125) * 0.125

    index_freqs = [0, 1, 4, 7, 10, 14, 17, 20, 24]
    index_freqs = range(0, 31, 1)
    freqList = []
    for index in index_freqs:
        storageList.append(data[index])
        freqList.append(round_to_nearest_125(frequencies[index]))

    freqList = [f"{freq:.3f} Hz" for freq in freqList]

    return storageList, freqList


def ratioElaVis(data_elastic, data_viscous):
    result = []

    for entry1 in data_elastic:
        for entry2 in data_viscous:
            if entry1['Sample'] == entry2['Sample']:
                if entry2["k'"] == 0:
                    new_entry = {
                        'Sample': entry1['Sample'],
                        "k'": None,
                        "± k'": None}
                else:
                    ratio_k_prime = entry1["k'"] / entry2["k'"]
                    uncertainty_k_prime = ratio_k_prime * sqrt(
                        (entry1["± k'"] / entry1["k'"]) ** 2 + (entry2["± k'"] / entry2["k'"]) ** 2)

                    rounded_ratio_k_prime = ceil(ratio_k_prime * 10) / 10
                    rounded_uncertainty_k_prime = ceil(uncertainty_k_prime * 10) / 10

                    new_entry = {
                        'Sample': entry1['Sample'],
                        "k'": rounded_ratio_k_prime,
                        "± k'": rounded_uncertainty_k_prime}
                result.append(new_entry)

    return result


def insertKey(keys):
    index = 1
    while index <= len(keys):
        keys.insert(index, 'Broken')
        index += 2
    return keys


def plotFreqSweeps(sampleName, axTop, axBottom, axTitle,
                   x, yP, yD, yPerr, yDerr,
                   yLabel, yLim, xLabel, xLim,
                   curveColor,
                   logScale=True,
                   startVal=0, endVal=16,
                   tableDataStor=None, tableDataLoss=None):
    def legendLabel(ax):
        legend = ax.legend(fancybox=False, frameon=True, framealpha=0.9, fontsize=9)
        legend.get_frame().set_facecolor('w')
        legend.get_frame().set_edgecolor('whitesmoke')

    def configPlot(ax, axisColor='#303030'):
        ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(0.75)
        ax.spines[['top', 'bottom', 'left', 'right']].set_color(axisColor)
        ax.tick_params(axis='both', which='both', colors=axisColor)

        ax.grid(True, which='major', axis='y', linestyle='-', linewidth=.75, color='lightgray', alpha=0.5)
        ax.grid(True, which='minor', axis='y', linestyle='-', linewidth=.5, color='lightgray', alpha=0.5)

        ax.set_xscale('log' if logScale else 'linear')
        ax.set_xlim(xLim)

        ax.set_ylabel(f'{yLabel}', color=axisColor)
        ax.set_yscale('log' if logScale else 'linear')
        ax.set_ylim(yLim)

    yTitleBottom = f"Viscous modulus $G''$ (Pa)"
    configPlot(axTop), configPlot(axBottom)
    axTop.set_title(axTitle, size=10, color='k'), axTop.set_xticklabels([])
    axBottom.set_xlabel(f'{xLabel}', color='#303030'), axBottom.set_ylabel(f'{yTitleBottom}', color='#303030')

    x_toFit_stor, y_toFit_stor = arraySplit(x, yP, startVal, endVal)
    params_stor, covariance_stor = curve_fit(powerLaw, x_toFit_stor, y_toFit_stor)  # p0=(y_mean[0], y_mean[-1], 100))
    errors_stor = np.sqrt(np.diag(covariance_stor))
    tableDataStor = exportFit(
        f'{sampleName}',
        params_stor, errors_stor,
        tableDataStor)

    axTop.errorbar(
        x[:-1], yP[:-1], yPerr[:-1],
        color=curveColor, alpha=.85,
        fmt='none', mfc=curveColor,
        capsize=2.5, capthick=1, linestyle='', lw=1,
        label=f'', zorder=2)
    axTop.errorbar(
        x[:-1], yP[:-1], 0,
        color=curveColor, alpha=.65,
        fmt='o', markersize=5.4,
        mfc=curveColor, mec='#383838', mew=.75,
        linestyle='',
        label=f'{sampleName}', zorder=3)

    x_toFit_loss, y_toFit_loss = arraySplit(x, yD, startVal, endVal)
    params_loss, covariance_loss = curve_fit(powerLaw, x_toFit_loss, y_toFit_loss)  # p0=(y_mean[0], y_mean[-1], 100))
    errors_loss = np.sqrt(np.diag(covariance_loss))
    tableDataLoss = exportFit(
        f'{sampleName}',
        params_loss, errors_loss,
        tableDataLoss)

    axBottom.errorbar(
        x[:-1], yD[:-1], yDerr[:-1],
        color=curveColor, alpha=.85,
        fmt='none', mfc=curveColor,
        capsize=2.5, capthick=1, linestyle='', lw=1,
        label=f'', zorder=2)
    axBottom.errorbar(
        x[:-1], yD[:-1], 0,
        color=curveColor, alpha=.75,
        fmt='o', markersize=5.4,
        mfc='w', mec=curveColor, mew=1,
        linestyle='',
        label=f'{sampleName}', zorder=3)

    if axTitle == 'After breakage':
        legendLabel(axTop)

    return tableDataStor, tableDataLoss


def plotBars(
        title, axes, lim,
        data_before, data_after,
        colors, dec, scale_correction=None,
        a=.9, h='', z=1):
    def configPlot(ax, yTitle, yLim, xLim):
        ax.tick_params(axis='y', labelsize=10, length=4)
        ax.tick_params(
            axis='x', which='both', labelsize=9, pad=1, length=0,
            labeltop=False, top=False,
            labelbottom=False, bottom=False)

        ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(.75)
        ax.spines[['top', 'bottom', 'left', 'right']].set_color('#303030')
        ax.set_yticks([]), ax.set_ylim(xLim), ax.set_xlim(yLim)
        ax.set_xlabel(yTitle, size=10, labelpad=5, loc='center')

        ax.xaxis.tick_top(), ax.xaxis.set_label_position('top')
        ax.xaxis.set_major_locator(MultipleLocator(yLim[1] / 5))
        ax.xaxis.set_minor_locator(MultipleLocator(yLim[1] / 20))

    samples = [d['Sample'] for d in data_before]
    key = "n'" if "n'" in title else "k'"
    height_bef, height_bef_err = [d[f"{key}"] for d in data_before], [d[f"± {key}"] for d in data_before]
    height_aft, height_aft_err = [d[f"{key}"] for d in data_after], [d[f"± {key}"] for d in data_after]

    bin_width, space_samples, space_break = 1, 2.5, 2
    x = np.arange(space_samples * len(data_before))
    posList, labelsList = [], []

    configPlot(axes, title, (.0, lim), (x.min() - bin_width - .5, x.max()))

    for i in range(len(height_bef)):
        text_scale_correction = 1
        if scale_correction is not None:
            height_bef[i] = height_bef[i] / 10 if i == scale_correction else height_bef[i]
            height_bef_err[i] = height_bef_err[i] / 10 if i == scale_correction else height_bef_err[i]
            text_scale_correction = 10 if i == scale_correction else 1

        axes.barh(
            space_samples * x[i] - bin_width / space_break,
            width=height_bef[i] if height_bef_err[i] < height_bef[i] else 0, xerr=0,
            color=colors[i], edgecolor='#383838', alpha=a,
            height=bin_width, hatch=h, linewidth=.5,
            zorder=z)
        axes.errorbar(
            y=space_samples * x[i] - bin_width / space_break,
            x=height_bef[i] if height_bef_err[i] < height_bef[i] else 0,
            xerr=height_bef_err[i] if height_bef_err[i] < height_bef[i] else 0,
            color='#383838', alpha=.9,
            linewidth=1, capsize=3, capthick=1.05, zorder=3)

        text = (f'{ceil(height_bef[i] * text_scale_correction * 100) / 100:.{dec}f} '
                f'± {ceil(height_bef_err[i] * text_scale_correction * 100) / 100:.{dec}f}')
        axes.text(
            lim * .975,
            space_samples * x[i] + 0.1 - bin_width / space_break,
            text if height_bef_err[i] < height_bef[i] else 'Not fitted',
            va='center_baseline', ha='left',
            color='#383838', fontsize=9)

    for i in range(len(height_aft)):
        text_scale_correction = 1
        if scale_correction is not None:
            height_aft[i] = height_aft[i] / 10 if i == scale_correction else height_aft[i]
            height_aft_err[i] = height_aft_err[i] / 10 if i == scale_correction else height_aft_err[i]
            text_scale_correction = 10 if i == scale_correction else 1

        axes.barh(
            space_samples * x[i] + bin_width / space_break,
            width=height_aft[i] if height_aft_err[i] < height_aft[i] else 0, xerr=0,
            color=colors[i], edgecolor='#383838', alpha=a,
            height=bin_width, hatch='////', linewidth=.5,
            zorder=2)
        axes.errorbar(
            y=space_samples * x[i] + bin_width / space_break,
            x=height_aft[i] if height_aft_err[i] < height_aft[i] else 0,
            xerr=height_aft_err[i] if height_aft_err[i] < height_aft[i] else 0,
            color='#383838', alpha=.99, linewidth=1, capsize=3, capthick=1.05,
            zorder=3)

        text = (f'{ceil(height_aft[i] * text_scale_correction * 100) / 100:.{dec}f} '
                f'± {ceil(height_aft_err[i] * text_scale_correction * 100) / 100:.{dec}f}')
        axes.text(
            lim * .975,
            space_samples * x[i] + 0.1 + bin_width / space_break,
            text if height_aft_err[i] < height_aft[i] else 'Not fitted',
            va='center_baseline', ha='left',
            color='#383838', fontsize=9)

        if i == 5:
            posList.append(space_samples * x[i] - bin_width / space_break)
            posList.append(space_samples * x[i] + bin_width / space_break)
            labelsList.append('Before'), labelsList.append('After')
        if scale_correction is not None and i == scale_correction:
            posList.append(space_samples * x[i]), labelsList.append('10×')

    axes.yaxis.tick_right()
    axes.set_yticks(posList)
    axes.set_yticklabels(labelsList)
    axes.invert_yaxis(), axes.invert_xaxis()


def plotHeatMap(
        title,
        data_map, frequencies, formulations,
):
    fonts('C:/Users/petrus.kirsten/AppData/Local/Microsoft/Windows/Fonts/')

    plt.rcParams['text.color'] = '#383838'  # Cor de todos os textos
    plt.rcParams['axes.labelcolor'] = '#383838'  # Cor dos labels dos eixos
    plt.rcParams['xtick.color'] = '#383838'  # Cor das ticks do eixo X
    plt.rcParams['ytick.color'] = '#383838'  # Cor das ticks do eixo Y
    plt.rcParams['axes.titlecolor'] = '#383838'  # Cor do título do gráfico

    plt.style.use('seaborn-v0_8-ticks')
    plt.figure(figsize=(18, 4.5), facecolor='snow')
    plt.gca().set_facecolor('snow')

    if title == "Elastic recovery (%)":
        for i in range(len(data_map)):
            for j in range(len(data_map[i])):
                if data_map[i][j] > 100:
                    data_map[i][j] = None
        colors = 'RdYlGn'
        decimal = '.0f'

    if title == 'Loss factor $tan(\delta)$':
        for i in range(len(data_map)):
            for j in range(len(data_map[i])):
                for k in range(len(data_map[i][j])):
                    if data_map[i][j][k] > 2.1:
                        data_map[i][j][k] = None
        colors = 'coolwarm'
        decimal = '.2f'
        data_map = np.array(data_map, dtype=float).flatten().reshape(12, 31)

    df = pd.DataFrame(data_map, index=formulations, columns=frequencies)
    sns.heatmap(
        df,
        annot=True, cmap=colors,
        fmt=decimal, linewidths=0.5, cbar_kws={'label': title})

    # plt.title("Elastic modulus $G'$ recovery across frequency.")
    plt.tick_params(axis='both', which='both', length=0)
    plt.xticks(rotation=45, ha='center', fontsize=10, color='#383838')
    plt.yticks(rotation=0, ha='right', fontsize=10, color='#383838')

    plt.subplots_adjust(
        wspace=0, hspace=0,
        top=0.97, bottom=0.14,
        left=0.05, right=1.0)

    dirSave = Path(*Path(filePath[0]).parts[:Path(filePath[0]).parts.index('data') + 1])
    plt.savefig(f'{dirSave}' + f'\\{title[0]}' + '.png', facecolor='w', dpi=600)


def main(dataPath, fileName):
    fonts('C:/Users/petrus.kirsten/AppData/Local/Microsoft/Windows/Fonts/')
    plt.style.use('seaborn-v0_8-ticks')
    fig = plt.figure(figsize=(18, 10), facecolor='snow')
    gs = GridSpec(4, 3, width_ratios=[1.5, 1.5, 1.2], height_ratios=[1, 1, 1, 1])

    axPreTop, axPostTop = fig.add_subplot(gs[:2, 0]), fig.add_subplot(gs[:2, 1])
    axPreBottom, axPostBottom = fig.add_subplot(gs[2:, 0]), fig.add_subplot(gs[2:, 1])
    axBar1, axBar2, axBar3, axBar4 = (fig.add_subplot(gs[0, 2]),
                                      fig.add_subplot(gs[1, 2]),
                                      fig.add_subplot(gs[2, 2]),
                                      fig.add_subplot(gs[3, 2]))

    fig.suptitle(f'Viscoelastic recovery by frequency sweeps assay.')
    yTitle, yLimits = f"Elastic modulus $G'$ (Pa)", (1 * 10 ** (-2), 1 * 10 ** 4)
    xTitle, xLimits = f'Frequency (Hz)', (.075, 100)

    nSamples, colorSamples = getSamplesInfos(
        3, 4,
        2, 2, 2, 3,
        '#fb7e8f', '#e30057',
        '#80ed99', '#57cc99', '#38a3a5', '#22577a')
    data, labels = getSamplesData(dataPath, nSamples)

    listBefore = {
        labels[0]: ([], [], [], []),
        labels[1]: ([], [], [], []),
        labels[2]: ([], [], [], []),
        labels[3]: ([], [], [], []),
        labels[4]: ([], [], [], []),
        labels[5]: ([], [], [], []),
    }
    listAfter = {
        labels[0]: ([], [], [], []),
        labels[1]: ([], [], [], []),
        labels[2]: ([], [], [], []),
        labels[3]: ([], [], [], []),
        labels[4]: ([], [], [], []),
        labels[5]: ([], [], [], []),
    }

    meanBefore, meanAfter = [], []
    meanBeforeErr, meanAfterErr = [], []
    dataFittingBef_stor, dataFittingAft_stor = [], []
    dataFittingBef_loss, dataFittingAft_loss = [], []
    recoveryPCT, tan_delta = [], []

    for key, (x, gP, gD, d) in listBefore.items():
        x.append(data[f'{key}_freq'])
        gP.append(data[f'{key}_storage'])
        gD.append(data[f'{key}_loss'])
        d.append(data[f'{key}_delta'])

    for key, (x, gP, gD, d) in listAfter.items():
        x.append(data[f'{key}_freq_broken'])
        gP.append(data[f'{key}_storage_broken'])
        gD.append(data[f'{key}_loss_broken'])
        d.append(data[f'{key}_delta_broken'])

    for key, color in zip(listBefore, colorSamples):
        recoveryBef, recoveryAft, delta_bef, delta_aft = [], [], [], []

        freqs = np.mean(listBefore[key][0], axis=1)[0]
        gP, gD = np.mean(listBefore[key][1], axis=1)[0], np.mean(listBefore[key][2], axis=1)[0]
        gPerr, gDerr = np.std(listBefore[key][1], axis=1)[0], np.std(listBefore[key][2], axis=1)[0]
        delta, deltaErr = np.mean(listBefore[key][3], axis=1)[0], np.std(listBefore[key][3], axis=1)[0]

        recoveryBef, freqsRecovery = getRecoveryByFreq(recoveryBef, gP, freqs)
        delta_bef, freqsRecovery = getRecoveryByFreq(delta_bef, delta, freqs)

        meanStorage, storageMeanErr, fitStart, fitEnd = getCteMean(gP)
        meanBefore.append(meanStorage)
        meanBeforeErr.append(storageMeanErr)

        dataFittingBef_stor, dataFittingBef_loss = plotFreqSweeps(  # Before axes
            sampleName=key,
            axTop=axPreTop, axBottom=axPreBottom,
            x=freqs, yP=gP, yD=gD, yPerr=gPerr, yDerr=gDerr,
            axTitle='Before breakage',
            yLabel=yTitle, yLim=yLimits, xLabel=xTitle, xLim=xLimits, curveColor=color,
            logScale=True,
            tableDataStor=dataFittingBef_stor, tableDataLoss=dataFittingBef_loss)

        gP = listAfter[key][1] if key != '0St + kCar/CL' else listAfter[key][1][0][0]
        gD = listAfter[key][2] if key != '0St + kCar/CL' else listAfter[key][2][0][0]

        gPerr = np.std(gP, axis=1)[0] if key != '0St + kCar/CL' else np.zeros(31)
        gDerr = np.std(gD, axis=1)[0] if key != '0St + kCar/CL' else np.zeros(31)
        gP = np.mean(gP, axis=1)[0] if key != '0St + kCar/CL' else gP
        gD = np.mean(gD, axis=1)[0] if key != '0St + kCar/CL' else gD
        delta, deltaErr = np.mean(listAfter[key][3], axis=1)[0], np.std(listAfter[key][3], axis=1)[0]

        recoveryAft, freqsRecovery = getRecoveryByFreq(recoveryAft, gP, freqs)
        delta_aft, freqsRecovery = getRecoveryByFreq(delta_aft, delta, freqs)
        tan_delta.append([delta_bef, delta_aft])

        dataFittingAft_stor, dataFittingAft_loss = plotFreqSweeps(  # After axes
            sampleName=key,
            axTop=axPostTop, axBottom=axPostBottom,
            x=freqs, yP=gP, yD=gD, yPerr=gPerr, yDerr=gDerr,
            axTitle='After breakage',
            yLabel=yTitle, yLim=yLimits, xLabel=xTitle, xLim=xLimits, curveColor=color,
            logScale=True,
            tableDataStor=dataFittingAft_stor, tableDataLoss=dataFittingAft_loss)
        axPostTop.set_ylabel(''), axPostBottom.set_ylabel('')
        axPostTop.set_yticklabels([]), axPostBottom.set_yticklabels([])

        recoveryPCT.append(((np.array(recoveryAft) / np.array(recoveryBef)) * 100).tolist())
        # recoveryPCT = recoveryPCT[0]

    ratioBef = ratioElaVis(dataFittingBef_stor, dataFittingBef_loss)
    ratioAft = ratioElaVis(dataFittingAft_stor, dataFittingAft_loss)

    plotBars(  # First table
        "$n'$", axBar1, 1.1,
        dataFittingBef_stor, dataFittingAft_stor, colorSamples, dec=2,
        scale_correction=0, z=1)

    plotBars(  # Second table
        "$G_0'$ (Pa)", axBar2, 110,
        dataFittingBef_stor, dataFittingAft_stor, colorSamples, dec=1,
        scale_correction=1, z=1)

    plotBars(  # Third table
        "$G_0''$ (Pa)", axBar3, 12,
        dataFittingBef_loss, dataFittingAft_loss, colorSamples, dec=1,
        scale_correction=1, z=1)

    plotBars(  # Fourth table
        "$G_0'\,/\,G_0''$", axBar4, 20,
        ratioBef, ratioAft, colorSamples, dec=1,
        scale_correction=None, z=1)

    plt.subplots_adjust(
        wspace=0.015, hspace=0.15,
        top=0.93, bottom=0.07,
        left=0.045, right=0.965)

    plotHeatMap(
        "Elastic recovery (%)",
        recoveryPCT, freqsRecovery, labels)

    brokenLabels = insertKey(labels)
    plotHeatMap(
        'Loss factor $tan(\delta)$',
        tan_delta, freqsRecovery, brokenLabels)
    plt.show()

    dirSave = Path(*Path(filePath[0]).parts[:Path(filePath[0]).parts.index('data') + 1])
    fig.savefig(f'{dirSave}' + f'\\{fileName}' + '.png', facecolor='w', dpi=600)
    print(f'\n\n· Charts saved at\n{dirSave}.')


if __name__ == '__main__':
    folderPath = "C:/Users/petrus.kirsten/PycharmProjects/RheometerPlots/data"  # CEBB
    # folderPath = "C:/Users/Petrus Kirsten/Documents/GitHub/RheometerPlots/data"  # Personal

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

    main(filePath, 'Car-ViscoelasticRecoveryWithViscous')
