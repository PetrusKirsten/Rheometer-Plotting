import numpy as np
import pandas as pd
import seaborn as sns

from pathlib import Path
from math import ceil, sqrt

from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt, ticker
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator
from matplotlib.font_manager import FontProperties


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


def arraySplit(xArr, yArr, startValue, endValue):
    startIndex, endIndex = np.where(xArr >= startValue)[0][0], np.where(xArr <= endValue)[0][-1]

    return xArr[startIndex:endIndex], yArr[startIndex:endIndex]


def exportFit(
        sample,
        data, err,
        table, mode=''
):
    if mode == 'thixo':
        keys = ('$\\tau_0$', '$\\tau_e$', '$\lambda$')
    elif mode == 'HB':
        keys = ("k'", "n'", 'sigma_zero')
    else:
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


def downsampler(array, n=200):
    if len(array) > n:
        step = len(array) // n  # Calculate step size
        return array[::step][:n]
    return array


class Recovery:
    def __init__(
            self,
            dataPath, fileName,
            names_samples, number_samples, colors_samples
    ):

        def getData():
            def getSegments(dataframe):
                freq = dataframe['f in Hz'].to_numpy()
                elastic = dataframe["G' in Pa"].to_numpy()
                loss = dataframe['G" in Pa'].to_numpy()
                delta = dataframe['tan(δ) in -'].to_numpy()

                if 'off' in path:
                    seg2, seg3 = (
                        dataframe.index[dataframe['SegIndex'] == seg].to_list()[0] for seg in ['2|1', '3|1'])
                    segments = lambda arr: (arr[seg2:seg3], arr[seg2:seg3])

                else:
                    seg2, seg3, seg5, seg6 = (
                        dataframe.index[dataframe['SegIndex'] == seg].to_list()[0] for seg in
                        ['2|1', '3|1', '5|1', '5|31'])
                    segments = lambda arr: (arr[seg2:seg3], arr[seg5:seg6 + 1])

                return {
                    'freq': segments(freq),
                    'storage': segments(elastic),
                    'loss': segments(loss),
                    'delta': segments(delta)
                }

            def dict_OscFreqSweeps(labels):
                return {label: ([], [], [], []) for label in labels}

            self.sample_keys = list(self.names_samples.keys())
            if len(self.sample_keys) != len(self.number_samples):
                raise ValueError('The length of "number_samples" must match the number of sample keys.')
            sample_labels = [
                key for key, count in zip(self.sample_keys, self.number_samples) for _ in range(count)
            ]

            for sample_type, path in zip(sample_labels, self.dataPath):
                df = pd.read_excel(path)
                segments = getSegments(df)
                self.names_samples[sample_type].append(segments)

            dict_data = {}
            for sample_type in self.names_samples:
                dict_data[f'{sample_type}_freq'] = [s['freq'][0] for s in self.names_samples[sample_type]]
                dict_data[f'{sample_type}_storage'] = [s['storage'][0] for s in self.names_samples[sample_type]]
                dict_data[f'{sample_type}_loss'] = [s['loss'][0] for s in self.names_samples[sample_type]]
                dict_data[f'{sample_type}_delta'] = [s['delta'][0] for s in self.names_samples[sample_type]]

                dict_data[f'{sample_type}_freq_broken'] = [s['freq'][-1] for s in self.names_samples[sample_type]]
                dict_data[f'{sample_type}_storage_broken'] = [s['storage'][-1] for s in self.names_samples[sample_type]]
                dict_data[f'{sample_type}_loss_broken'] = [s['loss'][-1] for s in self.names_samples[sample_type]]
                dict_data[f'{sample_type}_delta_broken'] = [s['delta'][-1] for s in self.names_samples[sample_type]]

            return dict_data, dict_OscFreqSweeps(self.sample_keys), dict_OscFreqSweeps(self.sample_keys)

        def appendData(
                inputList,
                isBroken=False
        ):
            broken = ''
            if isBroken:
                broken = '_broken'
            for key, (x, gP, gD, d) in inputList.items():
                x.append(self.data[f'{key}_freq' + f'{broken}'])
                gP.append(self.data[f'{key}_storage' + f'{broken}'])
                gD.append(self.data[f'{key}_loss' + f'{broken}'])
                d.append(self.data[f'{key}_delta' + f'{broken}'])

            return inputList

        # input vars
        self.dataPath = dataPath
        self.fileName = fileName
        self.names_samples = names_samples
        self.number_samples = number_samples
        self.colors_samples = colors_samples

        # data vars
        self.meanBefore, self.meanAfter = [], []
        self.meanBeforeErr, self.meanAfterErr = [], []
        self.dataFittingBef_stor, self.dataFittingAft_stor = [], []
        self.dataFittingBef_loss, self.dataFittingAft_loss = [], []
        self.recoveryPCT, self.freqsRecovery, self.tan_delta = [], [], []
        self.ratioBef, self.ratioAft, self.ratioReady = None, None, False
        # data reading
        self.data, self.listBefore, self.listAfter = getData()
        self.listBefore = appendData(self.listBefore)
        self.listAfter = appendData(self.listAfter, True)

        # chart config
        fonts('C:/Users/petrus.kirsten/AppData/Local/Microsoft/Windows/Fonts/')
        plt.style.use('seaborn-v0_8-ticks')

    def plotGmodulus(
            self,
            yTitle, yLimits,
            xTitle, xLimits,
            show=True, save=False
    ):
        def getValues(inputList, loopKey, condition=False):
            if not condition:
                frequencies = np.mean(inputList[loopKey][0], axis=1)[0]

                elasticMod, viscousMod = (np.mean(inputList[loopKey][1], axis=1)[0],
                                          np.mean(inputList[loopKey][2], axis=1)[0])

                elasticModErr, viscousModErr = (np.std(inputList[loopKey][1], axis=1)[0],
                                                np.std(inputList[loopKey][2], axis=1)[0])

                lossFactor, lossFactorErr = (np.mean(inputList[loopKey][3], axis=1)[0],
                                             np.std(inputList[loopKey][3], axis=1)[0])

            else:
                frequencies = np.mean(inputList[loopKey][0], axis=1)[0]

                elasticMod, viscousMod = (inputList[key][1] if key != condition else inputList[key][1][0][0],
                                          inputList[key][2] if key != condition else inputList[key][2][0][0])

                elasticModErr, viscousModErr = (np.std(elasticMod, axis=1)[0] if key != condition else np.zeros(31),
                                                np.std(viscousMod, axis=1)[0] if key != condition else np.zeros(31))

                elasticMod, viscousMod = (np.mean(elasticMod, axis=1)[0] if key != condition else elasticMod,
                                          np.mean(viscousMod, axis=1)[0] if key != condition else viscousMod)

                lossFactor, lossFactorErr = (np.mean(inputList[loopKey][3], axis=1)[0],
                                             np.std(inputList[loopKey][3], axis=1)[0])

            return frequencies, elasticMod, viscousMod, elasticModErr, viscousModErr, lossFactor, lossFactorErr

        def getValuesByFreq(storageList, data, frequencies, index_freqs=range(3, 23, 1)):
            def round_to_nearest(value):
                return round(value / 0.05) * 0.05

            freqList = []
            for index in index_freqs:
                storageList.append(data[index])
                freqList.append(round_to_nearest(frequencies[index]))

            freqList = [f"{freq:.2f} Hz" for freq in freqList]

            return storageList, freqList

        def cteRegionMean(values, tolerance=100):
            """
            :param values: to be analysed
            :param tolerance: the difference betweem two points data
            :return: the mean of the "cte" region and its indexes
            """
            diffs = np.abs(np.diff(values))
            constantRegions = diffs < tolerance

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

        def powerLaw(omega, kPrime, nPrime):
            return kPrime * (omega ** nPrime)

        def viscPerElas(data_elastic, data_viscous):
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

        def drawPlot(
                sampleName, axTop, axBottom, axTitle,
                x, yP, yD, yPerr, yDerr,
                yLabel, yLim, xLabel, xLim,
                curveColor,
                logScale=True,
                startVal=0, endVal=16,
                tableDataStor=None, tableDataLoss=None
        ):
            def legendLabel(ax):
                legend = ax.legend(loc='lower right', fancybox=False, frameon=True, framealpha=0.9, fontsize=9)
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
            params_stor, covariance_stor = curve_fit(powerLaw, x_toFit_stor, y_toFit_stor)
            # p0=(y_mean[0], y_mean[-1], 100))
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
            params_loss, covariance_loss = curve_fit(powerLaw, x_toFit_loss, y_toFit_loss)
            # p0=(y_mean[0], y_mean[-1], 100))
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

        fig = plt.figure(figsize=(13, 9), facecolor='snow')
        fig.canvas.manager.set_window_title(self.fileName + ' - Elastic and viscous moduli')
        gs = GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])

        fig.suptitle(f'')
        axPreTop, axPostTop = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
        axPreBottom, axPostBottom = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])

        for key, color in zip(self.listBefore, self.colors_samples):
            recoveryBef, recoveryAft, delta_bef, delta_aft = [], [], [], []

            freqs, gP, gD, gPerr, gDerr, delta, deltaErr = getValues(self.listBefore, key)

            recoveryBef, _ = getValuesByFreq(recoveryBef, gP, freqs)
            delta_bef, _ = getValuesByFreq(delta_bef, delta, freqs)

            meanStorage, meanStorageErr, fitStart, fitEnd = cteRegionMean(gP)
            self.meanBefore.append(meanStorage), self.meanBeforeErr.append(meanStorageErr)

            self.dataFittingBef_stor, self.dataFittingBef_loss = drawPlot(  # Before axes
                sampleName=key,
                axTop=axPreTop, axBottom=axPreBottom,
                x=freqs, yP=gP, yD=gD, yPerr=gPerr, yDerr=gDerr,
                axTitle='Before breakage',
                yLabel=yTitle, yLim=yLimits, xLabel=xTitle, xLim=xLimits, curveColor=color,
                logScale=True,
                tableDataStor=self.dataFittingBef_stor, tableDataLoss=self.dataFittingBef_loss
            )

            freqs, gP, gD, gPerr, gDerr, delta, deltaErr = getValues(self.listAfter, key, 'St + kCar/CL_7')

            recoveryAft, _ = getValuesByFreq(recoveryAft, gP, freqs)
            delta_aft, self.freqsRecovery = getValuesByFreq(delta_aft, delta, freqs)
            self.tan_delta.append([delta_bef, delta_aft])

            meanStorage, meanStorageErr, fitStart, fitEnd = cteRegionMean(gP)
            self.meanAfter.append(meanStorage), self.meanAfterErr.append(meanStorageErr)

            self.dataFittingAft_stor, self.dataFittingAft_loss = drawPlot(  # After axes
                sampleName=key,
                axTop=axPostTop, axBottom=axPostBottom,
                x=freqs, yP=gP, yD=gD, yPerr=gPerr, yDerr=gDerr,
                axTitle='After breakage',
                yLabel=yTitle, yLim=yLimits,
                xLabel=xTitle, xLim=xLimits,
                curveColor=color, logScale=True,
                tableDataStor=self.dataFittingAft_stor, tableDataLoss=self.dataFittingAft_loss
            )

            axPostTop.set_ylabel(''), axPostBottom.set_ylabel('')
            axPostTop.set_yticklabels([]), axPostBottom.set_yticklabels([])

            self.recoveryPCT.append(((np.array(recoveryAft) / np.array(recoveryBef)) * 100).tolist())

        self.ratioBef = viscPerElas(self.dataFittingBef_stor, self.dataFittingBef_loss)
        self.ratioAft = viscPerElas(self.dataFittingAft_stor, self.dataFittingAft_loss)
        self.ratioReady = True

        plt.subplots_adjust(
            wspace=0.015, hspace=0.060,
            top=0.970, bottom=0.070,
            left=0.060, right=0.985)
        if show:
            plt.show()
        if save:
            dirSave = Path(*Path(self.dataPath[0]).parts[:Path(self.dataPath[0]).parts.index('data') + 1])
            fig.savefig(
                f'{dirSave}' + f'\\{self.fileName}' + ' - Elastic and viscous moduli' + '.png',
                facecolor='w', dpi=600)
            print(f'\n\n· Elastic and viscous moduli chart saved at:\n{dirSave}.')

    def plotBars(
            self,
            limits,
            corrections,
            show=True, save=False
    ):

        def drawBars(
                title, axes, lim,
                data_before, data_after,
                colors, dec, scale_correction=None,
                a=.9, h='', z=1
        ):

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

            # data from before
            for i in range(len(height_bef)):
                text_scale_correction = 1
                if scale_correction is not None:
                    height_bef[i] = height_bef[i] / 10 if i == scale_correction else height_bef[i]
                    height_bef_err[i] = height_bef_err[i] / 10 if i == scale_correction else height_bef_err[i]
                    text_scale_correction = 10 if i == scale_correction else 1

                axes.barh(
                    space_samples * x[i] - bin_width / space_break,
                    width=height_bef[i] if height_bef_err[i] < height_bef[i] / 2 else 0, xerr=0,
                    color=colors[i], edgecolor='#383838', alpha=a,
                    height=bin_width, hatch=h, linewidth=.5,
                    zorder=z)
                axes.errorbar(
                    y=space_samples * x[i] - bin_width / space_break,
                    x=height_bef[i] if height_bef_err[i] < height_bef[i] / 2 else 0,
                    xerr=height_bef_err[i] if height_bef_err[i] < height_bef[i] / 2 else 0,
                    color='#383838', alpha=.9,
                    linewidth=1, capsize=3, capthick=1.05, zorder=3)

                text = (f'{ceil(height_bef[i] * text_scale_correction * 100) / 100:.{dec}f} '
                        f'± {ceil(height_bef_err[i] * text_scale_correction * 100) / 100:.{dec}f}')
                axes.text(
                    lim * .975,
                    space_samples * x[i] + 0.1 - bin_width / space_break,
                    text if height_bef_err[i] < height_bef[i] / 2 else 'Not fitted',
                    va='center_baseline', ha='left',
                    color='#383838', fontsize=9)
            # data from after
            for i in range(len(height_aft)):
                text_scale_correction = 1
                if scale_correction is not None:
                    height_aft[i] = height_aft[i] / 10 if i == scale_correction else height_aft[i]
                    height_aft_err[i] = height_aft_err[i] / 10 if i == scale_correction else height_aft_err[i]
                    text_scale_correction = 10 if i == scale_correction else 1

                axes.barh(
                    space_samples * x[i] + bin_width / space_break,
                    width=height_aft[i] if height_aft_err[i] < height_aft[i] / 2 else 0, xerr=0,
                    color=colors[i], edgecolor='#383838', alpha=a,
                    height=bin_width, hatch='////', linewidth=.5,
                    zorder=2)
                axes.errorbar(
                    y=space_samples * x[i] + bin_width / space_break,
                    x=height_aft[i] if height_aft_err[i] < height_aft[i] / 2 else 0,
                    xerr=height_aft_err[i] if height_aft_err[i] < height_aft[i] / 2 else 0,
                    color='#383838', alpha=.99, linewidth=1, capsize=3, capthick=1.05,
                    zorder=3)

                text = (f'{ceil(height_aft[i] * text_scale_correction * 100) / 100:.{dec}f} '
                        f'± {ceil(height_aft_err[i] * text_scale_correction * 100) / 100:.{dec}f}')
                axes.text(
                    lim * .975,
                    space_samples * x[i] + 0.1 + bin_width / space_break,
                    text if height_aft_err[i] < height_aft[i] / 2 else 'Not fitted',
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
            # axes.invert_yaxis()
            axes.invert_xaxis()

        # figure configs
        fig = plt.figure(figsize=(5, 9), facecolor='snow')
        fig.canvas.manager.set_window_title(self.fileName + ' - Bars plots')
        gs = GridSpec(4, 1, height_ratios=[1, 1, 1, 1], width_ratios=[1])

        axBar1, axBar2, axBar3, axBar4 = (
            fig.add_subplot(gs[0, 0]),
            fig.add_subplot(gs[1, 0]),
            fig.add_subplot(gs[2, 0]),
            fig.add_subplot(gs[3, 0]))

        fig.suptitle(f'')

        # bars data
        drawBars(  # First table
            "$n'$", axBar1, limits[0],
            self.dataFittingBef_stor, self.dataFittingAft_stor, self.colors_samples, dec=2,
            scale_correction=corrections[0], z=1)

        drawBars(  # Second table
            "$G_0'$ (Pa)", axBar2, limits[1],
            self.dataFittingBef_stor, self.dataFittingAft_stor, self.colors_samples, dec=1,
            scale_correction=corrections[1], z=1)

        drawBars(  # Third table
            "$G_0''$ (Pa)", axBar3, limits[2],
            self.dataFittingBef_loss, self.dataFittingAft_loss, self.colors_samples, dec=1,
            scale_correction=corrections[2], z=1)

        drawBars(  # Fourth table
            "Loss factor $G_0''\,/\,G_0'$", axBar4, limits[3],
            self.ratioBef, self.ratioAft, self.colors_samples, dec=1,
            scale_correction=corrections[3], z=1)

        plt.subplots_adjust(
            wspace=0.015, hspace=0.15,
            top=0.97, bottom=0.07,
            left=0.045, right=0.965)

        if show:
            plt.show()
        if save:
            dirSave = Path(*Path(self.dataPath[0]).parts[:Path(self.dataPath[0]).parts.index('data') + 1])
            fig.savefig(
                f'{dirSave}' + f'\\{self.fileName}' + ' - Bars plots' + '.png',
                facecolor='w', dpi=600)
            print(f'\n\n· Elastic and viscous moduli chart saved at:\n{dirSave}.')

    def plotHeatMap(self, show=True, save=False):

        def insertKey(keys):
            keys = list(keys)
            index = 1
            while index <= len(keys):
                keys.insert(index, 'Broken')
                index += 2
            return keys

        def drawMap(
                title,
                data_map, frequencies, formulations,
        ):
            fig = plt.figure(figsize=(12, 4.5), facecolor='whitesmoke')
            plt.gca().set_facecolor('gainsboro')
            fig.canvas.manager.set_window_title(self.fileName + f' - {title}')

            colors, decimal, minValue, maxValue = None, None, None, None
            if title == "Elastic recovery (%)":
                for i in range(len(data_map)):
                    for j in range(len(data_map[i])):
                        if data_map[i][j] > 100:
                            data_map[i][j] = None
                colors, decimal, minValue, maxValue = 'RdYlGn', '.0f', 0, 100

            if title == "Loss factor $\\tan(G_0''\,/\,G_0')$":
                for i in range(len(data_map)):
                    for j in range(len(data_map[i])):
                        for k in range(len(data_map[i][j])):
                            if data_map[i][j][k] > 50:
                                data_map[i][j][k] = None
                colors, decimal, minValue, maxValue = 'coolwarm', '.2f', 0.1, 1.9
                data_map = np.array(data_map, dtype=float).flatten().reshape(len(formulations), 20)

            df = pd.DataFrame(data_map, index=formulations, columns=frequencies)
            sns.heatmap(
                df,
                vmin=minValue, vmax=maxValue, cmap=colors,
                annot=True, annot_kws={'size': 8}, fmt=decimal, linewidths=.75, linecolor='whitesmoke')

            plt.title(f"{title} across frequency.")
            plt.tick_params(axis='both', which='both', length=0)
            plt.xticks(rotation=90, ha='center', fontsize=10, color='#383838')
            plt.yticks(rotation=0, ha='right', fontsize=10, color='#383838')

            plt.subplots_adjust(
                wspace=0, hspace=0,
                top=0.93, bottom=0.14,
                left=0.07, right=1.0)

            if show:
                plt.show()
            if save:
                dirSave = Path(*Path(self.dataPath[0]).parts[:Path(self.dataPath[0]).parts.index('data') + 1])
                fig.savefig(f'{dirSave}' + f'\\{title[:3]}' + '.png', facecolor='w', dpi=600)
                print(f'\n\n· Heatmap {title} chart saved at:\n{dirSave}.')

        drawMap(
            "Elastic recovery (%)",
            self.recoveryPCT, self.freqsRecovery, self.names_samples.keys())

        brokenLabels = insertKey(self.names_samples.keys())
        drawMap(
            "Loss factor $\\tan(G_0''\,/\,G_0')$",
            self.tan_delta, self.freqsRecovery, brokenLabels)


class Flow:
    def __init__(
            self,
            dataPath, fileName,
            names_samples, number_samples, colors_samples
    ):

        def getData():
            def getSegments(dataframe):
                time = dataframe['t in s'].to_numpy()
                shear_rate = dataframe['ɣ̇ in 1/s'].to_numpy()
                shear_stress = dataframe['τ in Pa'].to_numpy()

                seg3, seg4, seg5 = (
                    dataframe.index[dataframe['SegIndex'] == seg].to_list()[0] for seg in ['3|1', '4|1', '5|1'])

                segments34 = lambda arr: (arr[seg3:seg4])
                segments43 = lambda arr: (arr[seg4:seg5])
                t_cte = segments34(time)

                return {
                    'time': [t_cte - t_cte[0]],
                    'shear_stress_cte': segments34(shear_stress),
                    'shear_rate': segments43(shear_rate),
                    'shear_stress_step': segments43(shear_stress),
                }

            def dict_FlowShearing(labels):
                return {label: ([], []) for label in labels}

            self.sample_keys = list(self.names_samples.keys())
            if len(self.sample_keys) != len(self.number_samples):
                raise ValueError('The length of "number_samples" must match the number of sample keys.')
            sample_labels = [
                key for key, count in zip(self.sample_keys, self.number_samples) for _ in range(count)
            ]

            for sample_type, path in zip(sample_labels, self.dataPath):
                df = pd.read_excel(path)
                segments = getSegments(df)
                self.names_samples[sample_type].append(segments)

            dict_data_cte = {}
            for sample_type in self.names_samples:
                dict_data_cte[f'{sample_type} time'] = [s['time'] for s in self.names_samples[sample_type]]
                dict_data_cte[f'{sample_type} shear_stress_cte'] = [s['shear_stress_cte'] for s in
                                                                    self.names_samples[sample_type]]

            dict_data_steps = {}
            for sample_type in self.names_samples:
                dict_data_steps[f'{sample_type} shear_rate'] = \
                    [downsampler(s['shear_rate'], 15) for s in self.names_samples[sample_type]]
                dict_data_steps[f'{sample_type} shear_stress_step'] = \
                    [downsampler(s['shear_stress_step'], 15) for s in self.names_samples[sample_type]]

            return (dict_data_cte, dict_data_steps,
                    dict_FlowShearing(self.sample_keys), dict_FlowShearing(self.sample_keys))

        def appendData():
            for key, (t, ss) in self.cteShearRate.items():
                t.append(self.cteData[f'{key} time'])
                ss.append(self.cteData[f'{key} shear_stress_cte'])

            for key, (sr, ss) in self.stepShearRate.items():
                sr.append(self.stepData[f'{key} shear_rate'])
                ss.append(self.stepData[f'{key} shear_stress_step'])

        # input vars
        self.dataPath = dataPath
        self.fileName = fileName
        self.names_samples = names_samples
        self.number_samples = number_samples
        self.colors_samples = colors_samples

        # data vars
        self.tableCteSS, self.tableStepSS = [], []

        # data reading
        self.cteData, self.stepData, self.cteShearRate, self.stepShearRate = getData()
        appendData()

        # chart config
        fonts('C:/Users/petrus.kirsten/AppData/Local/Microsoft/Windows/Fonts/')
        plt.style.use('seaborn-v0_8-ticks')

    def plotShearFlow(
            self,
            cteTitle, cteLimits,
            stepTitle, stepLimits,
            show=True, save=False
    ):

        def getSplitMean():
            time, stress, stressErr = (
                self.cteShearRate[key][0][0],
                self.cteShearRate[key][1][0],
                self.cteShearRate[key][1][0])

            timeSplit, stressSplit, stressErrSplit = [], [], []
            for t in range(len(time)):
                index = np.where(time[t][0] <= 175)[-1][-1]
                timeSplit.append(downsampler(time[t][0][:index]))
                stressSplit.append(downsampler(stress[t][:index]))
                stressErrSplit.append(downsampler(stressErr[t][:index]))

            timeMean, stressMean, stressErrMean = (
                np.mean(timeSplit, axis=0),
                np.mean(stressSplit, axis=0),
                np.std(stressErrSplit, axis=0))

            return timeMean, stressMean, stressErrMean

        def getMean():

            return (
                np.mean(self.stepShearRate[key][0], axis=1)[0],
                np.mean(self.stepShearRate[key][1], axis=1)[0],
                np.std(self.stepShearRate[key][1], axis=1)[0]
            )

        def drawCteSS(listRows, sampleName, ax,
                      x, y, yErr,
                      axTitle, yLabel, yLim, xLabel, xLim,
                      curveColor,
                      fit='', logScale=False):

            def funcTransient(t, tau_0, tau_e, time_cte):
                return tau_e + (tau_0 - tau_e) * np.exp(- t / time_cte)

            def legendLabel():
                legend = ax.legend(
                    fancybox=False, frameon=True,
                    framealpha=0.9, fontsize=9, ncols=2,
                    loc='upper right'
                )
                legend.get_frame().set_facecolor('w')
                legend.get_frame().set_edgecolor('whitesmoke')

            def configPlot():
                ax.set_title(axTitle, size=10, color='k')
                ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(0.75)

                ax.set_xlabel(f'{xLabel}')
                ax.set_xscale('log' if logScale else 'linear')
                ax.set_xlim(xLim)
                ax.xaxis.set_major_locator(MultipleLocator(60))
                ax.xaxis.set_minor_locator(MultipleLocator(10))

                ax.set_ylabel(f'{yLabel}')
                ax.set_yscale('log' if logScale else 'linear')
                ax.set_ylim(yLim)
                ax.yaxis.set_major_locator(MultipleLocator(yLim[1] / 10))
                ax.yaxis.set_minor_locator(MultipleLocator(yLim[1] / 20))

            params, covariance = curve_fit(funcTransient, x, y)
            # p0=(x[0], y[-1], 100))  # method='trf')  # method='dogbox', maxfev=5000)
            errors = np.sqrt(np.diag(covariance))

            x_fit = np.linspace(0, 300, 300)
            y_fit = funcTransient(x_fit, *params)

            listRows = exportFit(
                f'{sampleName}',
                params, errors,
                listRows, mode='thixo')

            ax.errorbar(
                x[::2] if sampleName == 'St + iCar' else x[::4],
                y[::2] if sampleName == 'St + iCar' else y[::4],
                yerr=yErr[::2] if sampleName == 'St + iCar' else yErr[::4],
                color=curveColor, alpha=.85,
                fmt='none', mfc=curveColor,
                capsize=2.5, capthick=1, lw=1, linestyle='',
                label=f'', zorder=2)

            ax.errorbar(
                x[::2] if sampleName == 'St + iCar' else x[::4],
                y[::2] if sampleName == 'St + iCar' else y[::4],
                yerr=0,
                color=curveColor, alpha=.65,
                fmt='D' if 'CL' in sampleName else 'o',
                markersize=6.5 if 'CL' in sampleName else 7,
                mfc=curveColor, mec='#383838', mew=.75,
                linestyle='',
                label=f'{sampleName}', zorder=3)

            ax.plot(
                x_fit, y_fit, color=curveColor, linestyle='--', linewidth=1,
                zorder=4)

            configPlot()
            legendLabel()

            return listRows

        def drawStepSS(listRows, sampleName, ax,
                       x, y, yErr,
                       axTitle, yLabel, yLim, xLabel, xLim,
                       curveColor,
                       fit='', logScale=False):

            def funcHB(sigma, k, n, sigmaZero):
                return sigmaZero + k * (sigma ** n)

            def legendLabel():
                legend = ax.legend(
                    fancybox=False, frameon=False,
                    framealpha=0.9, fontsize=9, ncols=2,
                    loc='upper left')
                legend.get_frame().set_facecolor('w')
                legend.get_frame().set_edgecolor('whitesmoke')

            def configPlot():
                ax.set_title(axTitle, size=10, color='k')
                ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(0.75)

                ax.set_xlabel(f'{xLabel}')
                ax.set_xscale('log' if logScale else 'linear')
                ax.set_xlim(xLim)

                ax.set_ylabel(f'{yLabel}')
                ax.set_yscale('log' if logScale else 'linear')
                ax.set_ylim(yLim)
                ax.yaxis.set_major_locator(MultipleLocator(yLim[1] / 10))
                ax.yaxis.set_minor_locator(MultipleLocator(yLim[1] / 20))

            configPlot()

            if fit == 'HB':
                params, covariance = curve_fit(
                    funcHB, x, y,
                    p0=(2, 1, 0))
                errors = np.sqrt(np.diag(covariance))
                x_fit = np.linspace(.1, 1000, 1000)
                y_fit = funcHB(x_fit, *params)
                listRows = exportFit(
                    f'{sampleName}',
                    params, errors,
                    listRows, mode='HB')

                ax.plot(
                    x_fit, y_fit, color=curveColor, linestyle=':', linewidth=1,
                    zorder=2)

            ax.errorbar(
                x, y, yerr=yErr,
                color=curveColor, alpha=.85,
                fmt='none', mfc=curveColor,
                capsize=2.5, capthick=1, lw=1, linestyle='',
                label=f'', zorder=2)

            ax.errorbar(
                x, y, yerr=0,
                color=curveColor, alpha=.65,
                fmt='D' if 'CL' in sampleName else 'o',
                markersize=7 if 'CL' in sampleName else 6.5,
                mfc=curveColor, mec='#383838', mew=.75,
                linestyle='',
                label=f'{sampleName}', zorder=3)

            legendLabel()

            return listRows

        fig, axes = plt.subplots(
            figsize=(16, 7), ncols=2, nrows=1,
            gridspec_kw={'width_ratios': [1, 1]}, facecolor='snow')
        fig.canvas.manager.set_window_title(self.fileName + ' - Flow shearing')
        axCteSS, axStepSS = axes[0], axes[1]

        fig.suptitle(f'')
        timeTitle, timeLimits = (f'Time (s)', (0, 180))
        rateTitle, rateLimits = ('Shear rate ($s^{-1}$)', (0, 315))

        for key, color in zip(self.cteShearRate, self.colors_samples):
            time_plot, stress_plot, stressErr_plot = getSplitMean()
            self.tableCteSS = drawCteSS(
                listRows=self.tableCteSS, ax=axCteSS,
                x=time_plot, y=stress_plot, yErr=stressErr_plot,
                axTitle='Constant shear rate', yLabel=cteTitle, yLim=cteLimits,
                xLabel=timeTitle, xLim=timeLimits,
                curveColor=color,
                sampleName=f'{key}', fit='transient')

            shear_plot, stress_plot, stressErr_plot = getMean()
            self.tableStepSS = drawStepSS(
                listRows=self.tableStepSS, ax=axStepSS,
                x=shear_plot, y=stress_plot, yErr=stressErr_plot,
                axTitle='Steps shear rate', yLabel=stepTitle, yLim=stepLimits,
                xLabel=rateTitle, xLim=rateLimits,
                curveColor=color,
                sampleName=f'{key}', fit='HB')

        # plt.subplots_adjust(
        #     hspace=0, wspace=0.21,
        #     top=0.92, bottom=0.075,
        #     left=0.045, right=0.96)
        plt.tight_layout()
        if show:
            plt.show()
        if save:
            dirSave = Path(*Path(self.dataPath[0]).parts[:Path(self.dataPath[0]).parts.index('data') + 1])
            fig.savefig(
                f'{dirSave}' + f'\\{self.fileName}' + ' - Flow shearing' + '.png',
                facecolor='w', dpi=600)
            print(f'\n\n· Flow shearing charts saved at:\n{dirSave}.')

    def plotFits(
            self,
            cteLimits, stepLimits,
            show=True, save=False
    ):

        def drawThixoData(
                title,
                axes, data,
                colors, a, z
        ):
            def configPlot(ax, yTitle, yLim):
                ax.set_title(title, size=10, color='k')

                ax.tick_params(axis='x', labelsize=10, length=4)
                ax.tick_params(axis='y', which='both', labelsize=9, pad=1, length=0)

                ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(.75)
                ax.spines[['top', 'bottom', 'left', 'right']].set_color('#303030')

                ax.set_xticks([]), ax.set_xlim([x[0] - 2, x[-1]]), ax.set_yticks([]), ax.set_ylim(yLim)

            axes3 = axes.twinx()

            bin_width, space_samples = 0.8, 3
            x = np.arange(space_samples * len(data))

            configPlot(axes, "$\\tau_0$ (Pa) and $\\tau_e$ (Pa)", (0, cteLimits[0]))
            configPlot(axes3, "$\lambda$ (s)", (0, cteLimits[1]))

            samples = [d['Sample'] for d in data]
            kPrime, kPrime_err = [d["$\\tau_0$"] for d in data], [d["± $\\tau_0$"] for d in data]
            nPrime, nPrime_err = [d["$\\tau_e$"] for d in data], [d["± $\\tau_e$"] for d in data]
            sigmaZero, sigmaZero_err = [d["$\lambda$"] for d in data], [d["± $\lambda$"] for d in data]

            posList, labelsList = [], []
            for i in range(len(kPrime)):
                axes.bar(
                    space_samples * x[i] - bin_width,
                    height=kPrime[i], yerr=0,
                    color=colors[i], edgecolor='#383838',
                    width=bin_width, hatch='////', alpha=a, linewidth=.5,
                    zorder=z)
                axes.errorbar(
                    x=space_samples * x[i] - bin_width, y=kPrime[i], yerr=kPrime_err[i],
                    color='#383838', alpha=.99, linewidth=1, capsize=4, capthick=1.05,
                    zorder=3)
                axes.text(
                    space_samples * x[i] - bin_width - .15,
                    kPrime[i] + kPrime_err[i] + cteLimits[0] * .075,
                    f'{kPrime[i]:.{1}f} ± {kPrime_err[i]:.{1}f}',
                    va='center', ha='left', rotation=90,
                    color='#383838', fontsize=9)

                axes.bar(
                    space_samples * x[i],
                    height=nPrime[i], yerr=0,
                    color=colors[i], edgecolor='#383838',
                    width=bin_width, hatch='....', alpha=a, linewidth=.5,
                    zorder=z)
                axes.errorbar(
                    x=space_samples * x[i], y=nPrime[i], yerr=nPrime_err[i],
                    color='#383838', alpha=.99, linewidth=1, capsize=4, capthick=1.05,
                    zorder=3)
                axes.text(
                    space_samples * x[i] - .15,
                    nPrime[i] + nPrime_err[i] + cteLimits[0] * .075,
                    f'{nPrime[i]:.{1}f} ± {nPrime_err[i]:.{1}f}',
                    va='center', ha='left', rotation=90,
                    color='#383838', fontsize=9)

                axes3.bar(
                    space_samples * x[i] + bin_width,
                    height=sigmaZero[i], yerr=0,
                    color=colors[i], edgecolor='#383838',
                    width=bin_width, hatch='', alpha=a, linewidth=.5,
                    zorder=z)
                axes3.errorbar(
                    x=space_samples * x[i] + bin_width, y=sigmaZero[i], yerr=sigmaZero_err[i],
                    color='#383838', alpha=.99, linewidth=1, capsize=4, capthick=1.05,
                    zorder=3)
                axes3.text(
                    space_samples * x[i] + bin_width - .15,
                    sigmaZero[i] + sigmaZero_err[i] + cteLimits[1] * .075,
                    f'{sigmaZero[i]:.{1}f} ± {sigmaZero_err[i]:.{1}f}',
                    va='center', ha='left', rotation=90,
                    color='#383838', fontsize=9)

                posList.append(space_samples * x[i] - bin_width), posList.append(space_samples * x[i]), posList.append(
                    space_samples * x[i] + bin_width)
                labelsList.append("$\\tau_0$"), labelsList.append("$\\tau_e$"), labelsList.append("$\lambda$")

            axes.set_xticks(posList)
            axes.set_xticklabels(labelsList)

            return nPrime

        def drawHBdata(
                title, axes,
                data, colors, a, z
        ):
            def configPlot(ax, yTitle, yLim):
                ax.set_title(title, size=10, color='k')
                if yTitle == "$k'$":
                    ax.grid(which='major', axis='y', linestyle='-', linewidth=1, color='lightgray', alpha=0.5,
                            zorder=-1)
                    ax.grid(which='minor', axis='y', linestyle='--', linewidth=.75, color='lightgray', alpha=0.5,
                            zorder=-1)

                ax.tick_params(axis='x', labelsize=10, length=4)
                ax.tick_params(
                    axis='y', which='both', labelsize=9, pad=1, length=0,
                    labeltop=False, top=False,
                    labelbottom=False, bottom=False)

                ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(.75)
                ax.spines[['top', 'bottom', 'left', 'right']].set_color('#303030')

                ax.set_xticks([]), ax.set_xlim([x[0] - 2, x[-1]]), ax.set_yticks([]), ax.set_ylim(yLim)

            axes2, axes3 = axes.twinx(), axes.twinx()

            bin_width, space_samples = 0.8, 3
            x = np.arange(space_samples * len(data))

            configPlot(axes, "$k'$", (0, stepLimits[0]))
            configPlot(axes2, "$n'$", (0, stepLimits[1]))
            configPlot(axes3, "$\sigma_0$ (Pa)", (0, stepLimits[2]))

            samples = [d['Sample'] for d in data]
            kPrime, kPrime_err = [abs(d["k'"]) for d in data], [d["± k'"] for d in data]
            nPrime, nPrime_err = [abs(d["n'"]) for d in data], [d["± n'"] for d in data]
            sigmaZero, sigmaZero_err = [abs(d["sigma_zero"]) for d in data], [d["± sigma_zero"] for d in data]

            posList, labelsList = [], []

            for i in range(len(kPrime)):
                axes.bar(
                    space_samples * x[i] - bin_width,
                    height=kPrime[i], yerr=0,
                    color=colors[i], edgecolor='#383838',
                    width=bin_width, hatch='////', alpha=a, linewidth=.5,
                    zorder=z)
                axes.errorbar(
                    x=space_samples * x[i] - bin_width, y=kPrime[i], yerr=kPrime_err[i],
                    color='#383838', alpha=.99, linewidth=1, capsize=4, capthick=1.05,
                    zorder=3)
                axes.text(
                    space_samples * x[i] - bin_width - .15,
                    kPrime[i] + kPrime_err[i] + stepLimits[0] * .075,
                    f'{kPrime[i]:.{2}f} ± {kPrime_err[i]:.{2}f}',
                    va='center', ha='left', rotation=90,
                    color='#383838', fontsize=9)

                axes2.bar(
                    space_samples * x[i],
                    height=nPrime[i], yerr=0,
                    color=colors[i], edgecolor='#383838',
                    width=bin_width, hatch='....', alpha=a, linewidth=.5,
                    zorder=z)
                axes2.errorbar(
                    x=space_samples * x[i], y=nPrime[i], yerr=nPrime_err[i],
                    color='#383838', alpha=.99, linewidth=1, capsize=4, capthick=1.05,
                    zorder=3)
                axes2.text(
                    space_samples * x[i] - .15,
                    nPrime[i] + nPrime_err[i] + stepLimits[1] * .075,
                    f'{nPrime[i]:.{2}f} ± {nPrime_err[i]:.{2}f}',
                    va='center', ha='left', rotation=90,
                    color='#383838', fontsize=9)

                axes3.bar(
                    space_samples * x[i] + bin_width,
                    height=sigmaZero[i], yerr=0,
                    color=colors[i], edgecolor='#383838',
                    width=bin_width, hatch='', alpha=a, linewidth=.5,
                    zorder=z)
                axes3.errorbar(
                    x=space_samples * x[i] + bin_width, y=sigmaZero[i], yerr=sigmaZero_err[i],
                    color='#383838', alpha=.99, linewidth=1, capsize=4, capthick=1.05,
                    zorder=3)
                axes3.text(
                    space_samples * x[i] + bin_width - .15,
                    sigmaZero[i] + sigmaZero_err[i] + stepLimits[2] * .075,
                    f'{sigmaZero[i]:.{1}f} ± {sigmaZero_err[i]:.{1}f}',
                    va='center', ha='left', rotation=90,
                    color='#383838', fontsize=9)

                posList.append(space_samples * x[i] - bin_width), posList.append(space_samples * x[i]), posList.append(
                    space_samples * x[i] + bin_width)
                labelsList.append("$k'$"), labelsList.append("$n'$"), labelsList.append("$\sigma_0$")

            axes.set_xticks(posList)
            axes.set_xticklabels(labelsList)

            return nPrime

        fig, axs = plt.subplots(
            figsize=(16, 7), ncols=2, nrows=1,
            gridspec_kw={'width_ratios': [1, 1]}, facecolor='snow')
        fig.canvas.manager.set_window_title(self.fileName + ' - Flow shearing fit parameters')
        axCteSS, axStepSS = axs[0], axs[1]

        _ = drawThixoData(
            '',
            axCteSS, self.tableCteSS,
            self.colors_samples, a=.85, z=2)

        _ = drawHBdata(
            '',
            axStepSS, self.tableStepSS,
            self.colors_samples, a=.85, z=2)

        plt.tight_layout()
        if show:
            plt.show()
        if save:
            dirSave = Path(*Path(self.dataPath[0]).parts[:Path(self.dataPath[0]).parts.index('data') + 1])
            fig.savefig(
                f'{dirSave}' + f'\\{self.fileName}' + ' - Flow shearing fit parameters' + '.png',
                facecolor='w', dpi=600)
            print(f'\n\n· Flow fit parameters charts saved at:\n{dirSave}.')


class Compression:
    def __init__(
            self,
            dataPath, fileName,
            names_samples, number_samples, colors_samples
    ):

        def getData(cut=500):

            def getSegments(dataframe, segInit, segEnd):
                time = dataframe['t in s'].to_numpy()
                height = dataframe['h in mm'].to_numpy()
                force = dataframe['Fn in N'].to_numpy()

                indexInit, indexEnd = (
                    dataframe.index[dataframe['SegIndex'] == seg].to_list()[0]
                    for seg in [segInit, segEnd])
                segments = lambda arr: (arr[indexInit:indexEnd])

                return {
                    'time': segments(time) - np.min(segments(time)),
                    'height': segments(height),
                    'force': segments(force)}

            def dict_Compression(labels):
                return {label: ([], [], []) for label in labels}

            self.sample_keys = list(self.names_samples.keys())
            if len(self.sample_keys) != len(self.number_samples):
                raise ValueError('The length of "number_samples" must match the number of sample keys.')
            sample_labels = [
                key for key, count in zip(self.sample_keys, self.number_samples) for _ in range(count)
            ]

            for sample_type, path in zip(sample_labels, self.dataPath):
                df = pd.read_excel(path)
                segments = getSegments(
                    df,
                    segInit='2|1' if '171024' in path or not 'kC-compression-4' in path else '1|1',
                    segEnd='62|1' if '10St_CL_7' in path else '61|1')
                self.names_samples[sample_type].append(segments)

            dict_data_dynamic = {}
            for sample_type in self.names_samples:
                dict_data_dynamic[f'{sample_type} time'] = \
                    [downsampler(s['time'], cut) for s in self.names_samples[sample_type]]
                dict_data_dynamic[f'{sample_type} height'] = \
                    [downsampler(s['height'], cut) for s in self.names_samples[sample_type]]
                dict_data_dynamic[f'{sample_type} force'] = \
                    [downsampler(s['force'], cut) for s in self.names_samples[sample_type]]

            # dict_data_steps = {}
            # for sample_type in self.names_samples:
            #     dict_data_steps[f'{sample_type} shear_rate'] = \
            #         [downsampler(s['shear_rate'], 15) for s in self.names_samples[sample_type]]
            #     dict_data_steps[f'{sample_type} shear_stress_step'] = \
            #         [downsampler(s['shear_stress_step'], 15) for s in self.names_samples[sample_type]]

            return dict_data_dynamic, dict_Compression(self.sample_keys)

        def appendData():
            for key, (x, height, f) in self.listDynamic.items():
                x.append(self.dynamicData[f'{key} time'])
                height.append(self.dynamicData[f'{key} height'])
                f.append(self.dynamicData[f'{key} force'])

            # for key, (sr, ss) in self.stepShearRate.items():
            #     sr.append(self.stepData[f'{key} shear_rate'])
            #     ss.append(self.stepData[f'{key} shear_stress_step'])

        # input vars
        self.dataPath = dataPath
        self.fileName = fileName
        self.names_samples = names_samples
        self.number_samples = number_samples
        self.colors_samples = colors_samples

        # data vars
        self.tableDynamic, self.tableCompression = [], []
        self.means_hMax = []

        # data reading
        self.dynamicData, self.listDynamic = getData()
        appendData()

        # chart config
        fonts('C:/Users/petrus.kirsten/AppData/Local/Microsoft/Windows/Fonts/')
        plt.style.use('seaborn-v0_8-ticks')
        self.figGraphs = plt.figure(figsize=(16, 9), facecolor='snow')
        self.gsGraphs = GridSpec(2, 2, height_ratios=[1, 1.5], width_ratios=[1, 1])
        self.figGraphs.canvas.manager.set_window_title(self.fileName + ' - Compression')

    def plotDynamic(
            self,
            sLimits,
            show=True, save=False
    ):

        def getValues():
            time, height, stress = (self.listDynamic[key][0][0],
                                    self.listDynamic[key][1][0],
                                    self.listDynamic[key][2][0])

            hMax_mean = np.max(height, axis=1)
            self.means_hMax.append(hMax_mean.tolist())

            return (np.mean(time, axis=0),
                    np.mean(stress, axis=0) / (35 / 1000) + 7.5,
                    np.std(stress, axis=0) / (35 / 1000))

        def drawDynamicStress(
                sampleName,
                ax, x, y, yErr,
                axTitle, yLabel, yLim, xLabel, xLim, axisColor,
                curveColor, markerStyle, markerFColor, markerEColor, markerEWidth=0.5,
                strain=False, lineStyle='', logScale=False
        ):
            def legendLabel():
                legend = ax.legend(fancybox=False, frameon=True, framealpha=0.9, fontsize=9)
                legend.get_frame().set_facecolor('w')
                legend.get_frame().set_edgecolor('whitesmoke')

            def configPlot():
                ax.set_title(axTitle, size=10, color='crimson')
                ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(0.75)
                ax.grid(True, which='both', axis='x', linestyle='--', linewidth=0.5, color='silver', alpha=0.5)

                ax.set_xlabel(f'{xLabel}')
                # ax.set_xscale('log' if logScale else 'linear')
                ax.set_xlim(xLim)
                ax.xaxis.set_minor_locator(MultipleLocator(2))

                ax.set_ylabel(f'{yLabel}', color=axisColor)
                ax.set_yscale('log' if logScale else 'linear')
                ax.set_ylim(yLim)
                ax.tick_params(axis='y', colors=axisColor, which='both')

            configPlot()

            if strain:
                ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}%"))
                ax.plot(
                    x[0], y[0],
                    color=curveColor, alpha=0.8, lw=1.5, linestyle=':',
                    label=f'{sampleName}', zorder=3)

            else:
                configPlot()

                x_smooth = np.linspace(x.min(), x.max(), len(x) * 5)
                interp_yErr_lower, interp_yErr_upper = interp1d(x, y - yErr, kind='cubic'), interp1d(x, y + yErr,
                                                                                                     kind='cubic')
                yErr_lower_smooth, yErr_upper_smooth = interp_yErr_upper(x_smooth), interp_yErr_lower(x_smooth)

                ax.fill_between(
                    x_smooth, yErr_lower_smooth, yErr_upper_smooth,
                    color=curveColor, alpha=0.15, lw=0,
                    zorder=6 - self.colors_samples.index(color))
                ax.errorbar(
                    x, y, 0,
                    color=curveColor, alpha=.35,
                    fmt=markerStyle, markersize=3.5, mfc=markerFColor, mec=markerEColor, mew=markerEWidth,
                    capsize=0, lw=.5, linestyle=lineStyle,
                    label=f'{sampleName}',
                    zorder=7 - self.colors_samples.index(color))
                legendLabel()

        def drawCycles(
                sampleName,
                ax, y, yErr,
                axTitle, yLabel, yLim, xLabel, xLim, axisColor,
                curveColor, markerFColor, markerEColor, paramsList, markerEWidth=0.75,
                lineStyle='', logScale=False
        ):

            def funcRelaxation(t, sigma_0, sigma_e, time_cte):
                return sigma_e + (sigma_0 - sigma_e) * np.exp(- t / time_cte)

            def configPlot():
                ax.set_title(axTitle, size=10, color='crimson')
                ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(0.75)
                # ax.grid(True, which='both', axis='x', linestyle='--', linewidth=0.5, color='silver', alpha=0.5)

                ax.set_xlabel(f'{xLabel}')
                # ax.set_xscale('log' if logScale else 'linear')
                ax.set_xlim(xLim)
                ax.xaxis.set_minor_locator(MultipleLocator(1))

                ax.set_ylabel(f'{yLabel}', color=axisColor)
                ax.set_yscale('log' if logScale else 'linear')
                ax.set_ylim(yLim)
                ax.tick_params(axis='y', colors=axisColor, which='both')

            # reshape to each cycle
            target_size = 29 * (len(y) // 29)
            y, yErr = y[:target_size], yErr[:target_size]
            y, yErr = y.reshape(29, -1), yErr.reshape(29, -1)

            # find and append maxs and mins
            listYmax, listYmaxErr = [], []
            listYmin, listYminErr = [], []
            for cycle in range(len(y)):
                yMax, yMin = y[cycle][np.argmax(y[cycle])], y[cycle][np.argmin(y[cycle])]
                yMaxErr, yMinErr = yErr[cycle][np.argmax(y[cycle])], yErr[cycle][np.argmin(y[cycle])]

                listYmax.append(yMax), listYmaxErr.append(yMaxErr)
                listYmin.append(yMin), listYminErr.append(yMinErr)
            yMax_values, yMin_values = np.array(listYmax), np.array(listYmin)
            yMax_errors, yMin_errors = np.array(listYmaxErr), np.array(listYminErr)

            x_values = np.arange(1, len(listYmax) + 1)

            # Fitting
            params, covariance = curve_fit(
                funcRelaxation, x_values, yMax_values,
                p0=(yMax_values[0], yMax_values[-1], 5))  # method='trf')  # method='dogbox', maxfev=5000)
            errors = np.sqrt(np.diag(covariance))
            x_fit = np.linspace(0, 30, 120)
            y_fit = funcRelaxation(x_fit, *params)

            paramsList.append([params.tolist(), errors.tolist()])

            configPlot()

            ax.errorbar(
                x_values, yMax_values, yMax_errors,
                color=curveColor, alpha=.85,
                fmt='none', mfc=curveColor,
                capsize=2.5, capthick=1, linestyle='', lw=1,
                label=f'', zorder=5 - self.colors_samples.index(color))

            ax.errorbar(
                x_values, yMax_values, 0,
                color=curveColor, alpha=.65,
                fmt='^', markersize=6.5,
                mfc=curveColor, mec='#383838', mew=markerEWidth,
                linestyle='',
                label=f'{sampleName}', zorder=6 - self.colors_samples.index(color))

            ax.plot(
                x_fit, y_fit, color=curveColor, linestyle='-.', linewidth=.75,
                zorder=1)
            # ax.errorbar(
            #     x_values, yMin_values, yMin_errors,
            #     color=curveColor, alpha=.65,
            #     fmt='v', markersize=5, mfc=markerFColor, mec=markerEColor, mew=markerEWidth,
            #     capsize=2.5, lw=.75, linestyle=lineStyle,
            #     label=f'{sampleName}',
            #     zorder=4)
            return paramsList

        axOscillation, axCycles, axBars = (self.figGraphs.add_subplot(self.gsGraphs[0, :]),
                                           self.figGraphs.add_subplot(self.gsGraphs[1, 0]),
                                           self.figGraphs.add_subplot(self.gsGraphs[1, 1]))

        s1Title, s1Limits = f'Stress (Pa)', (0, sLimits)
        s2Title, s2Limits = f'Stress peak (Pa)', (0, 210)
        x1Title, x1Limits = f'Time (s)', (0, 75)
        x2Title, x2Limits = f'Cycle', (0, 30)

        for key, color in zip(self.listDynamic, self.colors_samples):
            timeMean, stressMean, stressErrMean = getValues()

            drawDynamicStress(
                ax=axOscillation, axisColor='#303030',
                x=timeMean, y=stressMean, yErr=stressErrMean,
                axTitle='', yLabel=s1Title, yLim=s1Limits, xLabel=x1Title, xLim=x1Limits,
                curveColor=color, markerStyle='o', markerFColor=color, markerEColor='#303030',
                sampleName=f'{key}')

            self.tableDynamic = drawCycles(
                ax=axCycles, axisColor='#303030',
                y=stressMean, yErr=stressErrMean,
                axTitle='', yLabel=s2Title, yLim=s2Limits, xLabel=x2Title, xLim=x2Limits,
                curveColor=color, markerFColor=color, markerEColor=color,
                paramsList=self.tableDynamic, sampleName=f'{key}')

            plt.subplots_adjust(
                wspace=0.015, hspace=0.060,
                top=0.970, bottom=0.070,
                left=0.060, right=0.985)

            if show:
                plt.show()
            if save:
                dirSave = Path(*Path(self.dataPath[0]).parts[:Path(self.dataPath[0]).parts.index('data') + 1])
                self.figGraphs.savefig(
                    f'{dirSave}' + f'\\{self.fileName}' + ' - Dynamic compression' + '.png',
                    facecolor='w', dpi=600)
                print(f'\n\n· Dynamic compression chart saved at:\n{dirSave}.')
