from datetime import datetime

import numpy as np
import pandas as pd
from pathlib import Path

from matplotlib.ticker import MultipleLocator
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


def exportFit(
        sample,
        data, err,
        table):
    keys = ('k', 'n', 'sigma_zero')
    values = (data, err)

    dictFit = {'Sample': sample}
    iParams = 0
    for key, value in zip(keys, range(len(values[0]))):
        dictFit[f'{key}'] = values[0][iParams]
        dictFit[f'± {key}'] = values[1][iParams]
        iParams += 1

    table.append(dictFit)

    return table


def funcHB(sigma, k, n, sigmaZero):
    return sigmaZero + k * (sigma ** n)


def getSamplesInfos(
        # quantity
        n_kc_0, n_kc_7, n_kc_14, n_kc_21, n_kc_28, n_kc_42,
        # colors
        color_kc_0, color_kc_7, color_kc_14, color_kc_21, color_kc_28, color_kc_42,
):
    number_samples = [n_kc_0, n_kc_7, n_kc_14, n_kc_21, n_kc_28, n_kc_42]

    colors_samples = [color_kc_0, color_kc_7, color_kc_14, color_kc_21, color_kc_28, color_kc_42]

    return number_samples, colors_samples


def getSamplesData(dataPath, number_samples):
    def getSegments(dataframe):
        shear_rate = dataframe['ɣ̇ in 1/s'].to_numpy()
        shear_stress = dataframe['τ in Pa'].to_numpy()
        viscosity = dataframe['η in mPas'].to_numpy()

        seg4, seg5 = (dataframe.index[dataframe['SegIndex'] == seg].to_list()[0] for seg in ['4|1', '5|1'])

        segments = lambda arr: (arr[seg4:seg5])

        return {
            'shear rate': segments(shear_rate),
            'shear stress': segments(shear_stress),
            'viscosity': segments(viscosity)
        }

    def downsampler(array, n=15):
        if len(array) > n:
            step = len(array) // n  # Calculate step size
            return array[::step][:n]
        return array

    samples = {
        'kCar': [], 'kCar/CL-7': [], 'kCar/CL-14': [],
        'kCar/CL-21': [], 'kCar/CL-28': [], 'kCar/CL-42': []
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

    dict_stepsRate = {}
    for sample_type in samples:
        dict_stepsRate[f'{sample_type} rate'] = [downsampler(s['shear rate']) for s in samples[sample_type]]
        dict_stepsRate[f'{sample_type} stress'] = [downsampler(s['shear stress']) for s in samples[sample_type]]
        dict_stepsRate[f'{sample_type} viscosity'] = [downsampler(s['viscosity']) for s in samples[sample_type]]

    return dict_stepsRate, list(samples.keys())


def plotFlow(listRows, sampleName,
             ax, x, y, yErr,
             axTitle, yLabel, yLim, xLabel, xLim,
             curveColor, markerStyle,
             fit='', logScale=False):
    def legendLabel():
        legend = ax.legend(
            fancybox=False, frameon=False,
            framealpha=0.9, fontsize=9, ncols=2,
            loc='upper left')
        legend.get_frame().set_facecolor('w')
        legend.get_frame().set_edgecolor('whitesmoke')

    def configPlot():
        ax.set_title(axTitle, size=9, color='crimson')
        ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(0.75)

        ax.set_xlabel(f'{xLabel}')
        ax.set_xscale('log' if logScale else 'linear')
        ax.set_xlim(xLim)

        ax.set_ylabel(f'{yLabel}')
        ax.set_yscale('log' if logScale else 'linear')
        ax.set_ylim(yLim)
        ax.yaxis.set_major_locator(MultipleLocator(10))
        ax.yaxis.set_minor_locator(MultipleLocator(2.5))

    configPlot()

    if fit == 'HB':
        params, covariance = curve_fit(
            funcHB, x, y,
            p0=(2, 1, 0))
        errors = np.sqrt(np.diag(covariance))
        K, n, sigmaZero = params
        x_fit = np.linspace(.1, 1000, 1000)
        y_fit = funcHB(x_fit, K, n, sigmaZero)
        listRows = exportFit(
            f'{sampleName}',
            params, errors,
            listRows)

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
        fmt='D' if '21' in sampleName else 'o',
        markersize=7 if '21' in sampleName else 6.5,
        mfc=curveColor, mec='#383838', mew=.75,
        linestyle='',
        label=f'{sampleName}', zorder=3)

    legendLabel()

    return listRows


def plotBars(title, axes, data, colors, a, z):
    def configPlot(ax, yTitle, yLim):
        ax.set_title(title, size=10, color='k')
        if yTitle == "$k'$":
            ax.grid(which='major', axis='y', linestyle='-', linewidth=1, color='lightgray', alpha=0.5, zorder=-1)
            ax.grid(which='minor', axis='y', linestyle='--', linewidth=.75, color='lightgray', alpha=0.5, zorder=-1)

        ax.tick_params(axis='x', labelsize=10, length=4)
        ax.tick_params(
            axis='y', which='both', labelsize=9, pad=1, length=0,
            labeltop=False, top=False,
            labelbottom=False, bottom=False)

        ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(.75)
        ax.spines[['top', 'bottom', 'left', 'right']].set_color('#303030')

        ax.set_xticks([]), ax.set_xlim([-2, 17]), ax.set_yticks([]), ax.set_ylim(yLim)

    axes2, axes3 = axes.twinx(), axes.twinx()

    lim1, lim2, lim3 = (0, 50), (0, 2), (0, 25)

    configPlot(axes, "$k'$", lim1)
    configPlot(axes2, "$n'$", lim2)
    configPlot(axes3, "$\sigma_0$ (Pa)", lim3)

    samples = [d['Sample'] for d in data]
    kPrime, kPrime_err = [abs(d["k"]) for d in data], [d["± k"] for d in data]
    nPrime, nPrime_err = [abs(d["n"]) for d in data], [d["± n"] for d in data]
    sigmaZero, sigmaZero_err = [abs(d["sigma_zero"]) for d in data], [d["± sigma_zero"] for d in data]

    bin_width, space_samples = 0.8, 3
    x = np.arange(space_samples * len(data))

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
            kPrime[i] + kPrime_err[i] + lim1[1]*.075,
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
            nPrime[i] + nPrime_err[i] + lim2[1]*.075,
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
            sigmaZero[i] + sigmaZero_err[i] + lim3[1]*.075,
            f'{sigmaZero[i]:.{1}f} ± {sigmaZero_err[i]:.{1}f}',
            va='center', ha='left', rotation=90,
            color='#383838', fontsize=9)

        posList.append(space_samples * x[i] - bin_width), posList.append(space_samples * x[i]), posList.append(space_samples * x[i] + bin_width)
        labelsList.append("$k'$"), labelsList.append("$n'$"), labelsList.append("$\sigma_0$'")

    axes.set_xticks(posList)
    axes.set_xticklabels(labelsList)

    return nPrime


def main(dataPath, fileName):
    fonts('C:/Users/petrus.kirsten/AppData/Local/Microsoft/Windows/Fonts/')
    plt.style.use('seaborn-v0_8-ticks')

    fig, axes = plt.subplots(
        figsize=(16, 7), ncols=2, nrows=1,
        gridspec_kw={'width_ratios': [1.2, 1]}, facecolor='snow')
    axFlow, axBars = axes[0], axes[1]

    fig.suptitle(f'Steps shear rate flow')
    xTitle, xLimits = ('Shear rate ($s^{-1}$)', (0, 315))
    yTitle, yLimits = (f'Shear stress (Pa)', (0, 160))

    nSamples, colorSamples = getSamplesInfos(
        3, 4, 2,
        3, 2, 4,
        'lightsteelblue', '#A773FF',  '#892F99',
        '#AB247B', '#E64B83', '#FF0831')
    data, labels = getSamplesData(dataPath, nSamples)

    dictData = {
        labels[0]: ([], [], []),
        labels[1]: ([], [], []),
        labels[2]: ([], [], []),
        labels[3]: ([], [], []),
        labels[4]: ([], [], []),
        labels[5]: ([], [], [])
    }

    for key, (x, s, v) in dictData.items():
        x.append(data[f'{key} rate'])
        s.append(data[f'{key} stress'])
        v.append(data[f'{key} viscosity'])

    fitModeStress, fitModeVisc = 'HB', ''

    tableFitting = []
    for key, color in zip(dictData, colorSamples):
        shear, stress, stressErr = (
            np.mean(dictData[key][0], axis=1)[0],
            np.mean(dictData[key][1], axis=1)[0],
            np.std(dictData[key][1], axis=1)[0])

        tableFitting = plotFlow(
            listRows=tableFitting,
            ax=axFlow, x=shear, y=stress, yErr=stressErr,
            axTitle='', yLabel=yTitle, yLim=yLimits, xLabel=xTitle, xLim=xLimits,
            curveColor=color, markerStyle='o',
            sampleName=f'{key}', fit=fitModeStress)

    _ = plotBars(
        '',
        axBars, tableFitting,
        colorSamples, a=.85, z=2)

    plt.subplots_adjust(
        hspace=0, wspace=0.120,
        top=0.940, bottom=0.095,
        left=0.050, right=0.950)
    plt.show()

    dirSave = Path(*Path(filePath[0]).parts[:Path(filePath[0]).parts.index('data') + 1])
    fig.savefig(f'{dirSave}' + f'\\{fileName}' + '.png', facecolor='w', dpi=600)

    # fitParams = pd.DataFrame(tableFitting)
    # fitParams.to_excel(f'{dirSave}' + f'\\{fileName}' + '.xlsx', index=False)

    print(f'\n\n· Chart saved at\n{dirSave}.')


if __name__ == '__main__':
    folderPath = "C:/Users/petrus.kirsten/PycharmProjects/RheometerPlots/data"
    # folderPath = "C:/Users/Petrus Kirsten/Documents/GitHub/RheometerPlots/data"

    filePath = [
        # kC
        folderPath + "/231024/kC/kC-viscoelasticRecovery-1.xlsx",
        folderPath + "/231024/kC/kC-viscoelasticRecovery-2.xlsx",
        folderPath + "/231024/kC/kC-viscoelasticRecovery-3.xlsx",

        # kC CL 7
        folderPath + "/231024/kC_CL/kC_CL-viscoelasticRecovery-1.xlsx",
        folderPath + "/231024/kC_CL/kC_CL-viscoelasticRecovery-2.xlsx",
        folderPath + "/231024/kC_CL/kC_CL-viscoelasticRecovery-3.xlsx",
        folderPath + "/231024/kC_CL/kC_CL-viscoelasticRecovery-4.xlsx",

        # kC CL 14
        folderPath + "/071124/kC_CL_14/kC_CL_14-viscoelasticRecovery-1.xlsx",
        folderPath + "/071124/kC_CL_14/kC_CL_14-viscoelasticRecovery-2.xlsx",

        # kC CL 21
        folderPath + "/071124/kC_CL_21/kC_CL_21-viscoelasticRecovery-1.xlsx",
        # folderPath + "/071124/kC_CL_21/kC_CL_21-viscoelasticRecovery-2.xlsx",
        folderPath + "/071124/kC_CL_21/kC_CL_21-viscoelasticRecovery-3.xlsx",

        # kC CL 28
        folderPath + "/071124/kC_CL_28/kC_CL_28-viscoelasticRecovery-1.xlsx",
        folderPath + "/071124/kC_CL_28/kC_CL_28-viscoelasticRecovery-2.xlsx",

        # kC CL 42
        folderPath + "/071124/kC_CL_42/kC_CL_42-viscoelasticRecovery-1.xlsx",
        folderPath + "/071124/kC_CL_42/kC_CL_42-viscoelasticRecovery-2.xlsx",
        # folderPath + "/071124/kC_CL_42/kC_CL_42-viscoelasticRecovery-3.xlsx",
        # folderPath + "/071124/kC_CL_42/kC_CL_42-viscoelasticRecovery-4.xlsx",

    ]

    main(dataPath=filePath, fileName='kappaCar-Flow')
