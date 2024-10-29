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
        """Applies consistent styling to legends in plots."""
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
        ax.yaxis.set_minor_locator(MultipleLocator(25))

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
        fmt=markerStyle, markersize=7, mec='k', mew=0.5,
        capsize=3, lw=.5, linestyle='',  # ecolor='k'
        label=f'{sampleName}', zorder=3)

    legendLabel()

    return listRows


def plotBars(title, axes, data, colors, a, z):
    def configPlot(ax, yTitle, yLim):
        ax.set_title(title, size=10, color='k')
        if yTitle == "$k'$":
            ax.grid(which='major', axis='y', linestyle='-', linewidth=1, color='lightgray', alpha=0.5, zorder=-1)
            ax.grid(which='minor', axis='y', linestyle='--', linewidth=.75, color='lightgray', alpha=0.5, zorder=-1)

        if yTitle == "$n'$":
            ax.yaxis.tick_left()
            ax.yaxis.set_label_position('left')
            ax.tick_params(axis='y', which='both', direction='in', pad=-25)
            ax.yaxis.set_label_coords(0.1, 0.5)

        ax.tick_params(axis='x', labelsize=10, length=4)
        # ax.tick_params(axis='y', which='both', direction='out', pad=1)

        ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(.75)
        ax.spines[['top', 'bottom', 'left', 'right']].set_color('#303030')

        ax.set_xticks([])
        ax.set_xlim([-3, 23])

        ax.set_ylabel(yTitle)
        ax.set_ylim(yLim)

        ax.yaxis.set_major_locator(MultipleLocator(yLim[1] / 5))
        ax.yaxis.set_minor_locator(MultipleLocator(yLim[1] / 20))

    axes2, axes3 = axes.twinx(), axes.twinx()

    configPlot(axes, "$k'$", (0, 75))
    configPlot(axes2, "$n'$", (0, 2))
    configPlot(axes3, "$\sigma_0 (Pa)$", (0, 120))

    samples = [d['Sample'] for d in data]
    kPrime, kPrime_err = [d["k"] for d in data], [d["± k"] for d in data]
    nPrime, nPrime_err = [d["n"] for d in data], [d["± n"] for d in data]
    sigmaZero, sigmaZero_err = [d["sigma_zero"] for d in data], [d["± sigma_zero"] for d in data]

    w, s = 0.8, 3
    x = np.arange(s * len(data))

    posList, labelsList = [], []

    for i in range(len(kPrime)):
        axes.bar(
            s * x[i] - w,
            height=kPrime[i], yerr=0,
            color=colors[i], edgecolor='#383838',
            width=w, hatch='////', alpha=a, linewidth=.5,
            zorder=z)
        axes.errorbar(
            x=s * x[i] - w, y=kPrime[i], yerr=kPrime_err[i],
            color='#383838', alpha=.99, linewidth=1, capsize=4, capthick=1.05,
            zorder=3)

        axes2.bar(
            s * x[i],
            height=nPrime[i], yerr=0,
            color=colors[i], edgecolor='#383838',
            width=w, hatch='....', alpha=a, linewidth=.5,
            zorder=z)
        axes2.errorbar(
            x=s * x[i], y=nPrime[i], yerr=nPrime_err[i],
            color='#383838', alpha=.99, linewidth=1, capsize=4, capthick=1.05,
            zorder=3)

        axes3.bar(
            s * x[i] + w,
            height=sigmaZero[i], yerr=0,
            color=colors[i], edgecolor='#383838',
            width=w, hatch='', alpha=a, linewidth=.5,
            zorder=z)
        axes3.errorbar(
            x=s * x[i] + w, y=sigmaZero[i], yerr=sigmaZero_err[i],
            color='#383838', alpha=.99, linewidth=1, capsize=4, capthick=1.05,
            zorder=3)

        posList.append(s * x[i] - w), posList.append(s * x[i]), posList.append(s * x[i] + w)
        labelsList.append("$k'$"), labelsList.append("$n'$"), labelsList.append("$\sigma_0$'")

    axes.set_xticks(posList)
    axes.set_xticklabels(labelsList)

    return nPrime


def main(dataPath):
    fonts('C:/Users/petrus.kirsten/AppData/Local/Microsoft/Windows/Fonts/')
    plt.style.use('seaborn-v0_8-ticks')

    fig, axes = plt.subplots(
        figsize=(16, 7), ncols=2, nrows=1,
        gridspec_kw={'width_ratios': [1.25, 1]}, facecolor='snow')
    axFlow, axBars = axes[0], axes[1]

    fig.suptitle(f'Steps shear rate flow')
    xTitle, xLimits = ('Shear rate ($s^{-1}$)', (-15, 315))
    yTitle, yLimits = (f'Shear stress (Pa)', (0, 500))

    nSamples, colorSamples = getSamplesInfos(
        2, 3, 3,
        3, 1, 3,
        3, 4,
        'lightgray', 'hotpink', 'deepskyblue',
        'darkgray', 'crimson', 'mediumblue',
        'mediumorchid', 'rebeccapurple')
    data, labels = getSamplesData(dataPath, nSamples)

    dictData = {
        labels[0]: ([], [], []),
        labels[1]: ([], [], []),
        labels[2]: ([], [], []),
        labels[3]: ([], [], []),
        labels[4]: ([], [], []),
        labels[5]: ([], [], []),
        labels[6]: ([], [], []),
        labels[7]: ([], [], []),
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
    # plt.tight_layout()
    plt.show()

    fileName = f'0St_Car_CL-Flow'
    dirSave = Path(*Path(filePath[0]).parts[:Path(filePath[0]).parts.index('data') + 1])
    fig.savefig(f'{dirSave}' + f'\\{fileName}' + '.png', facecolor='w', dpi=600)

    # fitParams = pd.DataFrame(tableFitting)
    # fitParams.to_excel(f'{dirSave}' + f'\\{fileName}' + '.xlsx', index=False)

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

    main(dataPath=filePath)
