import numpy as np
import pandas as pd
from pathlib import Path

from matplotlib.ticker import MultipleLocator
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


def exportFit(
        sample,
        data, err,
        rows):
    keys = ('$\\tau_0$', '$\\tau_e$', '$\lambda$')
    values = (data, err)

    dictData = {'Sample': sample}
    iParams = 0
    for key, value in zip(keys, range(len(values[0]))):
        dictData[f'{key}'] = values[0][iParams]
        dictData[f'± {key}'] = values[1][iParams]
        iParams += 1

    rows.append(dictData)

    return rows


def funcTransient(t, tau_0, tau_e, time_cte):
    return tau_e + (tau_0 - tau_e) * np.exp(- t / time_cte)


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
        """
        Extracts time, shear rate, shear stress, and viscosity segments from the dataframe.
        Returns tuples of constant and step segments.
        """
        time = dataframe['t in s'].to_numpy()
        shear_rate = dataframe['ɣ̇ in 1/s'].to_numpy()
        shear_stress = dataframe['τ in Pa'].to_numpy()
        viscosity = dataframe['η in mPas'].to_numpy()

        # Identifying segments in the data
        seg3, seg4 = (dataframe.index[dataframe['SegIndex'] == seg].to_list()[0] for seg in ['3|1', '4|1'])

        # Slice segments
        segments = lambda arr: (arr[seg3:seg4])  # Returns (constant segment, step segment)
        t_cte = segments(time)

        return {
            'time': [t_cte - t_cte[0]],
            'shear_stress': segments(shear_stress),
            'viscosity': segments(viscosity)
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
            [sample_keys[7]] * number_samples[7])

    for sample_type, path in zip(sample_labels, dataPath):
        df = pd.read_excel(path)
        segments = getSegments(df)
        samples[sample_type].append(segments)

    dict_cteRate = {}

    for sample_type in samples:
        dict_cteRate[f'{sample_type} time'] = [s['time'] for s in samples[sample_type]]
        dict_cteRate[f'{sample_type} stress'] = [s['shear_stress'] for s in samples[sample_type]]
        dict_cteRate[f'{sample_type} viscosity'] = [s['viscosity'] for s in samples[sample_type]]

    return dict_cteRate, list(samples.keys())


def downsampler(array, n=200):
    if len(array) > n:
        step = len(array) // n  # Calculate step size
        return array[::step][:n]
    return array


def plotFlow(listRows, sampleName,
             ax, x, y, yErr,
             axTitle, yLabel, yLim, xLabel, xLim,
             curveColor, markerStyle,
             fit='', logScale=False):
    def legendLabel():
        legend = ax.legend(
            fancybox=False, frameon=True,
            framealpha=0.9, fontsize=9, ncols=2,
            loc='upper left'
        )
        legend.get_frame().set_facecolor('w')
        legend.get_frame().set_edgecolor('whitesmoke')

    def configPlot():
        ax.set_title(axTitle, size=9, color='crimson')
        ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(0.75)

        ax.set_xlabel(f'{xLabel}')
        ax.set_xscale('log' if logScale else 'linear')
        ax.set_xlim(xLim)
        ax.xaxis.set_major_locator(MultipleLocator(50))
        ax.xaxis.set_minor_locator(MultipleLocator(10))

        ax.set_ylabel(f'{yLabel}')
        ax.set_yscale('log' if logScale else 'linear')
        ax.set_ylim(yLim)
        ax.yaxis.set_major_locator(MultipleLocator(100))
        ax.yaxis.set_minor_locator(MultipleLocator(25))

    params, covariance = curve_fit(
        funcTransient, x, y)
    # p0=(x[0], y[-1], 100))  # method='trf')  # method='dogbox', maxfev=5000)
    errors = np.sqrt(np.diag(covariance))

    tau_0, tau_e, t_cte = params
    x_fit = np.linspace(0, 300, 300)
    y_fit = funcTransient(x_fit, tau_0, tau_e, t_cte)

    listRows = exportFit(
        f'{sampleName}',
        params, errors,
        listRows)

    ax.errorbar(
        x, y, yerr=yErr, color=curveColor, alpha=.65,
        fmt=markerStyle, markersize=7, mec='k', mew=0.5,
        capsize=3, lw=.5, linestyle='',  # ecolor='k'
        label=f'{sampleName}', zorder=3)

    ax.plot(
        x_fit, y_fit, color=curveColor, linestyle='-.', linewidth=.75,
        zorder=2)

    configPlot()
    legendLabel()

    return listRows


def plotBars(title, axes, data, colors, a, z):
    def configPlot(ax, yTitle, yLim):
        ax.set_title(title, size=10, color='k')
        if yTitle == "$\\tau_0$ (Pa)":
            ax.grid(which='major', axis='y', linestyle='-', linewidth=1, color='lightgray', alpha=0.5, zorder=-1)
            ax.grid(which='minor', axis='y', linestyle='--', linewidth=.75, color='lightgray', alpha=0.5, zorder=-1)

        if yTitle == "$\\tau_e$ (Pa)":
            ax.yaxis.tick_left()
            ax.yaxis.set_label_position('left')
            ax.tick_params(axis='y', which='both', direction='in', pad=-28)
            ax.yaxis.set_label_coords(.125, .5)

        ax.tick_params(axis='x', labelsize=10, length=4)
        # ax.tick_params(axis='y', which='both', direction='out', pad=1)

        ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(.75)
        ax.spines[['top', 'bottom', 'left', 'right']].set_color('#303030')

        ax.set_xticks([])
        ax.set_xlim([-3.5, 23])

        ax.set_ylabel(yTitle)
        ax.set_ylim(yLim)

        ax.yaxis.set_major_locator(MultipleLocator(yLim[1] / 10))
        ax.yaxis.set_minor_locator(MultipleLocator(yLim[1] / 50))

    axes2, axes3 = axes.twinx(), axes.twinx()

    configPlot(axes, "$\\tau_0$ (Pa)", (0, 1600))
    configPlot(axes2, "$\\tau_e$ (Pa)", (0, 500))
    configPlot(axes3, "$\lambda$ (s)", (0, 55))

    samples = [d['Sample'] for d in data]
    kPrime, kPrime_err = [d["$\\tau_0$"] for d in data], [d["± $\\tau_0$"] for d in data]
    nPrime, nPrime_err = [d["$\\tau_e$"] for d in data], [d["± $\\tau_e$"] for d in data]
    sigmaZero, sigmaZero_err = [d["$\lambda$"] for d in data], [d["± $\lambda$"] for d in data]

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
        labelsList.append("$\\tau_0$"), labelsList.append("$\\tau_e$"), labelsList.append("$\lambda$")

    axes.set_xticks(posList)
    axes.set_xticklabels(labelsList)

    return nPrime


def main(dataPath):
    fonts('C:/Users/petrus.kirsten/AppData/Local/Microsoft/Windows/Fonts/')
    plt.style.use('seaborn-v0_8-ticks')

    fig, axes = plt.subplots(
        figsize=(16, 7), ncols=2, nrows=1,
        gridspec_kw={'width_ratios': [1.25, 1]}, facecolor='snow')
    axStress, axBars = axes[0], axes[1]

    fig.suptitle(f'Constant shear rate flow')
    xTitle, xLimits = (f'Time (s)', (0, 200))
    yTitle, yLimits = (f'Shear stressMean (Pa)', (0, 800))
    yTitleVisc, yLimitsVisc = f'Viscosity (mPa·s)', (yLimits[0] * 3.33, yLimits[1] * 3.33)

    axesVisc = axStress.twinx()
    axesVisc.spines[['top', 'bottom', 'left', 'right']].set_linewidth(0)
    axesVisc.set_ylabel(f'{yTitleVisc}')
    axesVisc.set_ylim(yLimitsVisc)
    axesVisc.yaxis.set_major_locator(MultipleLocator(200))
    axesVisc.yaxis.set_minor_locator(MultipleLocator(50))

    nSamples, colorSamples = getSamplesInfos(
        2, 3, 3,
        3, 1, 3,
        3, 4,
        'lightgray', 'hotpink', 'deepskyblue',
        'darkgray', 'crimson', 'mediumblue',
        'mediumorchid', 'rebeccapurple')
    data, labels = getSamplesData(dataPath, nSamples)

    fitModeStress, fitModeVisc = 'transient', ''

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

    tableStress = []

    for key, (t, s, v) in dictData.items():
        t.append(data[f'{key} time'])
        s.append(data[f'{key} stress'])
        v.append(data[f'{key} viscosity'])

    for key, color in zip(dictData, colorSamples):
        time, stress, stressErr = (
            dictData[key][0][0],
            dictData[key][1][0],
            dictData[key][1][0])

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

        tableStress = plotFlow(
            listRows=tableStress, ax=axStress,
            x=timeMean, y=stressMean, yErr=stressErrMean,
            axTitle='', yLabel=yTitle, yLim=yLimits, xLabel=xTitle, xLim=xLimits,
            curveColor=color, markerStyle='o',
            sampleName=f'{key}', fit=fitModeStress)

    _ = plotBars(
        '',
        axBars, tableStress,
        colorSamples, a=.85, z=2)

    plt.subplots_adjust(
        hspace=0, wspace=0.30,
        top=0.940, bottom=0.095,
        left=0.090, right=0.900)
    plt.show()

    fileName = 'St_Car_CL-Thixotropy'
    dirSave = Path(*Path(filePath[0]).parts[:Path(filePath[0]).parts.index('data') + 1])
    fig.savefig(f'{dirSave}' + f'\\{fileName}' + '.png', facecolor='w', dpi=600)

    # fitParams = pd.DataFrame(tableStress)
    # fitParams.to_excel(f'{dirSave}' + f'\\{fileName}' + '.xlsx', index=False)

    print(f'\n\n· Chart and tableStress with fitted parameters saved at\n{dirSave}.')


if __name__ == '__main__':
    # folderPath = "C:/Users/petrus.kirsten/PycharmProjects/RheometerPlots/data"
    folderPath = "C:/Users/Petrus Kirsten/Documents/GitHub/RheometerPlots/data"
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
