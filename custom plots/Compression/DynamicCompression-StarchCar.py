import numpy as np
import pandas as pd
from pathlib import Path

from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt, ticker, gridspec
from matplotlib.ticker import MultipleLocator
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


def funcRelaxation(t, sigma_0, sigma_e, time_cte):
    return sigma_e + (sigma_0 - sigma_e) * np.exp(- t / time_cte)


def getSamplesInfos(
        # quantites
        stCL_n, st_kcCL_n, st_icCL_n,
        kc_n, kcCL_n,
        # colors
        stCL_color, st_kcCL_color, st_icCL_color,
        kc_color, kcCL_color
):
    number_samples = [
        stCL_n, st_kcCL_n, st_icCL_n,
        kc_n, kcCL_n]

    colors_samples = [
        stCL_color, st_kcCL_color, st_icCL_color,
        kc_color, kcCL_color]

    return number_samples, colors_samples


def getSamplesData(
        dataPath,
        number_samples):
    def downsampler(array, n=1121):
        if len(array) > n:
            step = len(array) // n  # Calculate step size
            return array[::step][:n]
        return array

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

    samples = {
        '0St/CL': [], '0St + kCar/CL': [], '0St + iCar/CL': [],
        'kCar': [], 'kCar/CL': []
    }
    sample_keys = list(samples.keys())
    sample_labels = (
            [sample_keys[0]] * number_samples[0] +
            [sample_keys[1]] * number_samples[1] +
            [sample_keys[2]] * number_samples[2] +
            [sample_keys[3]] * number_samples[3] +
            [sample_keys[4]] * number_samples[4]
    )

    for sample_type, path in zip(sample_labels, dataPath):  # Read data and categorize based on sample type
        df = pd.read_excel(path)
        segments = getSegments(df, '2|1' if '171024' in path or not 'kC-compression-4' in path else '1|1', '62|1')
        samples[sample_type].append(segments)

    dict_data = {}

    for sample_type in samples:
        # Apply the downsampling to each data list
        dict_data[f'{sample_type} time'] = [downsampler(s['time']) for s in samples[sample_type]]
        dict_data[f'{sample_type} height'] = [downsampler(s['height']) for s in samples[sample_type]]
        dict_data[f'{sample_type} force'] = [downsampler(s['force']) for s in samples[sample_type]]

    return dict_data, sample_keys


def plotCompression(sampleName,
                    ax, x, y, yErr,
                    axTitle, yLabel, yLim, xLabel, xLim, axisColor,
                    curveColor, markerStyle, markerFColor, markerEColor, markerEWidth=0.5,
                    strain=False, lineStyle='', logScale=False):
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
        interp_yErr_lower, interp_yErr_upper = interp1d(x, y - yErr, kind='cubic'), interp1d(x, y + yErr, kind='cubic')
        yErr_lower_smooth, yErr_upper_smooth = interp_yErr_upper(x_smooth), interp_yErr_lower(x_smooth)

        ax.fill_between(
            x_smooth, yErr_lower_smooth, yErr_upper_smooth,
            color=curveColor, alpha=0.15, lw=0,
            zorder=3)
        ax.errorbar(
            x, y, 0,
            color=curveColor, alpha=.35,
            fmt=markerStyle, markersize=3.5, mfc=markerFColor, mec=markerEColor, mew=markerEWidth,
            capsize=0, lw=.5, linestyle=lineStyle,
            label=f'{sampleName}',
            zorder=4)
        legendLabel()


def plotCycles(sampleName,
               ax, y, yErr,
               axTitle, yLabel, yLim, xLabel, xLim, axisColor,
               curveColor, markerFColor, markerEColor, paramsList, markerEWidth=0.5,
               lineStyle='', logScale=False):
    def legendLabel():
        legend = ax.legend(fancybox=False, frameon=True, framealpha=0.9, fontsize=9)
        legend.get_frame().set_facecolor('w')
        legend.get_frame().set_edgecolor('whitesmoke')

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

    target_size = 29 * (1121 // 29)
    y, yErr = y[:target_size], yErr[:target_size]
    y, yErr = y.reshape(29, -1), yErr.reshape(29, -1)

    listYmax, listYmaxErr = [], []
    listYmin, listYminErr = [], []

    for cycle in range(len(y)):
        yMax, yMin = y[cycle][np.argmax(y[cycle])], y[cycle][np.argmin(y[cycle])]
        yMaxErr, yMinErr = yErr[cycle][np.argmax(y[cycle])], yErr[cycle][np.argmin(y[cycle])]

        listYmax.append(yMax), listYmaxErr.append(yMaxErr)
        listYmin.append(yMin), listYminErr.append(yMinErr)

    x_values = np.arange(1, len(listYmax) + 1)
    yMax_values, yMin_values = np.array(listYmax) / (35 / 1000) + 7.5, np.array(listYmin) / (35 / 1000) + 7.5
    yMax_errors, yMin_errors = np.array(listYmaxErr) / (35 / 1000), np.array(listYminErr) / (35 / 1000)

    # Fitting
    params, covariance = curve_fit(
        funcRelaxation, x_values, yMax_values,
        p0=(yMax_values[0], yMax_values[-1], 5))  # method='trf')  # method='dogbox', maxfev=5000)
    errors = np.sqrt(np.diag(covariance))

    stress_0, stress_e, t_cte = params
    x_fit = np.linspace(0, 30, 120)
    y_fit = funcRelaxation(x_fit, stress_0, stress_e, t_cte)

    paramsList.append([params.tolist(), errors.tolist()])

    configPlot()

    ax.errorbar(
        x_values, yMax_values, yMax_errors,
        color=curveColor, alpha=.85,
        fmt='none', mfc=curveColor,
        capsize=2.5, capthick=1, linestyle='', lw=1,
        label=f'', zorder=2)

    ax.errorbar(
        x_values, yMax_values, 0,
        color=curveColor, alpha=.65,
        fmt='^', markersize=6.5,
        mfc=curveColor, mec='#383838', mew=.75,
        linestyle='',
        label=f'{sampleName}', zorder=3)

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


def plotBars(title, axes, data, keys, colors, h, z, a=.9):
    def configPlot(ax, yTitle, yLim):
        if yTitle == 'Initial stress peak / Equilibrim stress peak (Pa)':
            ax.grid(which='major', axis='y', linestyle='-', linewidth=.75, color='lightgray', alpha=0.5, zorder=-1)
            ax.grid(which='minor', axis='y', linestyle='--', linewidth=.5, color='lightgray', alpha=0.5, zorder=-1)

        ax.tick_params(axis='x', labelsize=10, length=4)
        ax.tick_params(axis='y', which='both', labelsize=9, pad=1, length=0)

        ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(.75)
        ax.spines[['top', 'bottom', 'left', 'right']].set_color('#303030')
        ax.set_xticks([])
        ax.set_xlim([-2, 12.5])

        ax.set_ylabel(yTitle)
        ax.set_ylim(yLim)

        ax.yaxis.set_major_locator(MultipleLocator(yLim[-1] / 5))
        ax.yaxis.set_minor_locator(MultipleLocator(yLim[-1] / 20))

    axes2 = axes.twinx()

    lim1, lim2 = (0, 250), (0, 15)

    configPlot(axes, "Initial stress peak / Equilibrim stress peak (Pa)", lim1)
    configPlot(axes2, "Time decay constant (1/cycle)", lim2)

    posList, labelsList = [], []

    bin_width, space_samples = 1, 3.5
    samples = keys
    x = np.arange(space_samples * len(samples))

    for i in range(len(samples)-1):
        initialStress, initialStress_err = data[i][0][0], data[i][1][0]
        equilibStress, equilibStress_err = data[i][0][1], data[i][1][1]
        timeCte, timeCte_err = data[i][0][2], data[i][1][2]
        axes.bar(
            space_samples * x[i] - bin_width,
            height=initialStress, yerr=0,
            color=colors[i], edgecolor='#383838',
            width=bin_width, hatch='\\\\\\', alpha=a, linewidth=.5,
            zorder=z)
        axes.errorbar(
            x=space_samples * x[i] - bin_width, y=initialStress, yerr=initialStress_err,
            color='#383838', alpha=.99, linewidth=1, capsize=5, capthick=1.05,
            zorder=3)
        axes.text(
            space_samples * x[i] - bin_width - .1,
            initialStress + initialStress_err + lim1[1]*.075,
            f'{initialStress:.{1}f} ± {initialStress_err:.{1}f}',
            va='center', ha='left', rotation=90,
            color='#383838', fontsize=9)

        axes.bar(
            space_samples * x[i],
            height=equilibStress, yerr=0,
            color=colors[i], edgecolor='#383838',
            width=bin_width, hatch='....', alpha=a, linewidth=.5,
            zorder=z)
        axes.errorbar(
            x=space_samples * x[i], y=equilibStress, yerr=equilibStress_err,
            color='#383838', alpha=.99, linewidth=1, capsize=5, capthick=1.05,
            zorder=3)
        axes.text(
            space_samples * x[i] - .1,
            equilibStress + equilibStress_err + lim1[1]*.075,
            f'{equilibStress:.{1}f} ± {equilibStress_err:.{1}f}',
            va='center', ha='left', rotation=90,
            color='#383838', fontsize=9)

        axes2.bar(
            space_samples * x[i] + bin_width,
            height=timeCte, yerr=0,
            color=colors[i], edgecolor='#383838',
            width=bin_width, hatch=h, alpha=a, linewidth=.5,
            zorder=z)
        axes2.errorbar(
            x=space_samples * x[i] + bin_width, y=timeCte, yerr=timeCte_err,
            color='#383838', alpha=.99, linewidth=1, capsize=5, capthick=1.05,
            zorder=3)
        axes2.text(
            space_samples * x[i] + bin_width - .1,
            timeCte + timeCte_err + lim2[1]*.075,
            f'{timeCte:.{1}f} ± {timeCte_err:.{1}f}',
            va='center', ha='left', rotation=90,
            color='#383838', fontsize=9)

        if i != 4:
            posList.append(space_samples * x[i] - bin_width), posList.append(space_samples * x[i]), posList.append(space_samples * x[i] + bin_width)
            labelsList.append('Initial'), labelsList.append('Eq'), labelsList.append('$\\tau$')

    axes.set_xticks(posList)
    axes.set_xticklabels(labelsList, rotation=45)


def main(dataPath):
    fonts('C:/Users/petrus.kirsten/AppData/Local/Microsoft/Windows/Fonts/')
    plt.style.use('seaborn-v0_8-ticks')

    fig, gs = (plt.figure(figsize=(12, 10), facecolor='snow'),
               gridspec.GridSpec(2, 2, height_ratios=[1, 1.5], width_ratios=[1, 1]))

    axStress, axCycles, axBars = fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])

    fig.suptitle(f'Oscilatory compression')
    s1Title, s1Limits = f'Stress (Pa)', (0, 210)
    s2Title, s2Limits = f'Stress peak (Pa)', (0, 210)
    x1Title, x1Limits = f'Time (s)', (0, 75)
    x2Title, x2Limits = f'Cycle', (0, 30)

    nSamples, colorSamples = getSamplesInfos(
        3, 2, 5,
        4, 4,
        'grey', 'mediumvioletred', 'royalblue',
        '#fb7e8f', '#e30057')

    data, labels = getSamplesData(dataPath, nSamples)

    dataList = {
        labels[0]: ([], [], []),
        labels[1]: ([], [], []),
        labels[2]: ([], [], []),
        labels[3]: ([], [], []),
        labels[4]: ([], [], [])
    }

    means_hMax, fitParams = [], []

    for key, (x, height, f) in dataList.items():
        x.append(data[f'{key} time'])
        height.append(data[f'{key} height'])
        f.append(data[f'{key} force'])

    for key, c in zip(dataList, colorSamples):
        time, height, stress = dataList[key][0][0], dataList[key][1][0], dataList[key][2][0]
        hMax_mean, stressMean, stressErrMean = np.max(height, axis=1), np.mean(stress, axis=0), np.std(stress, axis=0)
        means_hMax.append(hMax_mean.tolist())

        plotCompression(
            ax=axStress, axisColor='#303030',
            x=np.mean(time, axis=0), y=stressMean / (35 / 1000) + 7.5, yErr=stressErrMean / (35 / 1000),
            axTitle='', yLabel=s1Title, yLim=s1Limits, xLabel=x1Title, xLim=x1Limits,
            curveColor=c, markerStyle='o', markerFColor=c, markerEColor=c,
            sampleName=f'{key}')

        fitParams = plotCycles(
            ax=axCycles, axisColor='#303030',
            y=stressMean, yErr=stressErrMean,
            axTitle='', yLabel=s2Title, yLim=s2Limits, xLabel=x2Title, xLim=x2Limits,
            curveColor=c, markerFColor=c, markerEColor=c,
            paramsList=fitParams, sampleName=f'{key}')

    plotBars(
        '',
        axBars, fitParams, labels,
        colorSamples, h='', z=2)

    plt.subplots_adjust(
        hspace=0.2, wspace=0.2,
        top=0.94, bottom=0.09,
        left=0.06, right=0.95)
    plt.show()

    fileName = '0St_Car-DynamicCompression'
    dirSave = Path(*Path(filePath[0]).parts[:Path(filePath[0]).parts.index('data') + 1])
    fig.savefig(f'{dirSave}' + f'\\{fileName}' + '.png', facecolor='w', dpi=600)
    print(f'\n\n· Chart saved at\n{dirSave}.')


if __name__ == '__main__':
    # folderPath = "C:/Users/petrus.kirsten/PycharmProjects/RheometerPlots/data"
    folderPath = "C:/Users/Petrus Kirsten/Documents/GitHub/RheometerPlots/data"
    # filePath = [
    #     folderPath + "/old/200924/7PSt_2_Compression.xlsx"
    # ]

    filePath = [

        # 0St/CL
        folderPath + "/171024/10_0St_CL/10_0St_CL-compression-1.xlsx",
        folderPath + "/171024/10_0St_CL/10_0St_CL-compression-2.xlsx",
        folderPath + "/171024/10_0St_CL/10_0St_CL-compression-3.xlsx",

        # 0St + kCar/CL
        folderPath + "/171024/10_0St_kC_CL/10_0St_kC_CL-compression-1.xlsx",
        folderPath + "/171024/10_0St_kC_CL/10_0St_kC_CL-compression-2.xlsx",

        # 0St + iCar/CL - TODO: identify the outilers samples
        folderPath + "/171024/10_0St_iC_CL/10_0St_iC_CL-compression-1.xlsx",
        folderPath + "/171024/10_0St_iC_CL/10_0St_iC_CL-compression-2b.xlsx",
        folderPath + "/171024/10_0St_iC_CL/10_0St_iC_CL-compression-3.xlsx",
        folderPath + "/171024/10_0St_iC_CL/10_0St_iC_CL-compression-4.xlsx",
        folderPath + "/171024/10_0St_iC_CL/10_0St_iC_CL-compression-5.xlsx",

        # kC
        folderPath + "/231024/kC/kC-compression-1.xlsx",
        folderPath + "/231024/kC/kC-compression-2.xlsx",
        folderPath + "/231024/kC/kC-compression-3.xlsx",
        folderPath + "/231024/kC/kC-compression-4.xlsx",

        # kC/CL - TODO: adjust the 299 points
        folderPath + "/231024/kC_CL/kC_CL-compression-1.xlsx",  # 299
        folderPath + "/231024/kC_CL/kC_CL-compression-1b.xlsx",  # 299
        folderPath + "/231024/kC_CL/kC_CL-compression-3.xlsx",  # 99
        folderPath + "/231024/kC_CL/kC_CL-compression-4.xlsx",  # 99
    ]

    main(dataPath=filePath)
