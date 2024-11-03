import numpy as np
import pandas as pd
from pathlib import Path

from matplotlib.patches import Rectangle
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt, ticker
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


def fitLinear(strain, slope, intercept):
    return slope * strain + intercept


def arraySplit(xArr, yArr, startValue, endValue):
    startIndex, endIndex = np.where(xArr >= startValue)[0][0], np.where(xArr <= endValue)[0][-1]

    return xArr[startIndex:endIndex], yArr[startIndex:endIndex]


def exportFit(
        sample,
        slope, slopeErr,
        tensile, tensileErr,
        rows
):
    dictData = {
        'Sample': sample,
        f'Slope (Pa)': slope, f'Slope (Pa) err': slopeErr,
        f'Tensile stress (Pa)': tensile, f'Tensile stress (Pa) err': tensileErr}

    rows.append(dictData)

    return rows


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
        number_samples
):
    def getSegments(dataframe):
        time = dataframe['t in s'].to_numpy()
        height = dataframe['h in mm'].to_numpy()
        force = dataframe['Fn in N'].to_numpy()

        seg_init, seg_end = (  # Identifying the job segments in the lists
            dataframe.index[dataframe['SegIndex'] == seg].to_list()[0] for seg in ['62|1', '62|98'])

        segmentsBreakage = lambda arr: (arr[seg_init:seg_end])  # Slice segments

        return {
            'time to break': segmentsBreakage(time) - segmentsBreakage(time)[0],
            'height to break': (1 - segmentsBreakage(height) / segmentsBreakage(height).max()) * 100,
            'force to break': segmentsBreakage(force)}

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
        segments = getSegments(df)
        samples[sample_type].append(segments)

    dict_data = {}  # Initialize dictionaries to hold the results

    for sample_type in samples:  # Populate dictionaries with consolidated sample data
        dict_data[f'{sample_type} time to break'] = [s['time to break'] for s in samples[sample_type]]
        dict_data[f'{sample_type} height to break'] = [s['height to break'] for s in samples[sample_type]]
        dict_data[f'{sample_type} force to break'] = [s['force to break'] for s in samples[sample_type]]

    return dict_data, sample_keys


def plotCompression(sampleName,
                    ax, x, y, yErr,
                    axTitle, yLabel, yLim, xLabel, xLim, axisColor,
                    curveColor, markerStyle, markerFColor, markerEColor, markerEWidth=0.5,
                    linearFitting=False, startVal=6, endVal=16, tableData=None):

    if sampleName == '0St/CL':
        startVal, endVal = 15, 25

    def legendLabel():
        legend = ax.legend(loc='upper left', fancybox=False, frameon=True, framealpha=0.9, fontsize=10)
        legend.get_frame().set_facecolor('w')
        legend.get_frame().set_edgecolor('lightsteelblue')
        legend.get_frame().set_linewidth(0.5)

    def configPlot():
        ax.set_title(axTitle, size=9, color='crimson')
        ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(0.75)

        ax.grid(True, which='major', axis='both', linestyle='-', linewidth=.75, color='lightgray', alpha=0.5)
        ax.grid(True, which='minor', axis='both', linestyle='--', linewidth=.5, color='lightgray', alpha=0.5)

        ax.set_xlabel(f'{xLabel}')
        ax.set_xlim(xLim)
        ax.xaxis.set_minor_locator(MultipleLocator(5 if not linearFitting else 1))
        ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}%"))

        ax.set_ylabel(f'{yLabel}', color=axisColor)
        ax.set_yscale('linear')
        ax.set_ylim(yLim)
        ax.tick_params(axis='y', colors=axisColor, which='both')
        ax.yaxis.set_major_locator(MultipleLocator(200 if not linearFitting else 100))
        ax.yaxis.set_minor_locator(MultipleLocator(50 if not linearFitting else 25))

    if linearFitting:
        x_toFit, y_toFit = arraySplit(x, y, startVal, endVal)
        (slope, intercept), covariance = curve_fit(fitLinear, x_toFit, y_toFit)  # p0=(y_mean[0], y_mean[-1], 100))
        (slopeErr, interceptErr) = np.sqrt(np.diag(covariance))
        tableData = exportFit(
            f'{sampleName}',
            slope * 100, slopeErr * 100,
            np.max(y), yErr[np.argmax(y)],
            tableData)

        xFit = np.linspace(startVal, endVal, 100)
        yFit = fitLinear(xFit, slope, intercept)
        ax.plot(  # plot fit curve
            xFit, yFit,
            color=markerFColor, alpha=0.8, lw=1, linestyle='--',
            label=f'Linear fitting', zorder=4)

        ax.errorbar(
            x, y, yErr,
            color=curveColor, alpha=.85,
            fmt='none', mfc=curveColor,
            capsize=2.5, capthick=1, linestyle='', lw=1,
            label=f'', zorder=2)

        ax.errorbar(
            x, y, 0,
            color=curveColor, alpha=.65,
            fmt='o', markersize=5.4,
            mfc=curveColor, mec='#383838', mew=.75,
            linestyle='',
            label=f'{sampleName}', zorder=3)
        # show text and rectangle at the linear region

        # DRAW DATA
        # textLabel, textCoord = (
        #     f'YM = ${slope * 100:.1f}$ $±$ ${slopeErr * 100:.1f}$ $Pa$', (xFit[-1] + 1, np.median(yFit) - 20))
        # textConfig = {'horizontalalignment': 'left', 'verticalalignment': 'top', 'color': 'k', 'size': 9}
        # ax.text(textCoord[0], textCoord[1], s=textLabel, **textConfig)

        # SHOW START AND END STRAIN VALUES
        # textLabel, textCoord = f'{startVal}%', (xFit[0] + 0.5, 3)
        # textConfig = {'horizontalalignment': 'left', 'verticalalignment': 'bottom', 'color': 'k', 'size': 8}
        # ax.text(textCoord[0], textCoord[1], s=textLabel, **textConfig, zorder=5)
        # textLabel, textCoord = f'{endVal}%', (xFit[-1] - 0.5, 3)
        # textConfig = {'horizontalalignment': 'right', 'verticalalignment': 'bottom', 'color': 'k', 'size': 8}
        # ax.text(textCoord[0], textCoord[1], s=textLabel, **textConfig, zorder=5)

        # DRAW RECTANGLE FOR LINEAR REGION
        # textLabel, textCoord = 'Linear elastic region', (xFit[-1] + 2, np.median(yFit))
        # textConfig = {'horizontalalignment': 'left', 'verticalalignment': 'top', 'color': 'crimson', 'size': 10}
        # rectConfig = [(xFit[0], 0), xFit[-1] - xFit[0], yFit[-1] + 2]
        # ax.text(textCoord[0], textCoord[1], s=textLabel, **textConfig)
        # rect = Rectangle(
        #     *rectConfig, linewidth=1, edgecolor='w', facecolor='crimson', alpha=0.1, zorder=1)
        # ax.add_patch(rect)

    else:
        ax.errorbar(
            x[::2], y[::2], yErr[::2],
            color=curveColor, alpha=.85,
            fmt='none', mfc=curveColor,
            capsize=2.5, capthick=1, linestyle='', lw=1,
            label=f'', zorder=2)

        ax.errorbar(
            x[::2], y[::2], 0,
            color=curveColor, alpha=.65,
            fmt='o', markersize=5.4,
            mfc=curveColor, mec='#383838', mew=.75,
            linestyle='',
            label=f'{sampleName}', zorder=3)
        legendLabel()

    configPlot()

    return tableData if linearFitting else None


def plotBars(axes, data, keys, colors, h, z, a=.9):
    def configPlot(ax, yTitle, yLim):
        if yTitle == 'Young modulus (Pa)':
            ax.grid(which='major', axis='y', linestyle='-', linewidth=.75, color='lightgray', alpha=0.5, zorder=-1)
            ax.grid(which='minor', axis='y', linestyle='--', linewidth=.5, color='lightgray', alpha=0.5, zorder=-1)

        ax.tick_params(axis='x', labelsize=10, length=4)
        ax.tick_params(axis='y', which='both', labelsize=9, pad=1, length=0)

        ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(.75)
        ax.spines[['top', 'bottom', 'left', 'right']].set_color('#303030')
        ax.set_xticks([])
        ax.set_xlim([-2, 13.5])

        ax.set_ylabel(yTitle)
        ax.set_ylim(yLim)

        ax.yaxis.set_major_locator(MultipleLocator(yLim[-1] / 5))
        ax.yaxis.set_minor_locator(MultipleLocator(yLim[-1] / 20))

    axes2 = axes.twinx()

    lim1, lim2 = (0, 2500), (0, 750)

    configPlot(axes, "Young modulus (Pa)", lim1)
    configPlot(axes2, "Stress peak (Pa)", lim2)

    posList, labelsList = [], []

    bin_width, space_samples = 1, 3
    samples = keys
    x = np.arange(space_samples * len(samples))
    for i in range(len(samples)):
        slope, slope_err = data[i]['Slope (Pa)'], data[i]['Slope (Pa) err']
        peak, peak_err = data[i]['Tensile stress (Pa)'], data[i]['Tensile stress (Pa) err']
        
        axes.bar(
            space_samples * x[i] - bin_width,
            height=slope, yerr=0,
            color=colors[i], edgecolor='#383838',
            width=bin_width, hatch='///', alpha=a, linewidth=.5,
            zorder=z)
        axes.errorbar(
            x=space_samples * x[i] - bin_width, y=slope, yerr=slope_err,
            color='#383838', alpha=.99, linewidth=1, capsize=5, capthick=1.05,
            zorder=3)
        axes.text(
            space_samples * x[i] - bin_width - .15,
            slope + slope_err + lim1[1]*.075,
            f'{slope:.{1}f} ± {slope_err:.{1}f}',
            va='center', ha='left', rotation=90,
            color='#383838', fontsize=9)

        axes2.bar(
            space_samples * x[i],
            height=peak, yerr=0,
            color=colors[i], edgecolor='#383838',
            width=bin_width, hatch='', alpha=a, linewidth=.5,
            zorder=z)
        axes2.errorbar(
            x=space_samples * x[i], y=peak, yerr=peak_err,
            color='#383838', alpha=.99, linewidth=1, capsize=5, capthick=1.05,
            zorder=3)
        axes2.text(
            space_samples * x[i] - .15,
            peak + peak_err + lim2[1]*.075,
            f'{peak:.{1}f} ± {peak_err:.{1}f}',
            va='center', ha='left', rotation=90,
            color='#383838', fontsize=9)

        posList.append(space_samples * x[i] - bin_width)
        posList.append(space_samples * x[i])
        labelsList.append("YM"), labelsList.append("Peak")

    axes.set_xticks(posList)
    axes.set_xticklabels(labelsList, rotation=45)


def main(dataPath, fileName):
    fonts('C:/Users/petrus.kirsten/AppData/Local/Microsoft/Windows/Fonts/')
    plt.style.use('seaborn-v0_8-ticks')

    fig, axes = plt.subplots(
        figsize=(18, 8), ncols=3,
        gridspec_kw={'width_ratios': [3, 2, 2]}, facecolor='snow')

    axComplete, axFit, axBars = axes[0], axes[1], axes[2]
    fig.suptitle(f'Compression modulus')

    nSamples, colorSamples = getSamplesInfos(
        3, 2, 5,
        4, 2,
        'grey', 'mediumvioletred', 'royalblue',
        '#fb7e8f', '#e30057')

    raw_data, labels = getSamplesData(dataPath, nSamples)

    data = {
        labels[0]: ([], []),
        labels[1]: ([], []),
        labels[2]: ([], []),
        labels[3]: ([], []),
        labels[4]: ([], []),
    }

    for k, (x, y) in data.items():
        x.append(raw_data[f'{k} height to break'])
        y.append(raw_data[f'{k} force to break'])

    dataFitting = []
    for k, c in zip(data, colorSamples):
        strain, stress = (np.mean(data[k][0], axis=1)[0]), np.abs(np.mean(data[k][1], axis=1)[0] / (35 / 1000) + 5)
        stressErr = np.abs(np.std(data[k][0], axis=1)[0] / (35 / 1000))

        plotCompression(
            sampleName=f'{k}',
            ax=axComplete, axisColor='k',
            x=strain, y=stress, yErr=stressErr,
            axTitle='', yLabel='Stress (Pa)', yLim=(0, 1100), xLabel='Strain', xLim=(0, 100),
            curveColor=c, markerStyle='o', markerFColor=c, markerEColor='k',
            linearFitting=False)

        dataFitting = plotCompression(
            sampleName=f'{k}',
            ax=axFit, axisColor='k',
            x=strain[:35], y=stress[:35], yErr=stressErr[:35],
            axTitle='', yLabel='Stress (Pa)', yLim=(0, 425), xLabel='Strain', xLim=(0, 30),
            curveColor=c, markerStyle='o', markerFColor=c, markerEColor='k',
            linearFitting=True, tableData=dataFitting)

    plotBars(
        axBars, dataFitting, labels,
        colorSamples, h='', z=2)

    plt.subplots_adjust(
        wspace=0.15, hspace=0.1,
        top=0.940, bottom=0.08,
        left=0.05, right=0.96)
    plt.show()

    dirSave = Path(*Path(filesPath[0]).parts[:Path(filesPath[0]).parts.index('data') + 1])
    fig.savefig(f'{dirSave}' + f'\\{fileName}' + '.png', facecolor='w', dpi=600)

    print(f'\n\n· Chart saved at\n{dirSave}.')


if __name__ == '__main__':
    # folderPath = "C:/Users/petrus.kirsten/PycharmProjects/RheometerPlots/data"
    folderPath = "C:/Users/Petrus Kirsten/Documents/GitHub/RheometerPlots/data"

    filesPath = [
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
        # folderPath + "/231024/kC_CL/kC_CL-compression-1.xlsx",   # 299
        # folderPath + "/231024/kC_CL/kC_CL-compression-1b.xlsx",  # 299
        folderPath + "/231024/kC_CL/kC_CL-compression-3.xlsx",  # 99
        folderPath + "/231024/kC_CL/kC_CL-compression-4.xlsx",  # 99
    ]

    main(filesPath, '0St_Car-CompressionToBreakage')
