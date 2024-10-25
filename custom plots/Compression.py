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
        data, err,
        rows
):
    dictData = {'Sample': sample, f'Slope (Pa)': data, f'Slope (Pa) err': err}
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
    def legendLabel():
        legend = ax.legend(loc='upper left', fancybox=False, frameon=True, framealpha=0.9, fontsize=10)
        legend.get_frame().set_facecolor('w')
        legend.get_frame().set_edgecolor('lightsteelblue')
        legend.get_frame().set_linewidth(0.5)

    def configPlot():
        ax.set_title(axTitle, size=9, color='crimson')
        ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(0.75)

        ax.grid(True, which='both', axis='both', linestyle='--', linewidth=0.5, color='lightsteelblue', alpha=0.5)

        ax.set_xlabel(f'{xLabel}')
        ax.set_xlim(xLim)
        ax.xaxis.set_minor_locator(MultipleLocator(5))
        ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}%"))

        ax.set_ylabel(f'{yLabel}', color=axisColor)
        ax.set_yscale('linear')
        ax.set_ylim(yLim)
        ax.tick_params(axis='y', colors=axisColor, which='both')
        ax.yaxis.set_major_locator(MultipleLocator(200 if not linearFitting else 100))
        ax.yaxis.set_minor_locator(MultipleLocator(50 if not linearFitting else 25))
        if linearFitting:
            ax.tick_params(axis='y', which='both', direction='in', pad=-28)
            y_tick_labels = ["" if tick == 0 else f"{int(tick)}" for tick in ax.get_yticks()]
            ax.set_yticklabels(y_tick_labels)

    if linearFitting:
        # configs to split the values at the linear elastic region - TODO: adjust linear region for 0St/CL
        x_toFit, y_toFit = arraySplit(x, y, startVal, endVal)
        (slope, intercept), covariance = curve_fit(fitLinear, x_toFit, y_toFit)  # p0=(y_mean[0], y_mean[-1], 100))
        (slopeErr, interceptErr) = np.sqrt(np.diag(covariance))
        tableData = exportFit(
            f'{sampleName}',
            slope * 100, slopeErr * 100,
            tableData)

        xFit = np.linspace(startVal, endVal, 100)
        yFit = fitLinear(xFit, slope, intercept)
        ax.plot(  # plot fit curve
            xFit, yFit,
            color=markerFColor, alpha=0.8, lw=.75, linestyle='-',
            label=f'Linear fitting', zorder=4)

        ax.errorbar(  # plot raw data
            x, np.abs(y) + 5, yErr,
            color=curveColor, alpha=0.65,
            fmt=markerStyle, markersize=4.5, mfc=markerFColor, mec=markerFColor, mew=markerEWidth,
            capsize=1.5, lw=1, linestyle='',
            label=f'',
            zorder=3)

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
            x, y, yErr,
            color=curveColor, alpha=0.65,
            fmt=markerStyle, markersize=4.5, mfc=markerFColor, mec=markerFColor, mew=markerEWidth,
            capsize=1.5, lw=1, linestyle='',
            label=f'{sampleName}',
            zorder=3)
        legendLabel()

    configPlot()

    return tableData if linearFitting else None


def plotBars(ax, title, data, colors):
    ax.set_title(title, size=9, color='#303030')
    # Extracting data for plotting
    samples = [d['Sample'] for d in data]
    slopes = [d['Slope (Pa)'] for d in data]
    slope_errs = [d['Slope (Pa) err'] for d in data]

    ax.tick_params(axis='both', labelsize=8, length=0)
    ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(0.75)
    ax.spines[['top', 'bottom', 'left', 'right']].set_color('dimgrey')

    ax.set_xticks([])
    # ax.set_xticklabels(samples)
    # ax_inset.yaxis.tick_right()
    # ax_inset.yaxis.set_label_position('right')
    ax.set_yticks([])
    ax.set_ylim(0, 2300)
    x = np.arange(len(data))

    for i in range(len(slopes)):
        ax.bar(x[i], height=slopes[i], yerr=0, width=0.4,
               color=colors[i], edgecolor='#303030', alpha=0.75, linewidth=0.75)

        ax.errorbar(x=x[i], y=slopes[i], yerr=slope_errs[i], alpha=.85,
                    color='#303030', linestyle='', capsize=4, linewidth=0.75)
        # TODO: add tensile stress
        ax.text(x[i], slopes[i] + slope_errs[i] + 50,
                f'{slopes[i]:.0f} ± {slope_errs[i]:.0f} Pa',
                size=8, ha='center', va='bottom', color='black')
    # ax.text(
    #     0.5, 1.1, "Average G' values (Pa)",
    #     ha='center', va='top', fontsize=9, transform=ax.transAxes)
    # ax.set_facecolor('snow')


def main(dataPath, fileName):
    fonts('C:/Users/petrus.kirsten/AppData/Local/Microsoft/Windows/Fonts/')
    plt.style.use('seaborn-v0_8-ticks')

    fig, axes = plt.subplots(figsize=(18, 8), ncols=3, gridspec_kw={'width_ratios': [3, 2, 2]}, facecolor='snow')
    axComplete, axFit, axBars = axes[0], axes[1], axes[2]
    fig.suptitle(f'Compression modulus')

    nSamples, colorSamples = getSamplesInfos(
        3, 2, 5,
        4, 2,
        'chocolate', 'mediumvioletred', 'steelblue',
        'lightcoral', 'crimson')

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
            x=strain[:20], y=stress[:20], yErr=stressErr[:20],
            axTitle='', yLabel='', yLim=(0, 425), xLabel='Strain', xLim=(.001, 20),
            curveColor=c, markerStyle='o', markerFColor=c, markerEColor='k',
            linearFitting=True, tableData=dataFitting)
    plotBars(axBars, 'Young modulus', dataFitting, colorSamples)

    plt.subplots_adjust(
        wspace=0, hspace=0,
        top=0.940, bottom=0.08,
        left=0.065, right=0.975)
    plt.show()

    dirSave = Path(*Path(filesPath[0]).parts[:Path(filesPath[0]).parts.index('data') + 1])
    fig.savefig(f'{dirSave}' + f'\\{fileName}' + '.png', facecolor='w', dpi=600)

    print(f'\n\n· Chart saved at\n{dirSave}.')


if __name__ == '__main__':
    folderPath = "C:/Users/petrus.kirsten/PycharmProjects/RheometerPlots/data"
    # folderPath = "C:/Users/Petrus Kirsten/Documents/GitHub/RheometerPlots/data"

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

    main(filesPath, '0St_Car_andCL-CompressionToBreakage')
