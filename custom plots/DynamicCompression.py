import numpy as np
import pandas as pd
from pathlib import Path

from scipy.interpolate import interp1d
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
        ax.xaxis.set_minor_locator(MultipleLocator(1))

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
            fmt=markerStyle, markersize=3.5, mfc=markerFColor, mec=markerFColor, mew=markerEWidth,
            capsize=0, lw=.5, linestyle=lineStyle,
            label=f'{sampleName}',
            zorder=4)
        legendLabel()


def main(dataPath):
    fonts('C:/Users/petrus.kirsten/AppData/Local/Microsoft/Windows/Fonts/')
    plt.style.use('seaborn-v0_8-ticks')
    fig, axForce = plt.subplots(
        figsize=(18, 8), ncols=1,
        gridspec_kw={'width_ratios': [1.]}, facecolor='snow')

    fig.suptitle(f'Oscilatory compression')
    fTitle, fLimits = f'Stress (Pa)', (0, 210)
    xTitle, xLimits = f'Time (s)', (0, 75)

    nSamples, colorSamples = getSamplesInfos(
        3, 2, 5,
        4, 4,
        'darkgray', 'crimson', 'mediumblue',
        'mediumorchid', 'rebeccapurple')

    data, labels = getSamplesData(dataPath, nSamples)

    dataList = {
        labels[0]: ([], [], []),
        labels[1]: ([], [], []),
        labels[2]: ([], [], []),
        labels[3]: ([], [], []),
        labels[4]: ([], [], [])
    }

    means_hMax = []

    for key, (x, height, f) in dataList.items():
        x.append(data[f'{key} time'])
        height.append(data[f'{key} height'])
        f.append(data[f'{key} force'])

    for key, c in zip(dataList, colorSamples):
        time, height, stress = dataList[key][0][0], dataList[key][1][0], dataList[key][2][0]
        hMax_mean, stressMean, stressErrMean = np.max(height, axis=1), np.mean(stress, axis=0), np.std(stress, axis=0)
        means_hMax.append(hMax_mean.tolist())

        plotCompression(
            ax=axForce, axisColor='#303030',
            x=np.mean(time, axis=0), y=stressMean / (35 / 1000) + 7.5, yErr=stressErrMean / (35 / 1000),
            axTitle='', yLabel=fTitle, yLim=fLimits, xLabel=xTitle, xLim=xLimits,
            curveColor=c, markerStyle='o', markerFColor=c, markerEColor='k',
            sampleName=f'{key}')
    print(means_hMax)
    plt.subplots_adjust(
        wspace=0.175,
        top=0.890, bottom=0.14,
        left=0.05, right=0.95)
    plt.show()

    fileName = '0St_Car_CL-DynamicCompression'
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
