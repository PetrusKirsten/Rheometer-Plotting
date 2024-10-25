import numpy as np
import pandas as pd
from pathlib import Path
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


def getSamplesData(dataPath, n5st=0, n10St=0, nIc=0, n5Kc=0, n10Kc=1):
    def getSegments(dataframe):
        time = dataframe['t in s'].to_numpy()
        height = dataframe['h in mm'].to_numpy()
        force = dataframe['Fn in N'].to_numpy()

        seg2, seg3, seg4 = (  # Identifying the job segments in the lists
            dataframe.index[dataframe['SegIndex'] == seg].to_list()[0] for seg in ['1|1', '11|1', '12|1'])
        segments = lambda arr: (arr[seg2:seg3])  # Slice segments
        segmentsBreakage = lambda arr: (arr[seg3:seg4])  # Slice segments

        return {
            'time': segments(time) - segments(time)[0],
            'height': (1 - segments(height) / segments(height).max())*100,
            'force': segments(force),
            'time to break': segmentsBreakage(time) - segmentsBreakage(time)[0],
            'height to break': (1 - segmentsBreakage(height) / segmentsBreakage(height).max())*100,
            'force to break': segmentsBreakage(force)}

    samples = {'10% WSt kCar': []}  # Store data for each sample type
    sample_labels = ['10% WSt kCar'] * n10Kc  # Determine sample types for each path

    for sample_type, path in zip(sample_labels, dataPath):  # Read data and categorize based on sample type
        df = pd.read_excel(path)
        segments = getSegments(df)
        samples[sample_type].append(segments)

    dict_tempSweeps = {}  # Initialize dictionaries to hold the results

    for sample_type in samples:  # Populate dictionaries with consolidated sample data
        dict_tempSweeps[f'{sample_type} time'] = [s['time'] for s in samples[sample_type]]
        dict_tempSweeps[f'{sample_type} height'] = [s['height'] for s in samples[sample_type]]
        dict_tempSweeps[f'{sample_type} force'] = [s['force'] for s in samples[sample_type]]
        dict_tempSweeps[f'{sample_type} time to break'] = [s['time to break'] for s in samples[sample_type]]
        dict_tempSweeps[f'{sample_type} height to break'] = [s['height to break'] for s in samples[sample_type]]
        dict_tempSweeps[f'{sample_type} force to break'] = [s['force to break'] for s in samples[sample_type]]

    return dict_tempSweeps


def plotCompression(nSamples, sampleName,
                    ax, x, y,
                    axTitle, yLabel, yLim, xLabel, xLim, axisColor,
                    curveColor, markerStyle, markerFColor, markerEColor, markerEWidth=0.5,
                    strain=False, lineStyle='', logScale=False):
    def legendLabel():
        """Applies consistent styling to legends in plots."""
        legend = ax.legend(fancybox=False, frameon=True, framealpha=0.9, fontsize=9)
        legend.get_frame().set_facecolor('w')
        legend.get_frame().set_edgecolor('whitesmoke')

    def configPlot():
        ax.set_title(axTitle, size=9, color='crimson')
        ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(0.75)
        ax.grid(True, which='both', axis='x', linestyle='-', linewidth=0.5, color='lightsteelblue', alpha=0.5)

        ax.set_xlabel(f'{xLabel}')
        # ax.set_xscale('log' if logScale else 'linear')
        ax.set_xlim(xLim)
        ax.xaxis.set_minor_locator(MultipleLocator(1))

        ax.set_ylabel(f'{yLabel}', color=axisColor)
        ax.set_yscale('log' if logScale else 'linear')
        ax.set_ylim(yLim)
        ax.tick_params(axis='y', colors=axisColor, which='both')

        if strain:
            ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}%"))
            ax.plot(
                x[0], y[0],
                color=curveColor, alpha=0.8, lw=1.5, linestyle=':',
                label=f'{sampleName}', zorder=3)

        else:
            ax.errorbar(
                x[0], y[0], 0,
                color=curveColor, alpha=0.35,
                fmt=markerStyle, markersize=6, mfc=markerFColor, mec=markerEColor, mew=markerEWidth,
                capsize=0, lw=1, linestyle=lineStyle,
                label=f'{sampleName}',
                zorder=4)

    # x = np.mean(x, axis=0)
    # yerr = np.std(y, axis=0)
    # y = np.mean(y, axis=0)
    configPlot()


def main(dataPath):
    fonts('C:/Users/petrus.kirsten/AppData/Local/Microsoft/Windows/Fonts/')

    fileName = '10pct_0WSt_and_Car-DynamicCompression'
    dirSave = Path(*Path(filePath[0]).parts[:Path(filePath[0]).parts.index('data') + 1])

    fig, axForce = plt.subplots(figsize=(12, 5), facecolor='w', ncols=1)
    axHeight = axForce.twinx()

    plt.style.use('seaborn-v0_8-ticks')
    fig.suptitle(f'Oscilatory compression')
    forceColor, strainColor = 'tomato', 'mediumseagreen'

    st5_nSamples, st10_nSamples, ic10_nSamples, kc10_nSamples = 0, 0, 0, 1
    data = getSamplesData(dataPath)

    x_10kc, f_10kc, s_10kc = (
        data['10% WSt kCar time'],
        data['10% WSt kCar force'],
        data['10% WSt kCar height'])

    (fTitle, fLimits,
     hTitle, hLimits,
     xTitle, xLimits) = (
        f'Force (N)', (0, 3),
        f'Strain', (9, -1),
        f'Time (s)', (0, x_10kc[0][-1]))

    plotCompression(
        nSamples=kc10_nSamples,
        ax=axForce, x=x_10kc, y=f_10kc, axisColor=forceColor,
        axTitle='', yLabel=fTitle, yLim=fLimits, xLabel=xTitle, xLim=xLimits,
        curveColor=forceColor, markerStyle='o', markerFColor=forceColor, markerEColor='k',
        sampleName=f'')
    plotCompression(
        nSamples=kc10_nSamples,
        ax=axHeight, x=x_10kc, y=s_10kc, axisColor=strainColor, strain=True,
        axTitle='', yLabel=hTitle, yLim=hLimits, xLabel=xTitle, xLim=xLimits,
        curveColor=strainColor, markerStyle='o', markerFColor=strainColor, markerEColor='k',
        sampleName=f'')

    # axForce.spines['left'].set_color(forceColor)
    # axForce.spines['right'].set_color(strainColor)

    # plt.subplots_adjust(wspace=0.175, top=0.890, bottom=0.14, left=0.05, right=0.95)
    plt.tight_layout()
    plt.show()
    fig.savefig(f'{dirSave}' + f'\\{fileName}' + '.png', facecolor='w', dpi=600)

    print(f'\n\n· Chart saved at\n{dirSave}.')


if __name__ == '__main__':
    folderPath = "C:/Users/petrus.kirsten/PycharmProjects/RheometerPlots/data"
    # folderPath = "C:/Users/Petrus Kirsten/Documents/GitHub/RheometerPlots/data"
    filePath = [
        folderPath + "/old/200924/7PSt_2_Compression.xlsx"
    ]

    main(dataPath=filePath)
