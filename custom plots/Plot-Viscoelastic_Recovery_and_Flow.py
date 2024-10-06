import numpy as np
import pandas as pd
from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle


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
    plt.rc('axes', titlesize=m)  # fontsize of the axes title
    plt.rc('axes', labelsize=m)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=m)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=m)  # fontsize of the tick labels
    plt.rc('legend', fontsize=s)  # legend fontsize
    plt.rc('figure', titlesize=m)  # fontsize of the figure title


def columnRead(dataframe, key):
    """
    :param dataframe: the Pandas dataframe to read
    :param key: the column label/title
    :return: the values from column in a numpy array
    """
    return dataframe[key].to_numpy()


def dataFreqSweep(dataframe, recovery=False):
    freq = columnRead(dataframe, 'f in Hz')
    gPrime, gDouble = columnRead(dataframe, "G' in Pa"), columnRead(dataframe, "G'' in Pa")
    gPrime, gDouble = gPrime[gPrime > 0], gDouble[gDouble > 0]

    nG = gPrime.shape[0] // 2
    nT = dataframe.index[dataframe['Seg'] == '2-1'].to_list()[0]
    if recovery:
        nT = dataframe.index[dataframe['Seg'] == '5-1'].to_list()[0]
        return freq[nT:nT + nG], gPrime[nG:], gDouble[nG:]
    elif not recovery:
        return freq[nT:nT + nG], gPrime[:nG], gDouble[:nG]


def dataFlow(dataframe):
    time = columnRead(dataframe, 't in s')
    shearRate, shearStress = columnRead(dataframe, "GP in 1/s"), columnRead(dataframe, "Tau in Pa")
    segInd1 = dataframe.index[dataframe['Seg'] == '3-1'].to_list()[0]
    segInd2 = dataframe.index[dataframe['Seg'] == '5-1'].to_list()[0]
    timeSeg, shearRate, shearStress = time[segInd1: segInd2], shearRate[segInd1: segInd2], shearStress[segInd1: segInd2]
    timeSeg = timeSeg - timeSeg[0]

    return timeSeg, shearRate, shearStress


def plotFreqSweep(
        ax, x, gP, gD, markerSize,
        axTitle, textConfig, textLabel, textCoord, rectConfig,
        yFlowLimits, tickRight=False, showRest=False
):
    """Plots the Oscillatory Frequency Sweep Assay."""
    ax2 = None
    yFreqLimits = [1, 50000]
    ax.set_title(axTitle, size=10, color='crimson')
    ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(0.75)

    ax.set_xlabel('Frequency (Hz)')
    # ax.set_xticks([])
    ax.set_xlim([x[0], x[-1] + 20])
    ax.set_xscale('log')

    if tickRight:
        ax.set_ylabel('Shear stress (Pa)')
        ax.set_ylim(yFlowLimits)
        # ax.set_yscale('log')

        ax2 = ax.twinx()
        ax2.spines[['top', 'bottom', 'left', 'right']].set_linewidth(0.75)
        ax2.set_ylabel('Storage and loss moduli (Pa)')
        ax2.set_ylim(yFreqLimits)
        ax2.set_yscale('log')
        ax2.plot(x[1:], gP[1:], lw=1, alpha=0.85,
                 c='dodgerblue', marker='v', markersize=markerSize / 7, mew=0.5, mec='k',
                 label="G'", zorder=2)
        ax2.plot(x[1:], gD[1:], lw=1, alpha=0.8,
                 c='lightskyblue', marker='^', markersize=markerSize / 7, mew=0.5, mec='k',
                 label='G"', zorder=2)
        legendLabel(ax2)

    else:
        ax.set_ylabel('Storage and loss moduli (Pa)')
        ax.set_ylim(yFreqLimits)
        ax.set_yscale('log')
        ax.plot(x[1:], gP[1:], lw=1, alpha=0.85,
                c='dodgerblue', marker='v', markersize=markerSize / 7, mew=0.5, mec='k',
                label="G'", zorder=2)
        ax.plot(x[1:], gD[1:], lw=1, alpha=0.8,
                c='lightskyblue', marker='^', markersize=markerSize / 7, mew=0.5, mec='k',
                label='G"', zorder=2)

    if showRest:
        ax.text(textCoord[0], textCoord[1], s=textLabel, **textConfig)
        rect = Rectangle(*rectConfig, linewidth=1, edgecolor='k', facecolor='whitesmoke', alpha=0.75, zorder=1)
        ax.add_patch(rect)

    # ax.grid(ls='--', color='mediumaquamarine', alpha=0.5, zorder=0)
    legendLabel(ax)


def plotFlow(
        ax, x, y,
        axTitle, textSize
):
    """Plots the Flow Shearing Assay."""
    cteShear = 300  # in seg
    indexCteShear = np.where(x == 300)[0][0]  # find the index where time = final cte shear
    timeCteShear, timeDecShear = np.split(x, [indexCteShear + 1])

    ax.set_title(axTitle, size=10, color='crimson')
    ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(0.75)

    ax.set_xlabel('Time')
    ax.set_xticks([cteShear])
    ax.set_xticklabels([f'{cteShear} s'])
    ax.set_xlim([x[0] - 10, x[-1] + 10])

    ax.set_yticks([])
    yFlowLim = [round((y.min() - 10) / 10) * 10, round((y.max() + 10) / 10) * 10]
    yFlowLim = [1, 490]
    ax.set_ylim(yFlowLim)
    # ax.set_yscale('log')

    # Plot cte strain rate rectangle/line
    textCoord = (np.median(timeCteShear), 5)
    textLabel = 'Constant strain rate\n$\dot{γ}=300 \,\,\, s^{-1}$'
    rectConfig = [(x[0] - cteShear, 0), 2*cteShear, y.max() + 100]
    ax.text(textCoord[0], textCoord[1], textLabel,
            horizontalalignment='center', verticalalignment='bottom', color='k', size=textSize)
    rect = Rectangle(*rectConfig, linewidth=1, linestyle='--', edgecolor='grey', facecolor='w', alpha=0.7, zorder=1)
    ax.add_patch(rect)
    # TODO: alterar tamanho das ticklabels e colocar ylim cte em plotFlow
    # Plot decreasing strain rate text
    textCoord = (np.median(timeDecShear), 5)
    textLabel = 'Step decrease of strain rate\n$300< \dot{γ} < 0.1 \,\,\, s^{-1}$'
    ax.text(textCoord[0], textCoord[1], textLabel,
            horizontalalignment='center', verticalalignment='bottom', color='k', size=textSize)

    ax.scatter(x, y, alpha=0.6,
               lw=.5, c='hotpink', edgecolor='k', s=30, marker='o',
               zorder=3)
    # legendLabel(ax)
    return yFlowLim  # to use it in plotFreqSweep function


def legendLabel(ax):
    """Applies consistent styling to legends in plots."""
    legend = ax.legend(loc=4, frameon=True, framealpha=1, fancybox=False, scatterpoints=3)
    legend.get_frame().set_facecolor('w')
    legend.get_frame().set_edgecolor('whitesmoke')


# Main Plotting Configuration
def mainPlot(dataPath):
    df = pd.read_excel(dataPath)
    sampleName = Path(filePath).stem
    dirSave = f'{Path(filePath).parent}' + f'{Path(filePath).stem}' + '.png'

    fonts('C:/Users/petrus.kirsten/AppData/Local/Microsoft/Windows/Fonts/')
    markerSize = 50
    plt.style.use('seaborn-v0_8-ticks')
    fig, axes = plt.subplots(figsize=(16, 6), facecolor='w', ncols=3)
    fig.suptitle(f'Rheometry assay protocol to evaluate viscoelastic recovery | {sampleName}\n\n')

    # Common parameters
    text_coord = (0, 100)
    text_label = 'Rest and\nSet $T=37\,^oC$\nfor $180\,\,s$'
    text_properties = {'horizontalalignment': 'center',
                       'verticalalignment': 'top',
                       'color': 'k',
                       'size': 9.2}
    rect_properties = [(0, 0), 10, 100000]

    # Plot 1: Oscillatory Frequency Sweep Assay
    timeAx1, storage, loss = dataFreqSweep(df)
    plotFreqSweep(
        axes[0], timeAx1, storage, loss, markerSize,
        'Oscillatory frequency sweep assay.',
        text_properties, text_label, text_coord, rect_properties,
        yFlowLimits=None, showRest=False)

    # Plot 2: Flow Shearing Assay
    timeAx2, rate, stress = dataFlow(df)
    yLim = plotFlow(
        axes[1], timeAx2, stress,
        axTitle='Shear flow assay.', textSize=9.2)

    # Plot 3: Oscillatory Frequency Sweep Assay Again
    timeAx3, storage, loss = dataFreqSweep(df, True)
    plotFreqSweep(
        axes[2], timeAx3, storage, loss, markerSize,
        'Oscillatory frequency sweep assay again.',
        text_properties, text_label, text_coord, rect_properties,
        yFlowLimits=yLim, tickRight=True, showRest=False)

    plt.subplots_adjust(wspace=0.0, top=0.890, bottom=0.14, left=0.05, right=0.95)
    fig.savefig(dirSave, facecolor='w', dpi=600)
    plt.show()


# Run the main plotting function
# filePath = ('C:/Users/petrus.kirsten/PycharmProjects/RheometerPlots/data/031024/10pct_0WSt/10pct_0WSt'
#             '-RecoveryAndFlow_1.xlsx') -> CEBB PC
filePath = ('C:/Users/Petrus Kirsten/Documents/GitHub/RheometerPlots/data/031024/10pct_0WSt/10pct_0WSt'
            '-RecoveryAndFlow_2.xlsx')  # personal PC

mainPlot(dataPath=filePath)
