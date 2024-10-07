import numpy as np
import pandas as pd
from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle


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


def constantMean(values, tolerance=100):
    """
    Encontra a região de valores praticamente constantes e calcula sua média.

    Parâmetros:
    - values (array-like): Array de valores numéricos a serem analisados.
    - tolerance (float): Tolerância para definir regiões constantes. Quanto menor, mais rigoroso.

    Retorno:
    - mean (float): Média da região constante encontrada.
    - start (int): Índice inicial da região constante.
    - end (int): Índice final da região constante.
    """
    # Calcular as diferenças entre valores consecutivos
    diffs = np.abs(np.diff(values))

    # Identificar regiões onde a diferença está abaixo do valor de tolerância
    constant_regions = diffs < tolerance

    # Encontrar os índices onde a condição é satisfeita
    indexStart, indexEnd = None, None
    max_length, current_length = 0, 0
    current_start = 0

    for i, is_constant in enumerate(constant_regions):
        if is_constant:
            if current_length == 0:
                current_start = i
            current_length += 1
        else:
            if current_length > max_length:
                max_length = current_length
                indexStart = current_start
                indexEnd = i
            current_length = 0

    # Checar se a última sequência é a maior constante
    if current_length > max_length:
        indexStart = current_start
        indexEnd = len(values) - 1

    # Se nenhuma região constante foi encontrada
    if indexStart is None or indexEnd is None:
        return None, None, None

    # Calcular a média da região constante encontrada
    mean = np.mean(values[indexStart:indexEnd + 1])

    return mean, indexStart, indexEnd


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
        freq, gPrime, gDouble = freq[nT:nT + nG], gPrime[nG:], gDouble[nG:]

        return freq, gPrime, gDouble, constantMean(gPrime), constantMean(gDouble)

    elif not recovery:
        freq, gPrime, gDouble = freq[nT:nT + nG], gPrime[:nG], gDouble[:nG]
        return freq, gPrime, gDouble, constantMean(gPrime), constantMean(gDouble)


def dataFlow(dataframe):
    time = columnRead(dataframe, 't in s')
    shearRate, shearStress = columnRead(dataframe, "GP in 1/s"), columnRead(dataframe, "Tau in Pa")
    segInd1 = dataframe.index[dataframe['Seg'] == '3-1'].to_list()[0]
    segInd2 = dataframe.index[dataframe['Seg'] == '5-1'].to_list()[0]
    timeSeg, shearRate, shearStress = time[segInd1: segInd2], shearRate[segInd1: segInd2], shearStress[segInd1: segInd2]
    timeSeg = timeSeg - timeSeg[0]

    return timeSeg, shearRate, shearStress


def plotModuli(
        ax, x, gP, gD, gPmean, gDmean, markerSize, axTitle, curveColor,
        textConfig, textLabel, textCoord, rectConfig,
        yModLimits=(1, 50000),
        showRest=False
):
    """Plots the Oscillatory Frequency Sweep Assay."""
    ax.set_title(axTitle, size=10, color='crimson')
    ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(0.75)

    ax.set_xlabel('Frequency (Hz)')
    # ax.set_xticks([])
    ax.set_xlim([x[0], x[-1] + 20])
    ax.set_xscale('log')

    ax.set_ylabel('Storage and loss moduli (Pa)')
    ax.set_ylim(yModLimits)
    ax.set_yscale('log')

    ax.plot(x[1:], gP[1:], lw=1, alpha=0.85,
            c=curveColor, marker='v', markersize=markerSize / 7, mew=0.5, mec='k',
            label=f"G' ~ {gPmean[0]:.0f} Pa", zorder=2)
    ax.plot(x[1:], gD[1:], lw=1, alpha=0.8,
            c=curveColor, marker='^', markersize=markerSize / 7, mew=0.75, mec=curveColor, mfc='w',
            label=f'G" ~ {gDmean[0]:.0f} Pa', zorder=2)

    if showRest:
        ax.text(textCoord[0], textCoord[1], s=textLabel, **textConfig)
        rect = Rectangle(*rectConfig, linewidth=1, edgecolor='k', facecolor='whitesmoke', alpha=0.75, zorder=1)
        ax.add_patch(rect)

    # ax.grid(ls='--', color='mediumaquamarine', alpha=0.5, zorder=0)
    legendLabel(ax)


def plotFlow(
        ax, x, y,
        axTitle, textSize, curveColor,
        cteShear,  # in seg
        yFlowLim=(0, 600)

):
    """Plots the Flow Shearing Assay."""
    indexCteShear = np.where(x >= cteShear)[0][0]  # find the index where time >= final cte shear
    timeCteShear, timeDecShear = np.split(x, [indexCteShear + 1])

    ax.set_title(axTitle, size=10, color='crimson')
    ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(0.75)

    ax.set_xlabel('Time')
    ax.set_xticks([cteShear])
    ax.set_xticklabels([f'{cteShear} s'])
    ax.set_xlim([x[0] - 10, x[-1] + 10])

    ax.set_ylabel('Shear stress (Pa)')
    # ax.set_yticks([])
    ax.set_ylim(yFlowLim)
    # ax.set_yscale('log')

    # Plot cte strain rate rectangle/line
    textCoord, textLabel = ((np.median(timeCteShear), 5),
                            'Constant strain rate\n$\dot{γ}=300 \,\,\, s^{-1}$')
    rectConfig = [(x[0] - cteShear, 0), 2*cteShear, yFlowLim[-1] + 50]
    ax.text(textCoord[0], textCoord[1], textLabel,
            horizontalalignment='center', verticalalignment='bottom',
            color='k', size=textSize)
    rect = Rectangle(*rectConfig, linewidth=1, linestyle='--', edgecolor='grey', facecolor='w',
                     alpha=0.7, zorder=1)
    ax.add_patch(rect)

    # Plot decreasing strain rate text
    textCoord, textLabel = ((np.median(timeDecShear), 5),
                            'Step decrease of strain rate\n$300< \dot{γ} < 0.1 \,\,\, s^{-1}$')
    ax.text(textCoord[0], textCoord[1], textLabel,
            horizontalalignment='center', verticalalignment='bottom',
            color='k', size=textSize)

    # Plot data
    ax.scatter(x, y, alpha=0.6,
               lw=.5, c=curveColor, edgecolor='k', s=30, marker='o',
               zorder=3)
    # legendLabel(ax)


def legendLabel(ax):
    """Applies consistent styling to legends in plots."""
    legend = ax.legend(frameon=True, framealpha=0.9, fancybox=False, scatterpoints=3)
    legend.get_frame().set_facecolor('w')
    legend.get_frame().set_edgecolor('whitesmoke')


# Main Plotting Configuration
def mainPlot(dataPath, color):
    df = pd.read_excel(dataPath)
    sampleName = Path(filePath).stem
    dirSave = f'{Path(filePath).parent}' + f'{Path(filePath).stem}' + '.png'

    fonts('C:/Users/petrus.kirsten/AppData/Local/Microsoft/Windows/Fonts/')
    markerSize = 50
    plt.style.use('seaborn-v0_8-ticks')
    fig, axes = plt.subplots(figsize=(18, 6), facecolor='w', ncols=3)
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
    timeAx1, storage, loss, storageMean, lossMean = dataFreqSweep(df)
    plotModuli(
        axes[0], timeAx1, storage, loss, storageMean, lossMean, markerSize,
        axTitle='Oscillatory frequency sweep assay.', curveColor=color,
        textConfig=text_properties, textLabel=text_label, textCoord=text_coord, rectConfig=rect_properties,
        showRest=False)

    # Plot 2: Flow Shearing Assay
    timeAx2, rate, stress = dataFlow(df)
    plotFlow(
        axes[1], timeAx2, stress,
        axTitle='Shear flow assay.', textSize=9.2, curveColor=color,
        cteShear=300)

    # Plot 3: Oscillatory Frequency Sweep Assay Again
    timeAx3, storage, loss, storageMean, lossMean = dataFreqSweep(df, True)
    plotModuli(
        axes[2], timeAx3, storage, loss, storageMean, lossMean, markerSize,
        'Oscillatory frequency sweep assay again.', curveColor=color,
        textConfig=text_properties, textLabel=text_label, textCoord=text_coord, rectConfig=rect_properties,
        showRest=False)

    plt.subplots_adjust(wspace=0.175, top=0.890, bottom=0.14, left=0.05, right=0.95)
    fig.savefig(dirSave, facecolor='w', dpi=600)
    plt.show()


# Run the main plotting function

filePath = ('C:/Users/petrus.kirsten/PycharmProjects/RheometerPlots/data/031024/10pct_0WSt_iCar/10pct_0WSt_iCar'
            '-RecoveryAndFlow_4.xlsx')  # -> CEBB PC
# filePath = ('C:/Users/petrus.kirsten/PycharmProjects/RheometerPlots/data/031024/10pct_0WSt/10pct_0WSt'
#             '-RecoveryAndFlow_2.xlsx')  # -> CEBB PC

# filePath = ('C:/Users/Petrus Kirsten/Documents/GitHub/RheometerPlots/data/031024/10pct_0WSt/10pct_0WSt'
#             '-RecoveryAndFlow_2.xlsx')  # personal PC

mainPlot(dataPath=filePath, color='orange')
