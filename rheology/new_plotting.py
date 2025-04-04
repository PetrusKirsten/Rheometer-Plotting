import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.gridspec import GridSpec
from scipy.optimize import curve_fit


def fonts(folder_path, s=12, m=14):

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


def powerLaw(x, k, n):
    return k * (x ** n)


class OORecovery:
    qntsMisc = ['SegIndex', 't in s', 'h in mm', 'T in °C', 't_seg in s']
    qntsOFS = ['f in Hz', "G' in Pa", 'G" in Pa',]  # '|G*| in Pa', 'tan(δ) in -', '|η*| in mPas']
    qntsOFSbroken = qntsOFS[1:] + [param + ' | broken' for param in qntsOFS[1:]]
    qntsFlow = ['t in s', 't_seg in s', 'τ in Pa', 'ɣ̇ in 1/s', 'η in mPas']

    def __init__(self, paths, label, color):
        self.paths, self.label, self.color = paths, label, color

        self.dataRaw = self._get_data()
        self.nReplicates = self.dataRaw['Sample'].nunique()

        self.dataMean = self.dataRaw.groupby('f in Hz')[self.qntsOFSbroken].agg(['mean', 'std']).reset_index()

    def _get_data(self):
        dfs = []

        for i, file in enumerate(self.paths, start=1):

            df = pd.read_excel(file)[self.qntsOFS].dropna().reset_index(drop=True)

            breakIndex = len(df) // 2
            df_1, df_2 = df.iloc[:breakIndex].reset_index(drop=True), df.iloc[breakIndex:].reset_index(drop=True)
            df_2.rename(columns={col: f"{col} | broken" for col in df_2.columns if col != 'f in Hz'}, inplace=True)

            df = df_1.merge(df_2, on='f in Hz', how='inner')
            df['Sample'] = str(i)

            dfs.append(df)

        return pd.concat(dfs, ignore_index=True)

    def _powerLaw(self, modulus):
        """
        :type modulus: str
        :rtype: object
        """

        x, y = self.dataMean['f in Hz'][:16], None

        if modulus == "G'":
            y = self.dataMean["G' in Pa", 'mean'][:16]

        elif modulus == 'G"':
            y = self.dataMean['G" in Pa', 'mean'][:16]

        self.fitting = curve_fit(powerLaw, x, y)

        return self.fitting[0], self.fitting[1], np.sqrt(np.diag(self.fitting[1]))

def plotOFS(samples):

    def configFigure(width=4000, heigth=4000, dpi=300):

        def configAxes(ax, xLabel, yLabel, xLim=(0.08, 140), yLim=(1, 10**4), axisColor='#303030'):

            ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(.75)
            ax.spines[['top', 'bottom', 'left', 'right']].set_color(axisColor)
            ax.tick_params(axis='both', which='both', colors=axisColor)

            ax.grid(True, which='major', axis='y', linestyle='-', linewidth=.75, color='lightgray', alpha=.5)
            ax.grid(True, which='minor', axis='y', linestyle='-', linewidth=.5, color='lightgray', alpha=.5)

            ax.set_xlabel(f'{xLabel}', color=axisColor), ax.set_ylabel(f'{yLabel}', color=axisColor)
            ax.set_xscale('log'), ax.set_yscale('log')
            ax.set_xlim(xLim), ax.set_ylim(yLim)

            if xLabel == '':
                ax.set_xticklabels([])
                ax.tick_params(axis='x', which='both', direction='in')

            if yLabel == '':
                ax.set_yticklabels([])
                ax.tick_params(axis='y', which='both', direction='in')

        fonts('C:/Users/petrus.kirsten/AppData/Local/Microsoft/Windows/Fonts/')
        plt.style.use('seaborn-v0_8-ticks')

        figure = plt.figure(figsize=(width / dpi, heigth / dpi), facecolor='snow')
        figure.canvas.manager.set_window_title(windowTitle)

        gs = GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])
        axes = [figure.add_subplot(gs[line, column]) for line in range(2) for column in range(2)]

        configAxes(axes[0], '', 'Elastic modulus')
        configAxes(axes[1], '', '')
        configAxes(axes[2], 'Frequency (Hz)', 'Viscous modulus')
        configAxes(axes[3], 'Frequency (Hz)', '')

        return figure, axes, dpi

    def addLegend(ax):

        legend = ax.legend(
            loc='lower right',
            fancybox=False,
            frameon=True,
            framealpha=0.9,
            fontsize=11,
            markerscale=1.3)

        legend.get_frame().set_facecolor('w')
        legend.get_frame().set_edgecolor('whitesmoke')

    def plotData(ax, label, color, x, y, yerr):

        lowError, highError = np.abs(yerr / y) <= 1.5, np.abs(yerr / y) > 1.5
        transp = .5

        # Plot error lines
        # markers with relative error <= 50%
        ax.errorbar(
            x[lowError], y[lowError], yerr[lowError],
            label='', fmt='none', markersize=0,
            color=color, mfc=color,
            alpha=transp,
            mec='w', mew=.75,
            capsize=2.5, capthick=1, linestyle='', lw=1,
            zorder=2)
        # markers with relative error > 50%
        ax.errorbar(
            x[highError], y[highError], yerr[highError],
            label='', fmt='none', markersize=0,
            color=color, mfc=color,
            alpha=transp / 10,
            mec='w', mew=.75,
            capsize=2.5, capthick=1, linestyle='', lw=1,
            zorder=2)

        # Plot markers
        # markers with relative error <= 50%
        ax.errorbar(
            x[lowError], y[lowError], 0,
            fmt='o', markersize=6,
            label=label,
            color=color, mfc=color,
            alpha=transp,
            mec=color, mew=.75,
            capsize=0, capthick=0, linestyle='', lw=0,
            zorder=3)
        # markers with relative error > 50%
        ax.errorbar(
            x[highError], y[highError], 0,
            fmt='o', markersize=6,
            label='',
            color=color, mfc=color,
            alpha=transp / 10,
            mec=color, mew=.75,
            capsize=0, capthick=0, linestyle='', lw=0,
            zorder=3)

    windowTitle = 'Viscoelastic recovery'
    fig, [axTL, axTR, axBL, axBR], resolution = configFigure()

    for sample in samples:
        df, l, c = sample.dataMean, sample.label, sample.color

        plotData(
            axTL, l, c,
            df['f in Hz'], df[("G' in Pa", 'mean')], df[("G' in Pa", 'std')])
        plotData(
            axTR, l, c,
            df['f in Hz'], df[("G' in Pa | broken", 'mean')], df[("G' in Pa | broken", 'std')])

        plotData(
            axBL, l, c,
            df['f in Hz'], df[('G" in Pa', 'mean')], df[('G" in Pa', 'std')])
        plotData(
            axBR, l, c,
            df['f in Hz'], df[('G" in Pa | broken', 'mean')], df[('G" in Pa | broken', 'std')])

    addLegend(axTL)

    plt.subplots_adjust(
        wspace=0.010, hspace=0.015,
        top=0.970, bottom=0.070,
        left=0.060, right=0.985)

    fig.savefig(
        f'{windowTitle} ' + f'{l.split('CL')[0].strip()}' + ' CL' + '.svg',
        facecolor='w', dpi=resolution)

def main():

    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)

    folderPath = "D:/Documents/GitHub/Rheometer-Plotting/data/by sample"

    filePath = [
        # St
        folderPath + "/10St/10_0WSt-viscRec_1.xlsx",
        folderPath + "/10St/10_0WSt-viscRec_2.xlsx",

        # St/CL_7
        folderPath + "/10St_CL_7/10_0St_CL-recovery-1.xlsx",
        # folderPath + "/10St_CL_7/10_0St_CL-recovery-2.xlsx",
        folderPath + "/10St_CL_7/10_0St_CL-recovery-3.xlsx",

        # St/CL_14
        folderPath + "/10St_CL_14/St_CL_14-viscoelasticRecovery-1.xlsx",
        folderPath + "/10St_CL_14/St_CL_14-viscoelasticRecovery-2.xlsx",
        folderPath + "/10St_CL_14/St_CL_14-viscoelasticRecovery-3.xlsx",
        folderPath + "/10St_CL_14/St_CL_14-viscoelasticRecovery-4.xlsx",

        # St/CL_28
        folderPath + "/10St_CL_28/St_CL_28-viscoelasticRecovery-1.xlsx",
        folderPath + "/10St_CL_28/St_CL_28-viscoelasticRecovery-2.xlsx",
        folderPath + "/10St_CL_28/St_CL_28-viscoelasticRecovery-3.xlsx",
        folderPath + "/10St_CL_28/St_CL_28-viscoelasticRecovery-4.xlsx",
    ]

    stCL = [
        OORecovery(filePath[:2], 'St CL 0', '#E1C96B'),
        OORecovery(filePath[2:4], 'St CL 7', '#FFE138'),
        OORecovery(filePath[4:8], 'St CL 14', '#F1A836'),
        OORecovery(filePath[8:], 'St CL 21', '#E36E34'),
    ]

    plotOFS(stCL)

    plt.show()


if __name__ == '__main__':
    main()
