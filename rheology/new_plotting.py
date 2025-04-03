import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.gridspec import GridSpec


class OORecovery:
    qntsMisc = ['SegIndex', 't in s', 'h in mm', 'T in °C', 't_seg in s']
    qntsOFS = ['f in Hz', "G' in Pa", 'G" in Pa',]  # '|G*| in Pa', 'tan(δ) in -', '|η*| in mPas']
    qntsFlow = ['t in s', 't_seg in s', 'τ in Pa', 'ɣ̇ in 1/s', 'η in mPas']

    def __init__(self, label, paths):
        self.name = label
        self.paths = paths
        self.dataRaw = self._get_data()
        self.dataMean = self.dataRaw.groupby('f in Hz')[self.qntsOFS[1:]].agg(['mean', 'std']).reset_index()
        self.nReplicates = self.dataRaw['Sample'].nunique()

    def _get_data(self):
        dfs = []

        for i, file in enumerate(self.paths, start=1):
            df = pd.read_excel(file)[self.qntsOFS].dropna().reset_index(drop=True)

            breakIndex = len(df) // 2
            df_1 = df.iloc[:breakIndex].reset_index(drop=True)
            df_2 = df.iloc[breakIndex:].reset_index(drop=True)

            df_2.rename(columns={col: f"{col} | broken" for col in df_2.columns if col != 'f in Hz'}, inplace=True)

            df = df_1.merge(df_2, on='f in Hz', how='inner')
            df['Sample'] = str(i)

            dfs.append(df)

        return pd.concat(dfs, ignore_index=True)

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

def plot_errorbar(dfs):

    def configPlot(ax, yLabel, axisColor='#303030'):

        ax.spines[['top', 'bottom', 'left', 'right']].set_linewidth(0.75)
        ax.spines[['top', 'bottom', 'left', 'right']].set_color(axisColor)
        ax.tick_params(axis='both', which='both', colors=axisColor)

        ax.grid(True, which='major', axis='y', linestyle='-', linewidth=.75, color='lightgray', alpha=0.5)
        ax.grid(True, which='minor', axis='y', linestyle='-', linewidth=.5, color='lightgray', alpha=0.5)

        ax.set_xscale('log')
        # ax.set_xlim(xLim)

        ax.set_ylabel(f'{yLabel}', color=axisColor)
        ax.set_yscale('log')
        # ax.set_ylim(yLim)

    colors = plt.cm.viridis(np.linspace(0, 1, len(dfs)))

    dpi = 300
    heigth, width = 1920*2 / dpi, 1080*2 / dpi
    fig = plt.figure(figsize=(heigth, width), facecolor='snow')

    fig.canvas.manager.set_window_title('Elastic modulus')

    gs = GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])
    axPreTop, axPostTop = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    axPreBottom, axPostBottom = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])

    configPlot(axPreTop, 'Elastic modulus')
    configPlot(axPreBottom, 'Viscous modulus')

    for i, df in enumerate(dfs):
        x, y, yerr = df['f in Hz'], df[("G' in Pa", 'mean')], df[("G' in Pa", 'std')]

        axPreTop.errorbar(x, y, yerr=yerr, fmt='o', capsize=3, label=f'St CL {i * 7}', color=colors[i])
        axPreBottom.errorbar(x, y, yerr=yerr, fmt='o', capsize=3, label=f'St CL {i * 7}', color=colors[i])

    plt.xscale('log')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel("G' (Pa)")
    plt.title("G' vs Frequency with Error Bars")
    plt.legend()
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)

    plt.subplots_adjust(
        wspace=0.015, hspace=0.060,
        top=0.970, bottom=0.070,
        left=0.060, right=0.985)

    fig.savefig(
        f'Elastic and viscous moduli' + '.png',
        facecolor='w', dpi=300)
    # print(f'\n\n· Elastic and viscous moduli chart saved at:\n{dirSave}.')

    plt.show()


def main():
    fonts('C:/Users/petrus.kirsten/AppData/Local/Microsoft/Windows/Fonts/')

    plt.style.use('seaborn-v0_8-ticks')

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
        OORecovery('St CL 0', filePath[:2]).dataMean,
        OORecovery('St CL 7', filePath[2:4]).dataMean,
        OORecovery('St CL 14', filePath[4:8]).dataMean,
        OORecovery('St CL 21', filePath[8:]).dataMean
    ]

    plot_errorbar(stCL)


if __name__ == '__main__':
    main()
