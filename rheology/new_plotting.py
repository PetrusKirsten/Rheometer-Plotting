import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

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


if __name__ == '__main__':

    def plot_errorbar(dfs):
        colors = plt.cm.viridis(np.linspace(0, 1, len(dfs)))
        plt.figure(figsize=(10, 6))

        for i, df in enumerate(dfs):
            x, y, yerr = df['f in Hz'], df[("G' in Pa", 'mean')], df[("G' in Pa", 'std')]

            plt.errorbar(x, y, yerr=yerr, fmt='o', capsize=3, label=f'St CL {i * 7}', color=colors[i])

        plt.xscale('log')
        plt.xlabel('Frequency (Hz)')
        plt.ylabel("G' (Pa)")
        plt.title("G' vs Frequency with Error Bars")
        plt.legend()
        plt.grid(True, which="both", linestyle="--", linewidth=0.5)

        plt.show()

    pd.set_option('display.max_rows', None), pd.set_option('display.max_columns', None), pd.set_option('display.width', None)

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
