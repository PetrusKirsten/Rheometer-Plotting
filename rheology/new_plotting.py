import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

class OORecovery:
    def __init__(
            self,
            label,
            paths
    ):
        """
        :type label: str
        """

        def getData():
            dfs = []
            for i, file in enumerate(paths, start=1):
                df = pd.read_excel(file)[qntsOFS].dropna().reset_index(drop=True)

                breakageIndex = len(df) // 2
                df_1 = df.iloc[:breakageIndex].reset_index(drop=True)  # Before break
                df_2 = df.iloc[breakageIndex:].reset_index(drop=True)  # After break

                df_2 = df_2.rename(
                    columns={
                        col: col + " | broken"
                        for col in df_2.columns
                        if col != 'f in Hz'})

                df = df_1.merge(df_2, on='f in Hz')
                df['Sample'] = f'{i}'

                dfs.append(df)

            df_replicates = pd.concat(dfs, ignore_index=True)

            return df_replicates

        qntsMisc = [
            'SegIndex',
            't in s',
            'h in mm',
            'T in °C',
            't_seg in s']
        qntsOFS = [
            'f in Hz',
            "G' in Pa",
            'G" in Pa',
            '|G*| in Pa',
            'tan(δ) in -',
            '|η*| in mPas']
        qntsFlow = [
            't in s',
            't_seg in s'
            'τ in Pa',
            'ɣ̇ in 1/s',
            'η in mPas']

        self.name = label
        self.paths = paths

        self.dataRaw = getData()
        self.dataMean = self.dataRaw.groupby('f in Hz')[qntsOFS[1:]].agg(['mean', 'std']).reset_index()

        self.nReplicates = int(self.dataRaw['Sample'].max())


if __name__ == '__main__':

    pd.set_option('display.max_rows', None), pd.set_option('display.max_columns', None), pd.set_option('display.width', None)

    folderPath = "D:/Documents/GitHub/Rheometer-Plotting/data/by sample"   # New Personal
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
        OORecovery('St CL 0', filePath[:2]),
        OORecovery('St CL 7', filePath[2:4]),
        OORecovery('St CL 14', filePath[4:8]),
        OORecovery('St CL 21', filePath[8:])
    ]

    noCL = stCL[0]
    plt.errorbar(noCL.dataMean['f in Hz'], noCL.dataMean["G' in Pa"]['mean'], noCL.dataMean["G' in Pa"]['std'])
    plt.show()

    print()