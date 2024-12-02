import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import MultipleLocator


def fonts(folder_path, small=10, medium=12):  # To config different fonts but it isn't working with these
    font_path = folder_path + 'HelveticaNeueThin.otf'
    helvetica_thin = FontProperties(fname=font_path)

    font_path = folder_path + 'HelveticaNeueLight.otf'
    helvetica_light = FontProperties(fname=font_path)

    font_path = folder_path + 'HelveticaNeueMedium.otf'
    helvetica_medium = FontProperties(fname=font_path)

    font_path = folder_path + 'HelveticaNeueBold.otf'
    helvetica_bold = FontProperties(fname=font_path)

    plt.rc('font', size=small)  # controls default text sizes
    plt.rc('axes', titlesize=medium)  # fontsize of the axes title
    plt.rc('axes', labelsize=medium)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=medium)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=medium)  # fontsize of the tick labels
    plt.rc('legend', fontsize=small)  # legend fontsize
    plt.rc('figure', titlesize=medium)  # fontsize of the figure title


def plotFormulas3D(
        axes, formulas, legend, axWidth=0.75, gridColor='w',
        markerColor='crimson', markerStyle='o', lineColor='dimgrey',
        paneColor='whitesmoke', edgeColor='w',
):
    def configPlot():
        axes.set_proj_type('ortho')
        axes.view_init(elev=10, azim=-18, roll=0)

        for axis in axes.xaxis, axes.yaxis, axes.zaxis:
            axis._axinfo['grid'].update({
                'color': gridColor,
                'linewidth': 1,
                'linestyle': '-'})
            axis.set_label_position('lower'), axis.set_ticks_position('lower')
            axis.set_pane_color(paneColor), axis.pane.set_edgecolor(edgeColor), axis.pane.set_linewidth(axWidth)
        axes.zaxis.pane.set_edgecolor('k')

        # x axis
        axes.set_xlabel('')
        axes.xaxis.line.set_linewidth(axWidth)
        axes.set_xlim([-0.5, 2.5])
        axes.set_xticks([0, 1, 2])
        axes.set_xticklabels(['No Car', '1% kCar', '1% iCar'])
        ax.invert_xaxis()
        axes.xaxis.pane.fill = False

        # y axis
        axes.set_ylabel('')
        axes.yaxis.line.set_linewidth(axWidth)
        axes.set_ylim([-.5, 2.5])
        axes.set_yticks([0, 1, 2])
        axes.set_yticklabels(['No St', '5% St', '10% St'])
        axes.yaxis.pane.fill = False

        # z axis
        axes.set_zlabel(f'Calcium ions concentration (mM)')
        axes.zaxis.line.set_linewidth(axWidth)
        axes.set_zlim([0, 45])
        axes.set_zticks([0, 7, 14, 21, 28, 42])
        axes.set_zticklabels(['0', '7', '14', '21', '28', '42'])

    def legendLabel():
        label = ax.legend(
            fancybox=False, frameon=False,
            framealpha=0.9, fontsize=11, ncols=1,
            loc='upper left',
            draggable=True)
        label.get_frame().set_facecolor('w')
        label.get_frame().set_edgecolor('whitesmoke')

    configPlot()
    for formulation in formulas:
        car, starch, cacl2 = formulation
        axes.plot(car, starch, cacl2,
                  c=markerColor, marker=markerStyle,
                  alpha=1, ms=8, mew=.5, mec='k', ls='',
                  label='' if formulation != formulas[-1] else legend,
                  zorder=2)

        axes.plot([car, car], [starch, starch], [0, cacl2 - 0.06],
                  color=lineColor, linestyle='-', alpha=1, lw=.6,
                  zorder=1)
    legendLabel()


def plot2D(
        axes, formulas, title,
        markerColor='crimson'
):
    def configPlot():
        axTitle = ''
        xLabel, yLabel = '', ''

        axes.set_title(axTitle, size=9, color='crimson')
        axes.spines[['top', 'bottom', 'left', 'right']].set_linewidth(0.75)

        axes.set_xlabel(f'{xLabel}')
        axes.set_xscale('linear')
        axes.set_xlim([])

        axes.set_ylabel(f'{yLabel}')
        axes.set_yscale('linear')
        axes.set_ylim([])
        axes.yaxis.set_major_locator(MultipleLocator(10))
        axes.yaxis.set_minor_locator(MultipleLocator(2.5))

    configPlot()

    # ax.errorbar(
    #     x, y, yerr=0,
    #     color=curveColor, alpha=.65,
    #     fmt='D' if '21' in sampleName else 'o',
    #     markersize=7 if '21' in sampleName else 6.5,
    #     mfc=curveColor, mec='#383838', mew=.75,
    #     linestyle='',
    #     label=f'{sampleName}', zorder=3)


# [Car, St, CaCl2]
finalForms = [
    [0, 2, 0], [0, 2, 7], [0, 2, 14], [0, 2, 28],  # St
    [1, 2, 0], [1, 2, 7], [1, 2, 14], [1, 2, 28],  # St + kappaCar
    [2, 2, 0], [2, 2, 7], [2, 2, 14], [2, 2, 28]]  # St + iotaCar

fivePct = [
    [0, 1, 7], [0, 1, 14], [0, 1, 28],  # St
    [1, 1, 7], [1, 1, 14], [1, 1, 28],  # St + kappaCar
    [2, 1, 7], [2, 1, 14], [2, 1, 28]]  # St + iotaCar

noCaFivePct = [
    [0, 1, 0],  # St
    [1, 1, 0],  # St + kappaCar
    [2, 1, 0]]  # St + iotaCar

cars = [
    [1, 0, 0], [1, 0, 7], [1, 0, 14], [1, 0, 21], [1, 0, 28], [1, 0, 42],  # St + kappaCar
    [2, 0, 0], [2, 0, 7], [2, 0, 14], [2, 0, 21], [2, 0, 28], [2, 0, 42]]  # St + iotaCar

plt.style.use('seaborn-v0_8-ticks')
fig = plt.figure(figsize=(8.8, 9), constrained_layout=False)
fig.suptitle('')

ax = fig.add_subplot(111, projection='3d')
fig.tight_layout()


plotFormulas3D(
    ax, cars, legend='Car crosslinking evaluation',
    markerColor='lightskyblue')

plotFormulas3D(
    ax, noCaFivePct, legend='Suspended formulations',
    markerColor='lightcoral')
plotFormulas3D(
    ax, fivePct, legend='',
    markerColor='lightcoral', markerStyle='x')

plotFormulas3D(
    ax, finalForms, legend='Chosen formulations',
    markerColor='lightgreen')

filename = 'Hydrogels-Formulations-Plot-NEW'
plt.show()
fig.savefig(f'{filename}.png', dpi=600, bbox_inches='tight')
