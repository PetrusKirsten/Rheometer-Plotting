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
        markerColor='crimson', markerStyle='o', markerSize=8, edgeWidth=.5, edgeColor='k',
        lineColor='dimgrey', paneColor='whitesmoke', paneEdgeColor='w',
):
    def configPlot():
        # axes.set_proj_type('ortho')
        axes.view_init(elev=9, azim=-130, roll=0)
        axes.set_box_aspect([1., 1.45, 1.5])

        for axis in axes.xaxis, axes.yaxis, axes.zaxis:
            axis._axinfo['grid'].update({
                'color': gridColor,
                'linewidth': 1,
                'linestyle': '-'})
            axis.set_label_position('lower'), axis.set_ticks_position('lower')
            axis.set_pane_color(paneColor), axis.pane.set_edgecolor(paneEdgeColor), axis.pane.set_linewidth(axWidth)
        axes.zaxis.pane.set_edgecolor('k')
        axes.zaxis.set_label_position('upper'), axes.zaxis.set_ticks_position('upper')

        # x axis
        axes.set_xlabel('')
        axes.xaxis.line.set_linewidth(axWidth)
        axes.set_xlim([-0.5, 2.5])
        axes.set_xticks([0, 1, 2])
        axes.set_xticklabels(['No St', '5% St', '10% St'])
        axes.xaxis.pane.fill = False

        # y axis
        axes.set_ylabel('')
        axes.yaxis.line.set_linewidth(axWidth)
        axes.set_ylim([-1, 5])
        axes.set_yticks([0, 1, 2, 3, 4])
        axes.set_yticklabels(['No Car', '0.5% kCar', '0.5% iCar', '1% kCar', '1% iCar'])
        # ax.invert_xaxis()
        axes.yaxis.pane.fill = False

        # z axis
        axes.set_zlabel(f'Calcium ions concentration (mM)')
        axes.zaxis.line.set_linewidth(axWidth)
        axes.set_zlim([1, 45])
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
        carWeight, carType, cacl2 = formulation
        axes.plot(carType, carWeight, cacl2,
                  c=markerColor, marker=markerStyle,
                  alpha=1, ms=markerSize, mew=edgeWidth, mec=edgeColor, ls='',
                  label='' if formulation != formulas[-1] else legend,
                  zorder=2)

        axes.plot([carType, carType], [carWeight, carWeight], [0, cacl2 - 0.06],
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


# [Car type, Starch, CaCl2]

# fivePct = [
#     [3, 1, 0], [4, 1, 0], [0, 1, 7], [0, 1, 14], [0, 1, 28],
#     [3, 1, 7], [3, 1, 14], [3, 1, 28],
#     [4, 1, 7], [4, 1, 14], [4, 1, 28]]

propForms = [
    [1, 2, 0], [1, 2, 7], [1, 2, 14], [1, 2, 28],  # St + kappaCar
    [2, 2, 0], [2, 2, 7], [2, 2, 14], [2, 2, 28]]  # St + iotaCar

doneForms = [
    [0, 1, 0], [1, 1, 0], [2, 1, 0],
    [0, 2, 0], [0, 2, 7], [0, 2, 14], [0, 2, 28],
    [3, 2, 0], [3, 2, 7], [3, 2, 14], [3, 2, 28],
    [4, 2, 0], [4, 2, 7], [4, 2, 14], [4, 2, 28],
    [3, 0, 0], [3, 0, 7], [3, 0, 14], [3, 0, 21], [3, 0, 28], [3, 0, 42],
    [4, 0, 0], [4, 0, 7], [4, 0, 14], [4, 0, 21], [4, 0, 28], [4, 0, 42],]

plt.style.use('seaborn-v0_8-ticks')
fig = plt.figure(figsize=(11, 9), constrained_layout=False)
fig.suptitle('')

ax = fig.add_subplot(111, projection='3d')

# plotFormulas3D(
#     ax, screen, legend='Car crosslinking evaluation',
#     markerColor='lightskyblue')

# plotFormulas3D(
#     ax, fivePct, legend='Suspended',
#     edgeColor='tomato', markerStyle='x',
#     markerSize=9, edgeWidth='1.3')
# plotFormulas3D(
#     ax, propForms, legend='Proposed',
#     markerColor='mediumslateblue')
plotFormulas3D(
    ax, doneForms, legend='Done',
    markerColor='coral')


plt.subplots_adjust(
    hspace=0, wspace=0.120,
    top=1, bottom=0,
    left=0, right=1)

filename = 'Hydrogels-Formulations-Plot_ScientificReport'
plt.show()
fig.savefig(f'{filename}.png', dpi=600, bbox_inches='tight')
