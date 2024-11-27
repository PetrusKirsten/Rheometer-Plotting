import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties


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


def plotFormulas(
        axes, formulas, title, typeCar=None,
        markerColor='crimson', paneColor='whitesmoke', edgeColor='gainsboro'
):
    axes.set_proj_type('ortho')
    axes.view_init(elev=16, azim=155, roll=0)
    axes.grid(False)
    for axis in axes.xaxis, axes.yaxis, axes.zaxis:
        axis.set_label_position('lower'), axis.set_ticks_position('lower')
        axis.set_pane_color(paneColor), axis.pane.set_edgecolor(edgeColor)

    axes.set_title(title)

    # x axis
    axes.set_xlabel('')
    axes.set_xlim([-0.5, 2.5])
    axes.set_xticks([0, 1, 2])
    axes.set_xticklabels(['No Car', 'kCar', 'iCar'])
    # ax.invert_xaxis()

    # y axis
    axes.set_ylabel('Starch')
    axes.set_ylim([-0.5, 2.5])
    axes.set_yticks([0, 1, 2])
    axes.set_yticklabels(['No St', '5%', '10%'])
    axes.invert_yaxis()

    # z axis
    axes.set_zlabel(f'Calcium ions concentration (mM)')
    axes.set_zlim([-2, 45])
    axes.set_zticks([0, 7, 14, 21, 28, 42])
    axes.set_zticklabels(['0', '7', '14', '21', '28', '42'])

    for formulation in formulas:
        car, starch, cacl2 = formulation
        axes.plot(car, starch, cacl2,
                  c=markerColor, alpha=1, ms=12, marker='o', mew=0.5, mec='k',
                  label='', zorder=2)

        axes.plot([car, car], [starch, starch], [0, cacl2 - 0.06],
                  color='#383838', linestyle='-', alpha=0.8, lw=0.6,
                  zorder=1)


finalForms = [
    [0, 2, 0], [0, 2, 7], [0, 2, 14], [0, 2, 28],  # St
    [1, 2, 0], [1, 2, 7], [1, 2, 14], [1, 2, 28],  # St + kappaCar
    [2, 2, 0], [2, 2, 7], [2, 2, 14], [2, 2, 28],  # St + iotaCar
]

fivePct = [
    [0, 1, 0],  # St
    [1, 1, 0],  # St + kappaCar
    [2, 1, 0],  # St + iotaCar
]

cars = [
    [1, 0, 0], [1, 0, 7], [1, 0, 14], [1, 0, 21], [1, 0, 28], [1, 0, 42],  # St + kappaCar
    [2, 0, 0], [2, 0, 7], [2, 0, 14], [2, 0, 21], [2, 0, 28], [2, 0, 42],  # St + iotaCar
]

plt.style.use('seaborn-v0_8-ticks')
fig = plt.figure(figsize=(8, 7), constrained_layout=False)
fig.suptitle('')

ax = fig.add_subplot(111, projection='3d')
fig.tight_layout()

plotFormulas(
    ax, finalForms, '', markerColor='orange')
plotFormulas(
    ax, fivePct, '', markerColor='gainsboro')
plotFormulas(
    ax, cars, '',  markerColor='mediumspringgreen')

filename = 'Hydrogels-Formulations-Plot-NEW'
plt.show()
fig.savefig(f'{filename}.png', dpi=600, bbox_inches='tight')
