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


def plotFormulas(ax, formulas, title, typeCar=None, paneColor='whitesmoke', edgeColor='dimgrey'):
    ax.set_proj_type('ortho')
    ax.view_init(elev=29, azim=134, roll=0)
    ax.grid(False)
    for axis in ax.xaxis, ax.yaxis, ax.zaxis:
        axis.set_label_position('lower'), axis.set_ticks_position('lower')
        axis.set_pane_color(paneColor), axis.pane.set_edgecolor(edgeColor)

    ax.set_title(title)

    # x axis
    ax.set_xlabel('')
    ax.set_xlim([3, -1])
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(['No Car', 'kCar', 'iCar'])
    # ax.invert_xaxis()

    # y axis
    ax.set_ylabel('Starch')
    ax.set_ylim([12, -2] if typeCar else [2, 13])
    ax.set_yticks([0, 5.0, 10.0] if typeCar else [5.0, 10.0])
    ax.set_yticklabels(['0%', '5%', '10%'] if typeCar else ['5%', '10%'])
    # ax.invert_yaxis()

    # z axis
    ax.set_zlabel(f'Calcium ions concentration (mM)')
    ax.set_zlim([0, 30])
    ax.set_zticks([0, 7, 14, 28])
    ax.set_zticklabels(['0', '7', '14', '28'])

    for formulation in formulas:
        car, starch, cacl2 = formulation
        color = 'whitesmoke'
        markerFill = 'left'
        ax.plot(car, starch, cacl2,
                c=color, alpha=(0.4 + starch * 0.06), fillstyle=markerFill, ms=12, marker='o', mew=0.5, mec='k',
                label='', zorder=2)

        ax.plot([car, car], [starch, starch], [0, cacl2 - 0.06],
                color='#383838', linestyle='-', alpha=(0.4 + starch * 0.06), lw=0.6,
                zorder=1)


formulations = [
    [0, 10, 0], [0, 10, 7], [0, 10, 14], [0, 10, 28],  # St
    [1, 10, 0], [1, 10, 7], [1, 10, 14], [1, 10, 28],  # St + kappaCar
    [2, 10, 0], [2, 10, 7], [2, 10, 14], [2, 10, 28],  # St + iotaCar
]

plt.style.use('seaborn-v0_8-ticks')
fig = plt.figure(figsize=(8, 7), constrained_layout=False)
fig.suptitle('Hydrogels formulations. 40 in total')

ax1 = fig.add_subplot(111, projection='3d')

fig.subplots_adjust(
    left=0.01, right=0.99,
    top=0.87, bottom=0.02,
    wspace=0.05, hspace=0.
)

plotFormulas(
    ax1, formulations,
    'Starch-only | 12 formulations.',
    paneColor='whitesmoke', edgeColor='gainsboro')

# plotFormulas(
#     ax2, kCar,
#     'Kappa Carrageenan | 14 formulations.', 'Kappa Carrageenan',
#     paneColor='floralwhite', edgeColor='antiquewhite')
# plotFormulas(
#     ax3, iCar,
#     'Iota Carrageenan | 14 formulations.', 'Iota Carrageenan',
#     paneColor='azure', edgeColor='paleturquoise')

filename = 'Hydrogels-Formulations-Plot-NEW'
plt.show()
fig.savefig(f'{filename}.png', dpi=600, bbox_inches='tight')
