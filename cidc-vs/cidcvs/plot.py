#!/usr/bin/env python3

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from collections import OrderedDict
from adjustText import adjust_text
import pandas as pd
#################### heatmap ####################


def mapColor(df, colors=None):
    legend_names = df.columns
    legend_dirs = OrderedDict()
    for legend_title in legend_names:
        _factors = list(df[legend_title].unique())
        _factors.sort()
        if isinstance(colors, dict) and (legend_title in colors):
            sub_colors = colors[legend_title]
        else:
            sub_colors = ['indianred', 'steelblue', 'gold', 'darkgreen']
        # truncate number of colors to fix factors
        sub_colors = sub_colors[:len(_factors)]
        _color_wrapper = dict(zip(_factors, sub_colors))
        df[legend_title] = df[legend_title].map(_color_wrapper)
        legend_dirs[legend_title] = _color_wrapper

    return df, legend_dirs

def genLegend(color_dict,ax):
    legend_list = []
    legend_y = 1
    for legend_title,panel in color_dict.items():
        _legend = [mpatches.Patch(color=value, label=key)
                   for key, value in panel.items()]
        cur_legend = ax.legend(
            loc=(1.1, legend_y), handles=_legend, frameon=False, title=legend_title)
        cur_legend._legend_box.align = 'left'
        legend_list.append(cur_legend)
        legend_y -= .1 * len(panel.keys())

    # Update legend to the figure
    for x in legend_list:
        ax.add_artist(x)


def heatMap(df, cmap='RdBu_r', annot=True, center=.3, fmt=".1f", figsize=(12, 12), 
            col_colors=None, col_color_dict=None, row_colors=None, row_color_dict=None):
    
    if not col_colors is None:
        col_colors, col_legend = mapColor(col_colors, colors=col_color_dict)
    if not row_colors is None:
        row_colors, row_legend = mapColor(row_colors, colors=row_color_dict)

    ax = sns.clustermap(df, cmap=cmap, annot=annot,
                        center=center, fmt=fmt, figsize=figsize,
                        col_colors=col_colors,row_colors=row_colors).ax_heatmap
    # remove x tickes
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xticks([])
    ax.set_yticks([])

    # add legend
    if not col_colors is None:
        genLegend(color_dict=col_legend, ax=ax)
    if not row_colors is None:
        genLegend(color_dict=row_legend, ax=ax)

    return ax


#################### scatter plot #################### 
def label_point(df,x, y, label, ax):
    x_v = df[x].tolist()
    y_v = df[y].tolist()
    label_v = df[label].tolist()
    texts = [ax.text(x_v[i], y_v[i], '%s' % label_v[i],
                     ha='center', va='center') for i in range(len(x_v))]
    adjust_text(texts, arrowprops=dict(arrowstyle='->', color='k'), expand_text=(1.01, 0.5), expand_points=(2, 2),
                force_text=(0.5, 2), force_points=(0.06, 0.25))

def add_identity(axes, *line_args, **line_kwargs):
    identity, = axes.plot([], [], *line_args, **line_kwargs)

    def callback(axes):
        low_x, high_x = axes.get_xlim()
        low_y, high_y = axes.get_ylim()
        low = max(low_x, low_y)
        high = min(high_x, high_y)
        identity.set_data([low, high], [low, high])
    callback(axes)
    axes.callbacks.connect('xlim_changed', callback)
    axes.callbacks.connect('ylim_changed', callback)
    return axes

def qqPlot(df, x, y, factor=None, ax=None, label=None, **kwargs):
    if ax is None:
        ax = plt.gca()
    axis_min = min(df[x].min(), df[y].min())
    axis_max = max(df[x].max(), df[y].max())
    axis_min = axis_min * 0.8 if axis_min > 0 else axis_min*1.2
    axis_max = axis_max * 1.2 if axis_max > 0 else axis_max*0.8

    if factor is None:
        ax.scatter(x=x, y=y, data=df, **kwargs)
    else:
        for f in df[factor].unique().tolist():
            ax.scatter(
                x=x, y=y, data=df[df[factor] == f], label=f, **kwargs)
        ax.legend()

    ax.set_xlim(axis_min, axis_max)
    ax.set_ylim(axis_min, axis_max)
    ax.set(xlabel=x, ylabel=y)
    add_identity(ax, color='r', ls='--')
    if not label is None:
        label_point(x=x, y=y, label=label, df=df, ax=ax)
    
   
#################### dendrogram ####################   
import scipy.spatial.distance as ssd
from scipy.cluster import hierarchy
from scipy.spatial import distance


def dendrogramPlot(df, linkage_method='single', threads=0, label=None, ax=None):

    data_linkage = hierarchy.linkage(ssd.squareform(df), method=linkage_method)
    if ax is None:
        ax = plt.gca()

    hierarchy.dendrogram(data_linkage, labels=df.index, leaf_rotation=0, orientation='left',
                         above_threshold_color='k', color_threshold=threads, distance_sort='descending',ax=ax)

    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    if label is not None:
        xlbls = ax.get_ymajorticklabels()
        for lbl in xlbls:
            lbl.set_color(label[lbl.get_text()])

    ax.set(xlabel=r'$ 1 - Jaccard\ Index$')
    return ax



#################### barplot ####################


def barPlot(df, x, y, hue, hue_order,title,ax=None):
    if ax == None:
        ax = plt.gca()

    sns.barplot(x=x, y=y, hue=hue, data=df, ax=ax,
                palette='pastel', hue_order=hue_order)
    ax.set(xlabel='', ylabel='', title=title)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    return ax


def plot96Mtrx(df, height='Seq', neighbor='Neighbor', mut='Alt', title='', rows=None):
    tri_n = OrderedDict([('A_A', 0), ('A_C', 0), ('A_G', 0), ('A_T', 0),
                         ('C_A', 0), ('C_C', 0), ('C_G', 0), ('C_T', 0),
                         ('G_A', 0), ('G_C', 0), ('G_G', 0), ('G_T', 0),
                         ('T_A', 0), ('T_C', 0), ('T_G', 0), ('T_T', 0),
                         ])

    mut_panel = ['C>A', 'C>G', 'C>T', 'T>C', 'T>G', 'T>A']
    colors = ["deepskyblue", "black", "red",
              "lightgray", 'springgreen', "pink"]
    if not rows is None:
        row_panel = sorted(df[rows].unique())
    else:
        row_panel = [title]
    num_row = len(row_panel)
    fig, axarr = plt.subplots(num_row, 6, sharex=True,
                              sharey=True, figsize=(24, 2*num_row))
    for j, r in enumerate(row_panel):
        for i, c in enumerate(mut_panel):
            tmp_tri = pd.Series(tri_n)
            determiner = (df[mut] == c) 
            if num_row > 1:
                determiner = determiner & (df[rows] == r)
                
            neighbor_h = pd.Series(
                df.loc[determiner, [neighbor, height]].set_index(neighbor).to_dict()[height])
            neighbor_h /= neighbor_h.sum()
            tmp_tri.update(neighbor_h)

            ax = axarr[j, i]
            if j == 0:
                ax.set_title(c, weight='bold')
                ax.spines['top'].set_visible(False)
            if i < len(mut_panel) - 1:
                ax.spines['right'].set_visible(False)

            if i == 0 and num_row == 1:
                ax.set_ylabel('Relative Contribution', weight="bold")
            if i == 0 and num_row > 1:
                ax.set_ylabel(r, weight="bold", fontsize=12)
            tmp_tri.plot(kind='bar', ax=ax, color=colors[i], width=.9)
#     plt.suptitle('Relative Contribution', x=0.1, y=.5,va='center',rotation=90, fontsize = 20,fontweight='bold')
    fig.subplots_adjust(wspace=0, hspace=0.2)
    return fig

        
    
                         
        
