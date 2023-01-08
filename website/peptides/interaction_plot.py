import pathlib
from bokeh.embed import file_html
from bokeh.io import output_file
from bokeh.layouts import column
from bokeh.models import (
    ColumnDataSource,
    Whisker,
    HoverTool,
    Div
)
from bokeh.palettes import Category20
from bokeh.plotting import figure, show
import bokeh.resources as bkr
# from bokeh.transform import jitter, factor_cmap
import numpy as np
import pandas as pd

CUR_DIR = pathlib.Path(__file__).parent

def process_csv(buffer):
    df = pd.read_csv(buffer)
    piv = df.melt(id_vars='patient', 
        var_name='accnum_type', value_name='intensity'
    )
    accnum_type = piv.accnum_type.str.split('_', expand=True)
    piv['acc_num'] = accnum_type.iloc[:, 0]
    piv['is_cancer'] = accnum_type.iloc[:, 1] == 'C'
    del piv['accnum_type']
    return piv

def suptitle_and_primary_iso(groups: list):
    '''add a supertitle to the plots, to display over all the figures
    also return the accession number of the primary isoform'''
    base_acc_num = groups[0][0]
    for acc_num, _ in groups[1:]:
        if '.' not in acc_num: # primary isoform
            base_acc_num = acc_num
            break
    if '.' in base_acc_num:
        base_acc_num = base_acc_num[:base_acc_num.index('.')]
    return ([Div(text=('<h2 style="text-align:center">'
                       f'Interaction plot for {base_acc_num} isoforms</h2>'))],
        base_acc_num)

def histograms(df: pd.DataFrame, bins=20, show_figs=False):
    '''
    For each isoform, make two histograms,
    one red for cancer and one blue for non-cancer.
    Also make vertical lines at the mean MS intensities
    for cancer and non-cancer.
    '''
    # take the long data, and get one column for cancer and one for non-cancer
    wide = (pd.merge(df[df.is_cancer], df[~df.is_cancer], on=['patient', 'acc_num'])
        .loc[:, ['acc_num', 'intensity_x', 'intensity_y']]
        .rename({
            'intensity_x': 'cancer', 
            'intensity_y': 'non-cancer'
        }, axis = 1)
    )
    groups = list(wide.groupby('acc_num').groups.items())
    figs, base_acc_num = suptitle_and_primary_iso(groups)
    # add supertitle
    for ii, (acc_num, inds) in enumerate(groups):
        p = figure(width=1000, height=180, title=acc_num)
        subdf = wide.iloc[inds, :]
        # add histograms for cancer and non-cancer
        hist_dict = {typ: np.histogram(subdf[typ], bins=bins) 
                    for typ in ['cancer', 'non-cancer']}
        # this will be the height of the vertical bar for the means
        counts_max = max(hist[0].max() for hist in hist_dict.values())
        for typ, color in zip(['cancer', 'non-cancer'], ['red', 'blue']):
            counts, edges = hist_dict[typ]
            src = ColumnDataSource({
                'count': counts,
                'left': edges[:-1],
                'right': edges[1:],
                'type': [typ] * bins # use this to display the type in tooltips
            })
            # add the histogram
            glyph = p.quad(source=src, top='count', bottom=0, left='left', right='right',
                fill_color=color, line_color=color, alpha=0.5,
                legend_label=typ)
            # add tooltips for each box of the histograms
            glyph_hover = HoverTool(
                renderers=[glyph], # add this tool only for this glyph
                tooltips=[
                ('intensity', '@left-@right'),
                ('count', '@count'),
                ('type', '@type')
            ])
            p.add_tools(glyph_hover)
            # add vertical lines at the mean intensity of both types
            mean = subdf[typ].mean()
            mean_line_ys = [0, counts_max * 1.1]
            p.line(x=[mean, mean], y=mean_line_ys, line_width=4, color=color)
        p.yaxis.axis_label = 'Frequency'
        if ii == len(groups) - 1:
            p.xaxis.axis_label = 'MS intensity'
        p.legend.click_policy = 'hide'
        # 'mute' is also a click_policy option, but we want to be able to see
        # the tooltip for exactly one histogram if two overlap, and 'mute'
        # allows two tooltips to show up at once
        # 'mute' has the advantage of not completely hiding the histogram
        figs.append(p)
    fig_rows = column(figs, sizing_mode='stretch_width')
    title = f'{base_acc_num} isoforms interaction plot histograms'
    if show_figs:
        # create a file locally and open the plot
        output_file(filename=CUR_DIR/f'static/peptides/interaction plot example/Interaction plot histograms bokeh {base_acc_num}.html', 
            title=title)
        show(fig_rows)
        return
    return file_html(fig_rows, resources=bkr.CDN, title=title)
    


def points_with_error_bars(df: pd.DataFrame, show_figs=False):
    '''For each accession number, create a plot
    where cancer and non-cancer each have error bar
    of +/- 1 standard deviation and a big point at mean intensity.'''
    df.is_cancer = df.is_cancer.map({True: 'cancer', False: 'non-cancer'})
    groups = list(df.groupby('acc_num').groups.items())
    figs, base_acc_num = suptitle_and_primary_iso(groups)
    classes =  ['cancer', 'non-cancer']
    title = f'{base_acc_num} isoforms interaction plot whisker'
    # get +/-1 standard deviation data for each isoform
    ymin = float('inf')
    ymax = 0
    uppers = []
    lowers = []
    for _, inds in groups:
        g = df.iloc[inds, :].groupby('is_cancer')
        stddev = g.intensity.std()
        mean_ = g.intensity.mean()
        upper = mean_ + stddev
        lower = mean_ - stddev
        plus_1_std = max(upper)
        minus_1_std = max(lower)
        if minus_1_std < ymin:
            ymin = minus_1_std
        if plus_1_std > ymax:
            ymax = plus_1_std
        lowers.append(lower)
        uppers.append(upper)
    ypad = (ymax - ymin) * 0.1
    ymin -= ypad
    ymax += ypad
    # make a separate plot for each accession number
    p = figure(width=700, height=300, x_range=classes, y_range=(ymin, ymax))
    ngroups = len(groups)
    colors = ['red', 'blue'] if ngroups == 2 else Category20[min(ngroups, 20)]
    for ii, (acc_num, inds), upper, lower in zip(range(ngroups), groups, uppers, lowers):
        color = colors[ii % 20]
        p.xgrid.grid_line_color = None
        subdf = df.iloc[inds, :]
        g = subdf.groupby('is_cancer')
        # add +/-1 standard deviation error bars
        err_src = ColumnDataSource(data = {
            'base': classes, 'upper': upper, 'lower': lower
        })
        error = Whisker(base='base', upper='upper', lower='lower', source=err_src,
            level='annotation', line_width=2, line_color=color)
        error.upper_head.size = 20
        error.upper_head.line_color = color
        error.lower_head.size = 20
        error.lower_head.line_color = color
        p.add_layout(error)
        # add points at the mean intensities for each of cancer and non-cancer
        means = [
            subdf[subdf.is_cancer == 'cancer'].intensity.mean(),
            subdf[subdf.is_cancer == 'non-cancer'].intensity.mean()
        ]
        p.line(x=['cancer', 'non-cancer'], y=means, color=color,
            legend_label=acc_num)
        p.scatter(x=['cancer', 'non-cancer'], y=means, size=8, color=color)
        # add a label and a title indicating the accession number
        p.yaxis.axis_label = 'MS intensity'
        figs.append(p)
    p.add_layout(p.legend[0], 'right')
    fig_rows = column(figs, sizing_mode='stretch_width')
    if show_figs:
        # create a file locally and open the plot
        output_file(filename=CUR_DIR/f'static/peptides/interaction plot example/Interaction plot whisker bokeh {base_acc_num}.html', 
            title=title)
        show(fig_rows)
        return
    return file_html(fig_rows, resources=bkr.CDN, title=title)


if __name__ == '__main__':
    import sys
    FDIR = CUR_DIR/'static'/'peptides'/'isoform abundance cancer vs not'
    acc_num = sys.argv[1]
    fname = FDIR/f'{acc_num}wide.csv'
    with fname.open() as f:
        df = process_csv(f)
        # histograms(df, show_figs=True)
        points_with_error_bars(df, show_figs=True)
