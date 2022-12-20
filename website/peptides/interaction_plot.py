import pathlib
from bokeh.embed import file_html
from bokeh.io import output_file
from bokeh.layouts import column
from bokeh.models import ColumnDataSource, Whisker, HoverTool, Div
from bokeh.plotting import figure, show
import bokeh.resources as bkr
# from bokeh.transform import jitter, factor_cmap
import numpy as np
import pandas as pd

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

def suptitle_and_primary_iso(groups: list[list]):
    '''add a supertitle to the plots, to display over all the figures
    also return the accession number of the primary isoform'''
    base_acc_num = groups[0][0]
    for acc_num, _ in groups[1:]:
        if '.' not in acc_num: # primary isoform
            base_acc_num = acc_num
            break
    return ([Div(text=('<h2 style="text-align:center">'
                       f'Interaction plot for {base_acc_num} isoforms</h2>'))],
        base_acc_num)

def histograms(df: pd.DataFrame, bins=20, show=False):
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
    title=f'{base_acc_num} isoforms interaction plot histograms'
    if show:
        # create a file locally and open the plot
        output_file(filename=f'static/peptides/interaction plot/Interaction plot histograms bokeh {base_acc_num}.html', 
            title=title)
        show(fig_rows)
        return
    return file_html(fig_rows, resources=bkr.CDN, title=title)
    


# def points_with_error_bars(df: pd.DataFrame):
#     df.is_cancer = df.is_cancer.map({True: 'cancer', False: 'non-cancer'})
#     groups = list(df.groupby('acc_num').groups.items())
#     figs, base_acc_num = suptitle_and_primary_iso(groups)
#     classes =  ['cancer', 'non-cancer']
#     for acc_num, inds in groups:
#         p = figure(width=450, height=160, title=acc_num, x_range=classes)
#         p.xgrid.grid_line_color = None
#         subdf = df.iloc[inds, :]
#         g = subdf.groupby('is_cancer')
#         # add error bars
#         upper = g.intensity.quantile(0.8)
#         lower = g.intensity.quantile(0.2)
#         err_src = ColumnDataSource(data = {
#             'base': classes, 'upper': upper, 'lower': lower
#         })
#         error = Whisker(base='base', upper='upper', lower='lower', source=err_src,
#             level='annotation', line_width=2)
#         error.upper_head.size = 20
#         error.lower_head.size = 20
#         p.add_layout(error)
#         # add points for each value, with some jitter
#         p.circle(jitter('is_cancer', width=0.18, range=p.x_range),
#             'intensity', source=subdf,
#             alpha=0.5, size=13,
#             line_color='white', color=factor_cmap('is_cancer', 'HighContrast3', classes))
#         p.yaxis.axis_label = 'MS intensity'
#         p.title = acc_num
#         figs.append(p)
#     output_file(filename=f'static/peptides/interaction plot/Interaction plot whisker bokeh {base_acc_num}.html', 
#         title=f'{base_acc_num} isoforms interaction plot whisker')
#     show(column(figs, sizing_mode='stretch_width'))


if __name__ == '__main__':
    import sys
    FDIR = pathlib.Path(__file__).parent/'static'/'peptides'/'isoform abundance cancer vs not'
    acc_num = sys.argv[1]
    fname = FDIR/f'{acc_num}wide.csv'
    with fname.open() as f:
        df = process_csv(f)
        histograms(df, show=True)
        # points_with_error_bars(df)
