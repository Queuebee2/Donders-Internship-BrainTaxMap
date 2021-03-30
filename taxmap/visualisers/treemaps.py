import os

import numpy as np
import pandas
import plotly.express as px
from plotly.express import data
from taxmap import *

TREEMAP_OUTPUT_PATH = str(Path('./data/output/graphs/treemaps/').resolve())


def save_all_figs(fig, figname, directory=TREEMAP_OUTPUT_PATH):
    """take a plotly figure and save it in all available formats under a specific name
    without specific formatting"""
    fig.write_json(
        os.path.join(TREEMAP_OUTPUT_PATH, f'fig-{figname}-json.json'))

    # html is quite big!
    fig.write_html(
        os.path.join(TREEMAP_OUTPUT_PATH, f'fig-{figname}-html.html'))

    # idea : dynamically create list as ENABLED_EXTENTIONS from some config.
    for ext in ['jpeg', 'png', 'webp', 'pdf', 'svg']:
        fig.write_image(os.path.join(directory, f'fig-{figname}.{ext}'))

    # fig.write_image(os.path.join(directory,f'fig-{figname}.jpeg'))
    # fig.write_image(os.path.join(directory,f'fig-{figname}.png'))
    # fig.write_image(os.path.join(directory,f'fig-{figname}.webp'))
    # fig.write_image(os.path.join(directory,f'fig-{figname}.pdf'))
    # fig.write_image(os.path.join(directory,f'fig-{figname}.svg'))


if __name__ == '__main__':
    path1 = str(Path('./data/output/hits/trimmed-annotation.txv').resolve())

    df = pandas.read_csv(path1, sep="\t")
    df["total"] = "total"  # in order to have a single root node
    fig = px.treemap(df,
                     path=['total', 't1', 't2', 't3', 't4'],
                     values='amt',
                     hover_data=['verb', 'pmid', 'struct'],
                     color_continuous_scale='RdBu')

    save_all_figs(fig, 'annotated hits')

    df = df[~df.stack().str.contains('unidentified').any(level=0)]
    fig = px.treemap(df,
                     path=['total', 't1', 't2', 't3', 't4'],
                     values='amt',
                     hover_data=['verb', 'pmid', 'struct'],
                     color_continuous_scale='RdBu')

    save_all_figs(fig, 'annotated-stripped-hits')

    path2 = str(
        Path('./data/output/hits/treemap-overview - Sheet1.tsv').resolve())

    df = pandas.read_csv(path2, sep="\t")
    df.dropna()
    df["total"] = "total"  # in order to have a single root node

    fig = px.treemap(df,
                     path=[
                         'Total', 'Database', 'found-status',
                         'analysis-status', 'filtered', 'hit-type'
                     ],
                     values='count',
                     hover_data=['count'],
                     color_continuous_scale='RdBu',
                     color_continuous_midpoint=np.average(df['count']))

    #color_continuous_midpoint=np.average(df['lifeExp'], weights=df['pop']))

    save_all_figs(fig, 'overview')
