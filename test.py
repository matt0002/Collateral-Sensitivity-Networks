import plotly.express as px
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import plotly.graph_objects as go
import numpy as np
from itertools import combinations
from optimization import optimization_


def makeAxis(title, tickangle):
    return {
        'title': title,
        'titlefont': {'size': 20},
        'tickangle': tickangle,
        'tickfont': {'size': 15},
        'tickcolor': 'rgba(0,0,0,0)',
        'ticklen': 5,
        'showline': True,
        'showgrid': True
    }


def ternary(df, num_drugs_, index):
    df_cut = df
    columns = df_cut.columns
    if index is not False:
        df_cut = df_cut.iloc[index, index]
        columns = df_cut.columns
    list_of_dict_for_ternary = []
    #           WT, CS, CR

    print(df_cut.shape)
    for row in range(df_cut.shape[0]):
        all_row_vals = list(df_cut.iloc[row, :])
        wt_count = 0
        cr_count = 0
        cs_count = 0
        for r_v in all_row_vals:
            if r_v == 0:
                wt_count += 1
            elif r_v > 0:
                cr_count += 1
            else:
                cs_count += 1

        list_of_dict_for_ternary.append(
            {'WT': wt_count, 'CS': cs_count, 'CR': cr_count, 'label': str(columns[row])})
    print(list_of_dict_for_ternary)
    final_traces = []
    d = [i for i in map(lambda x: x['WT'], list_of_dict_for_ternary)]
    print(d)
    for i in list_of_dict_for_ternary:
        final_traces.append(go.Scatterternary(mode='markers',
                                              a=[i['CR']],
                                              c=[i['CS']],
                                              b=[i['WT']],
                                              name=i['label'],
                                              marker={
                                                  'symbol': 100,
                                                  'color': '#DB7365',
                                                  'size': 14,
                                                  'line': {'width': 2}},
                                              showlegend=True))

    fig = go.Figure(data=final_traces)

    fig.update_layout({
        'ternary': {
            'sum': df_cut.shape[0],
            'aaxis': makeAxis('Cross Resistance', 0),
            'baxis': makeAxis('<br>Insensitive', 45),
            'caxis': makeAxis('<br>Collateral Sensitive', -45)
        },
        'annotations': [{
            'showarrow': False,
            'text': 'Collateral Sensitivity to all drugs',
            'x': 0.5,
            'y': 1.3,
            'font': {'size': 15}
        }]
    })

    return fig
