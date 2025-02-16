import plotly.graph_objects as go
import numpy as np


def construct_heatmap(df_, columns):
    print(list(columns))
    all_values = []
    for col in columns:
        all_values.append(list(df_.loc[:, col]))
    print(all_values)

    fig = go.Figure(data=go.Heatmap(z=all_values, x=list(columns), y=list(columns), colorscale='rdbu_r'))

    return fig
