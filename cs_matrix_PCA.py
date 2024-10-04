import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import plotly.graph_objects as go
import numpy as np


def pca_on_data(df, columns):
    colors = ['aqua',
              'black', 'blue',
              'blueviolet', 'brown', 'burlywood', 'cadetblue',
              'chocolate', 'coral', 'cornflowerblue',
              'crimson', 'cyan', 'darkblue', 'darkcyan',
              'darkgoldenrod', 'darkgray', 'darkgrey', 'darkgreen',
              'darkmagenta', 'darkolivegreen', 'darkorange',
              'magenta', 'maroon', 'mediumaquamarine',
              'mediumblue', 'mediumorchid', 'mediumpurple',
              'mediumseagreen', 'mediumslateblue', 'mediumspringgreen',
              'mediumturquoise', 'mediumvioletred', 'midnightblue',
              'mistyrose', 'moccasin', 'navy',
              'orange', 'orangered']

    labels = columns
    df = StandardScaler().fit_transform(df)
    pca = PCA(n_components=3)
    pca_df = pca.fit_transform(df)
    pca_df = pd.DataFrame(pca_df)
    graph_data = []
    for i, e in enumerate(labels):
        print(pca_df.iloc[i, 0])
        graph_data.append(
            go.Scatter3d(x=np.array(pca_df.iloc[i, 0]), y=np.array(pca_df.iloc[i, 1]), z=np.array(pca_df.iloc[i, 2]),
                         marker_color=colors[i], name=e))
    fig = go.Figure(data=graph_data)
    fig.update_layout(legend_title_text='Antibiotic')
    pc_1 = '{:2f}'.format(pca.explained_variance_ratio_[0] * 100)
    pc_2 = '{:2f}'.format(pca.explained_variance_ratio_[1] * 100)
    pc_3 = '{:2f}'.format(pca.explained_variance_ratio_[2] * 100)
    fig.update_layout(scene=dict(
        xaxis_title=f'PC 1 {pc_1} %',
        yaxis_title=f'PC 2 {pc_2} %',
        zaxis_title=f'PC 3 {pc_3} %'))
    return fig
