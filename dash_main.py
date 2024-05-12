from dash import Dash
from dash import dcc, html, Input, Output, State, exceptions
import plotly.graph_objects as go
import io
import base64
from network_final import build_network
from cs_matrix_PCA import pca_on_data
import pandas as pd
from optimization import optimization_
from making_matrix_profiles import matrix_profile
from test import ternary
import numpy as np
from heatmap_for_tool import construct_heatmap

alph = ['A', 'B', 'C']
final_graphs = {}
app = Dash(__name__)
app.layout = (html.Div(children=[
    html.H1(children=['Collateral Sensitivity']),
    html.Div(id='Enter Dataset Div',
             children=[html.Div(className='input_file_div',
                                id='input_file_div',
                                children=[dcc.Upload(id='input_file',
                                                     children='Enter your .csv or .xlsx file here by clicking or '
                                                              'dropping: ',
                                                     style={
                                                         'width': '100%',
                                                         'height': '60px',
                                                         'lineHeight': '60px',
                                                         'borderWidth': '1px',
                                                         'borderStyle': 'dashed',
                                                         'borderRadius': '5px',
                                                         'textAlign': 'center',
                                                         'margin': '10px'})])]),
    html.Div(id='pca_div'),
    html.Div(id='ternary'),
    html.Div(id='Heatmap'),
    html.Div(id='optimal_div', children=html.Button('Optimal?', id='optimal_run', n_clicks=0)),
    html.Div(id='drugs_for_optimization', children=dcc.Input(id='drug_num', placeholder='Number of Drugs to Cycle')),
    html.Div(id='div_tern_optimal'),

    html.Div([
        dcc.Dropdown(id='input_phenotypes_1', placeholder="Enter Phenotype after exposure to drug 1"),
        dcc.Dropdown(id='input_phenotypes_2', placeholder="Enter Phenotype after exposure to drug 2"),
        dcc.Dropdown(id='input_phenotypes_3', placeholder="Enter Phenotype after exposure to drug 3"),
        dcc.Dropdown(id='input_phenotypes_4', placeholder="Enter Phenotype after exposure to drug 4"),
        dcc.Input(id='mic_threshold', placeholder='Enter the MIC Threshold to consider resistant'),
        html.P(),
        html.Button('Run', id='button_to_run', n_clicks=0)]),

    html.Div(id='div_for_sankey'),

]))


@app.callback(
    Output(component_id='input_phenotypes_1', component_property='options'),
    Output(component_id='input_phenotypes_2', component_property='options'),
    Output(component_id='input_phenotypes_3', component_property='options'),
    Output(component_id='input_phenotypes_4', component_property='options'),
    Output(component_id='input_file', component_property='children'),
    Output(component_id='pca_div', component_property='children'),
    Output(component_id='ternary', component_property='children'),
    Output(component_id='Heatmap', component_property='children'),
    Input(component_id='input_file', component_property='contents'),
    State('input_file', 'filename'))
def get_file_and_info(contents, filename):
    if contents is None:
        raise exceptions.PreventUpdate

    elif contents is not None:
        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)
        try:
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))

        except Exception as e:
            return html.Div([
                'There was an error processing this file.'
            ])

        df = df.drop(columns='Unnamed: 0')
        columns = df.columns
        df = df.to_numpy()
        df_ = np.sign(df)
        df = pd.DataFrame(df_, columns=columns)
        print(df)
        final_graphs['DataFrame'] = df
        fig_ = pca_on_data(df, df.columns)
        fig__ternary = ternary(df, 0, False)
        heat_map = construct_heatmap(df, columns)

        return (
            list(df.columns), list(df.columns), list(df.columns), list(df.columns), filename, dcc.Graph(figure=fig_),
            dcc.Graph(figure=fig__ternary), None)

    else:
        pass


@app.callback(
    Output(component_id='input_phenotypes_1', component_property='value'),
    Output(component_id='input_phenotypes_2', component_property='value'),
    Output(component_id='input_phenotypes_3', component_property='value'),
    Output(component_id='input_phenotypes_4', component_property='value'),
    Output('div_tern_optimal', 'children'),
    Input('optimal_run', 'n_clicks'),
    Input('drug_num', 'value')

)
def run_optimal(optimal, n_drugs):
    if optimal == 0:
        raise exceptions.PreventUpdate
    else:
        df = final_graphs['DataFrame']
        col = df.columns
        drug_names, indecies = optimization_(df_in=df, columns=col, num_drugs=int(n_drugs))
        tern_optimal = ternary(df, len(drug_names), list(indecies))

        return drug_names[0], drug_names[1], drug_names[2],drug_names[3], dcc.Graph(figure=tern_optimal)


@app.callback(
    Output('div_for_sankey', 'children'),
    Input('input_phenotypes_1', 'value'),
    Input('input_phenotypes_2', 'value'),
    Input('input_phenotypes_3', 'value'),
    Input('input_phenotypes_4', 'value'),
    Input('mic_threshold', 'value'),
    Input('button_to_run', 'n_clicks'))
def run_code(*vals):
    drug1, drug2, drug3, drug_4, mic_cutoff, run = vals
    if mic_cutoff is None:
        mic_cutoff = 0
    if run == 0:
        raise exceptions.PreventUpdate
    else:
        print('MIC ', mic_cutoff)
        dict_of_phenotypes, start_pheno, total_R = matrix_profile(final_graphs, drug1, drug2, drug3, drug_4, mic_cutoff)

        sankey = build_network(start_pheno=start_pheno, dict_of_phenotypes=dict_of_phenotypes, total_R=total_R)

        return dcc.Graph(figure=sankey)


if __name__ == '__main__':
    app.run_server(port=5555, debug=True)
