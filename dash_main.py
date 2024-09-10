from dash import Dash
from dash import dcc, html, Input, Output, State, exceptions
import plotly.graph_objects as go
import io
import base64
from network_final import build_network
from cs_matrix_PCA import pca_on_data
import pandas as pd
from optimization_ import optimization_
from making_matrix_profiles import matrix_profile
from test import ternary
import numpy as np
from heatmap_for_tool import construct_heatmap
from Drug_matrix_for_ODE import diff_eq
from differential_evo import set_up_diff_evo

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
        dcc.Dropdown(id='input_phenotypes_1', placeholder="Enter drug 1"),
        dcc.Dropdown(id='input_phenotypes_2', placeholder="Enter drug 2"),
        dcc.Dropdown(id='input_phenotypes_3', placeholder="Enter drug 3"),
        dcc.Dropdown(id='input_phenotypes_4', placeholder="Enter drug 4"),
        dcc.Input(id='mic_threshold', placeholder='Enter the MIC Threshold to consider resistant'),
        html.P(),
        html.Button('Run', id='button_to_run', n_clicks=0)]),

    html.Div(id='div_for_sankey'),
    html.Div(children=[
        html.H4(children='Enter growth rate for antibiotic resistant population'),
        dcc.Input(id='Mutant_growth', value=0.29),
        html.H4(children='Enter growth rate for susceptible bacteria'),
        dcc.Input(id='Susceptible_growth', value=0.21),
        html.H4(children='Enter death rate for bacteria'),
        dcc.Input(id='Death_rate', value=0.25),
        html.H4(children='Enter Carrying Capacity'),
        dcc.Input(id='K', value=10 ** 10),
        html.H4(children='Enter mutation rate'),
        dcc.Input(id='mu_rate', value=10 ** -7),
        html.H4(children='Enter initial population size'),
        dcc.Input(id='x_i', value=10 ** 9),
        html.H4(children='Enter hours to dose with one antibiotic'),
        dcc.Input(id='cycle_time', value=72),
        html.P(),
        html.Button('Run ODE', id='button_to_run_ODE', n_clicks=0),
        html.Button('Run Differential Evolution', id='diff_evo_button', n_clicks=0)
    ], style={'width': '20%', 'display': 'inline-block'}),
    html.Div(),
    html.Div(id='div_for_ODE', style={'width': '50%', 'display': 'inline-block'}),
    html.Div(id='graph_for_diff_evo', style={'width': '50%', 'display': 'inline-block'}),
    html.Div(id='sequential_sequence', style={'width': '50%', 'display': 'inline-block'}),
    html.Div(id='best_sequence', style={'width': '50%', 'display': 'inline-block'})

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
        while len(drug_names) < 4:
            drug_names.append(None)

        return drug_names[0], drug_names[1], drug_names[2], drug_names[3], dcc.Graph(figure=tern_optimal)


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


@app.callback(Output('div_for_ODE', 'children'),
              Output('sequential_sequence', 'children'),
              Input('input_phenotypes_1', 'value'),
              Input('input_phenotypes_2', 'value'),
              Input('input_phenotypes_3', 'value'),
              Input('input_phenotypes_4', 'value'),
              Input('Mutant_growth', 'value'),
              Input('Susceptible_growth', 'value'),
              Input('Death_rate', 'value'),
              Input('cycle_time', 'value'),
              Input('K', 'value'),
              Input('mu_rate', 'value'),
              Input('x_i', 'value'),
              Input('button_to_run_ODE', 'n_clicks'))
def ODE(*vals):
    drug1, drug2, drug3, drug_4, mutant_growth, susceptible_growth, death_rate, cycle_time, K, mu_rate, initial, run = vals
    if run == 0:
        raise exceptions.PreventUpdate
    else:
        ODE_graph, sequence = diff_eq(final_graphs, drug1, drug2, drug3, drug_4, float(mutant_growth), float(susceptible_growth),
                            float(death_rate), int(cycle_time), float(K), float(mu_rate), float(initial))
        return dcc.Graph(figure=ODE_graph), dcc.Graph(figure=sequence)


@app.callback(Output('graph_for_diff_evo', 'children'),
              Output('best_sequence', 'children'),
              Input('input_phenotypes_1', 'value'),
              Input('input_phenotypes_2', 'value'),
              Input('input_phenotypes_3', 'value'),
              Input('input_phenotypes_4', 'value'),
              Input('Mutant_growth', 'value'),
              Input('Susceptible_growth', 'value'),
              Input('Death_rate', 'value'),
              Input('cycle_time', 'value'),
              Input('K', 'value'),
              Input('mu_rate', 'value'),
              Input('x_i', 'value'),
              Input('diff_evo_button', 'n_clicks'))
def diff_evo(*vals):
    drug1, drug2, drug3, drug_4, mutant_growth, susceptible_growth, death_rate, cycle_time, K, mu_rate, initial, run = vals
    if run == 0:
        raise exceptions.PreventUpdate
    else:
        diff_evo, sequence_fig = set_up_diff_evo(final_graphs, drug1, drug2, drug3, drug_4, float(mutant_growth),
                                                 float(susceptible_growth),
                                                 float(death_rate), int(cycle_time), float(K), float(mu_rate),
                                                 float(initial))
        return [dcc.Graph(figure=diff_evo), dcc.Graph(figure=sequence_fig)]


if __name__ == '__main__':
    app.run_server(port=5555, debug=True)
