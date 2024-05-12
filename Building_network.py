import os.path
import plotly.graph_objects as go
from itertools import permutations
import matplotlib
import numpy as np
import networkx as nx
import pandas as pd
import my_network as my_nx
from networkx.drawing.nx_pydot import graphviz_layout

from dash import Dash, html, dcc

# need to break up dict of phenotypes by individual drug

color_link = ['#000000', '#FFFF00', '#1CE6FF', '#FF34FF', '#FF4A46',
              '#008941', '#006FA6', '#A30059', '#FFDBE5', '#7A4900',
              '#0000A6', '#63FFAC', '#B79762', '#004D43', '#8FB0FF',
              '#997D87', '#5A0007', '#809693', '#FEFFE6', '#1B4400',
              '#4FC601', '#3B5DFF', '#4A3B53', '#FF2F80', '#61615A',
              '#BA0900', '#6B7900', '#00C2A0', '#FFAA92', '#FF90C9',
              '#B903AA', '#D16100', '#DDEFFF', '#000035', '#7B4F4B',
              '#A1C299', '#300018', '#0AA6D8', '#013349', '#00846F',
              '#372101', '#FFB500', '#C2FFED', '#A079BF', '#CC0744',
              '#C0B9B2', '#C2FF99', '#001E09', '#00489C', '#6F0062',
              '#0CBD66', '#EEC3FF', '#456D75', '#B77B68', '#7A87A1',
              '#788D66', '#885578', '#FAD09F', '#FF8A9A', '#D157A0',
              '#BEC459', '#456648', '#0086ED', '#886F4C', '#34362D',
              '#B4A8BD', '#00A6AA', '#452C2C', '#636375', '#A3C8C9',
              '#FF913F', '#938A81', '#575329', '#00FECF', '#B05B6F',
              '#8CD0FF', '#3B9700', '#04F757', '#C8A1A1', '#1E6E00',
              '#7900D7', '#A77500', '#6367A9', '#A05837', '#6B002C',
              '#772600', '#D790FF', '#9B9700', '#549E79', '#FFF69F',
              '#201625', '#72418F', '#BC23FF', '#99ADC0', '#3A2465',
              '#922329', '#5B4534', '#FDE8DC', '#404E55', '#0089A3',
              '#CB7E98', '#A4E804', '#324E72', '#6A3A4C'
              ]

color_dict = {'A': '#575329', 'B': '#7A4900', 'C': '#0CBD66', 'D': '#1CE6FF'}


class Tree:
    def __init__(self, node, parent):
        self.node = node
        self.parent = parent
        self.children = []
        self.drugs_seen = []
        self.children_values = []
        self.drug_to_get_pheno = None


def build_network_collateral_sensitivity(start_pheno, dict_of_phenotypes_after_exposure):
    root = Tree(start_pheno, None)
    root.drugs_seen = [None]
    inst_i = [root]
    for drug in dict_of_phenotypes_after_exposure:
        inst_ = Tree(node=dict_of_phenotypes_after_exposure[drug], parent=root.node)
        root.children.append(inst_)
        root.children_values.append(inst_.node)
        inst_.drugs_seen.append(drug)
        inst_.drug_to_get_pheno = drug
        inst_i.append(inst_)

    for drug_index, drug in enumerate(dict_of_phenotypes_after_exposure):
        after_treatment_phenotype = dict_of_phenotypes_after_exposure[drug]
        for i in inst_i:
            if drug in i.drugs_seen:
                pass
            else:
                i.drugs_seen.append(drug)
                new_pheno = []
                for phenotype_index, phenotype_element in enumerate(i.node):
                    if phenotype_index == drug_index:  # assuming they are in the same position key / value index
                        new_pheno.append(after_treatment_phenotype[drug_index])  # take abx pheno
                    else:
                        new_pheno.append(after_treatment_phenotype[phenotype_index])

                new_inst = Tree(node=new_pheno, parent=i.node)
                i.children.append(new_inst)
                i.children_values.append(new_pheno)
                new_inst.drugs_seen = [drug]
                new_inst.drug_to_get_pheno = drug
                new_inst.parent = i.node
                inst_i.append(new_inst)

    return inst_i


def run_network(start_pheno, dict_of_phenotypes_after_exposure):
    from_ = []
    to_ = []
    by_ = []
    net = build_network_collateral_sensitivity(start_pheno, dict_of_phenotypes_after_exposure)
    for i in net:
        for ind, c in enumerate(i.children_values):
            if i.node == c:
                pass
            else:
                node_concat = ''
                to_concat = ''

                for l in range(len(i.node)):
                    node_concat += i.node[l]
                    to_concat += c[l]

                from_.append(node_concat)
                to_.append(to_concat)
                by_.append(i.children[ind].drug_to_get_pheno)

    df_out = pd.DataFrame({'Source': from_, 'Target': to_, 'by': by_})
    print(df_out)
    source__ = list(df_out.loc[:, 'Source'])
    target = list(df_out.loc[:, 'Target'])
    unique_nodes = []
    [unique_nodes.append(t) for t in target if t not in unique_nodes]
    [unique_nodes.append(x) for x in source__ if x not in unique_nodes]
    make_matrix(df_out, unique_nodes, 0.1, -0.15)
    to_dict_for_sanky = {}
    num_id = 0
    for n in unique_nodes:
        to_dict_for_sanky[n] = num_id
        num_id += 1

    source_out = []
    target_out = []
    values = []

    [source_out.append(to_dict_for_sanky[x]) for x in source__]
    [target_out.append(to_dict_for_sanky[x]) for x in target]
    key_list = list(to_dict_for_sanky.keys())
    values_list = list(to_dict_for_sanky.values())

    pairs = []
    unique_source = []
    unique_target = []
    color_ind = ['A', 'B', 'C', 'D']
    drug_order = []
    color_ = []
    for s, t in zip(source_out, target_out):
        pair = [s, t]
        if pair in pairs:
            print('Cycle')
        else:
            pairs.append(pair)
            unique_source.append(s)
            unique_target.append(t)
            values.append(1)
            position = values_list.index(s)
            position_t = values_list.index(t)
            s_key = key_list[position]
            t_key = key_list[position_t]
            drug_ = df_out.loc[df_out['Source'] == s_key]
            print(drug_)
            drug_ = list(drug_.loc[drug_['Target'] == t_key, 'by'])
            print(drug_[-1])
            color_.append(color_dict[drug_[-1]])
            drug_order.append(drug_[-1])

    print(target_out)
    print(source_out)
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=unique_nodes,
            color="blue",
            customdata=unique_nodes,
            hovertemplate='%{source.customdata}'
        ),
        link=dict(
            source=unique_source,  # indices correspond to labels, eg A1, A2, A1, B1, ...
            target=unique_target,
            value=values,
            color=color_,
            customdata=drug_order,
            hovertemplate='Exposure to %{customdata}'
        ))])
    fig.update_layout(title_text="Basic Sankey Diagram", font_size=10)
    return fig


def make_matrix(df, unique_nodes, alpha, delta):
    # todo: alpha, delta, mu, tw_min, tw_max, t_a, t_b, t_c
    all_drug_matrix = []
    all_drugs = list(df.loc[:, 'by'])
    unique_drugs = []
    [unique_drugs.append(x) for x in all_drugs if x not in unique_drugs]
    df = df.drop_duplicates()
    converted_states = []  
    for state_ in unique_nodes:
        new_rep = []
        for n in state_:
            if n == 'S':
                new_rep.append(delta)
            elif n == 'R':
                new_rep.append(alpha)
            else:
                pass
        converted_states.append(new_rep)
    print(unique_nodes)
    print(converted_states)
    # matrix of growth and clearance
    all_g_c_matrix = []
    for d_ in range(len(unique_drugs)):
        g_c_matrix = []
        for convert_states in converted_states:
            g_c_matrix.append(convert_states[d_])
        diag_ = np.diag(g_c_matrix)
        all_g_c_matrix.append(diag_)

    # autonomous linear system


    # start mutation matrix generation
    for drug in unique_drugs:
        matrix = np.zeros(shape=(len(unique_nodes), len(unique_nodes)))
        df_cut_drug = df.loc[df['by'] == drug]
        for i, state in enumerate(unique_nodes):
            change_by_drug = df_cut_drug.loc[df_cut_drug['Source'] == state]
            if change_by_drug.empty:
                pass
            else:
                state_changed = list(change_by_drug.loc[:, 'Target'])
                val_column_wise = unique_nodes.index(state_changed[0])
                val_row_wise = unique_nodes.index(state)
                matrix[val_column_wise, val_row_wise] = 1

        all_drug_matrix.append(matrix)

    # all drug combinations
    treatment_options = list(permutations(unique_drugs))

    # todo: need to ask ale about this
    # tw = np.ones()

