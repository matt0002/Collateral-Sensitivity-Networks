import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.offline as ply


def build_network(start_pheno, dict_of_phenotypes, total_R):
    drugs = list(dict_of_phenotypes.keys())
    all_vals = [dict_of_phenotypes[key] for key in dict_of_phenotypes]
    arr_post_exposure = np.array(all_vals)
    print('Array')
    print(arr_post_exposure)

    class nodes:
        def __init__(self, node, parent):
            self.node = node
            self.parent = parent
            self.drugs_seen = []
            self.children = []
            self.insensitive = {}

    all_combinations = []

    for after_drug_exposure in range(arr_post_exposure.shape[1]):  # make initial phenotypes
        node_post_exposure = []
        node_working = list(arr_post_exposure[after_drug_exposure, :])
        # drugs in column col
        insensitivities = {}  # drug exposed : drug insensitive
        for i, e in enumerate(node_working):
            if '_0' in e:
                node_post_exposure.append(start_pheno[i])
                insensitivities[after_drug_exposure] = i
            else:
                node_post_exposure.append(e)

        new_node = nodes(node=node_post_exposure, parent=start_pheno)
        new_node.drugs_seen.append(after_drug_exposure)
        new_node.insensitive = insensitivities
        all_combinations.append(new_node)

    all_combinations_secondary = []
    for times_repeat in range(0, len(dict_of_phenotypes.keys()) - 1):  # iterate through the drugs
        for index, phenotype_pre_exposure in enumerate(all_combinations):  # drugs after initial exposure
            for drug_being_used in range(arr_post_exposure.shape[0]):
                if phenotype_pre_exposure.node == total_R:
                    pass
                else:
                    if drug_being_used in phenotype_pre_exposure.drugs_seen:
                        pass
                    else:
                        node_post_exposure = list(arr_post_exposure[drug_being_used, :])

                        phenotype_being_exposed = phenotype_pre_exposure.node

                        new_node = []
                        for state_index, state_element in enumerate(phenotype_being_exposed):
                            if '_0' in node_post_exposure[state_index]:
                                new_node.append(state_element)
                            else:
                                new_node.append(node_post_exposure[state_index])

                        node_post = nodes(node=new_node, parent=phenotype_being_exposed)
                        phenotype_pre_exposure.children.append(node_post)
                        phenotype_pre_exposure.drugs_seen.append(drug_being_used)
                        node_post.drugs_seen.append(drug_being_used)
                        all_combinations_secondary.append(node_post)

    all_combinations_round_2 = all_combinations + all_combinations_secondary

    for times_repeat in range(0, len(dict_of_phenotypes.keys()) - 1):  # iterate through the drugs
        for index, phenotype_pre_exposure in enumerate(all_combinations_round_2):  # drugs after initial exposure
            for drug_being_used in range(arr_post_exposure.shape[1]):
                if phenotype_pre_exposure.node == total_R:
                    pass
                else:
                    if drug_being_used in phenotype_pre_exposure.drugs_seen:
                        pass
                    else:
                        node_post_exposure = list(arr_post_exposure[drug_being_used, :])

                        phenotype_being_exposed = phenotype_pre_exposure.node
                        insensitive = phenotype_pre_exposure.insensitive
                        new_node = []
                        for state_index, state_element in enumerate(phenotype_being_exposed):
                            if '_0' in node_post_exposure[state_index]:
                                new_node.append(state_element)
                            else:
                                new_node.append(node_post_exposure[state_index])

                        node_post = nodes(node=new_node, parent=phenotype_being_exposed)
                        phenotype_pre_exposure.children.append(node_post)
                        node_post.drugs_seen.append(drug_being_used)
                        phenotype_pre_exposure.drugs_seen.append(drug_being_used)
                        all_combinations_secondary.append(node_post)

    df_ = pd.DataFrame()
    from_ = []
    to_ = []
    by_ = []

    for i, e in enumerate(all_combinations_round_2):
        if e.parent == e.node:
            pass
        else:
            from_.append(e.parent)
            to_.append(e.node)

            print('Parent')
            print(e.node)

            print(e.drugs_seen)
            by_.append(e.drugs_seen[0])
            print('popped', e.drugs_seen)
            for ind, c in enumerate(e.children):
                if c.node == e.node:
                    pass
                else:
                    from_.append(e.node)
                    print(c.node)
                    print(c.drugs_seen)
                    to_.append(c.node)
                    by_.append(e.drugs_seen[ind + 1])

    df_['From'] = from_
    df_['To'] = to_
    df_['by'] = by_
    df_ = df_[~df_.astype(str).duplicated()]
    fig_ = make_sankey(df_, dict_of_phenotypes, start_pheno)
    return fig_


def make_sankey(df, dict_of_phenotypes, start_pheno):
    from_ = list(df.loc[:, 'From'])
    to_ = list(df.loc[:, 'To'])
    by_ = list(df.loc[:, 'by'])
    color_by = ['Black', 'Blue', 'Red', 'Green']
    all_values = from_ + to_
    source = []
    target_ = []
    color_ = []
    for i, c in zip(from_, by_):
        print(i)
        source.append(from_.index(i))
        color_.append(color_by[c])
    for t in to_:
        if t in from_:
            target_.append(from_.index(t))
        else:
            target_.append(to_.index(t) + len(source))

    values = [1 for x in range(len(all_values))]
    print(all_values)
    new_values = []
    for lab in all_values:
        to_append = []
        for names in lab:
            to_append.append(f'{names[0]}{names[-1]}')
        new_values.append(to_append)
    sankey = go.Sankey(
        valueformat=".0f",
        node={
            'pad': 15,
            'thickness': 10,
            'label': new_values
        },
        link={
            'arrowlen': 30,
            'source': source,
            'target': target_,
            'value': values,
            'color': color_
        }

    )

    fig = go.Figure(data=sankey)
    title_colors = ['black', 'blue', 'red', 'green']
    final_title = ''
    for i, e in enumerate(dict_of_phenotypes['A']):
        final_title += f"<span style='color:{title_colors[i]}'> {e[:-2]}({e[0]})</span> & "

    fig.update_layout(
        hovermode='x',
        title_text=f"Flow diagram of Antibiotic States for: {final_title[:-2]}",
        font_size=15
    )
    return fig
