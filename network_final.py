import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.offline as ply



def build_network(start_pheno, dict_of_phenotypes, total_R):
    drugs = list(dict_of_phenotypes.keys())
    all_vals = [dict_of_phenotypes[key] for key in dict_of_phenotypes]
    arr_ = np.array(all_vals)

    class nodes:
        def __init__(self, node, parent):
            self.node = node
            self.parent = parent
            self.drugs_seen = []
            self.children = []

    all_combinations = []

    for col in range(arr_.shape[1]):  # make initial phenotypes
        new_node = nodes(node=list(arr_[:, col]), parent=start_pheno)
        new_node.drugs_seen.append(col)
        all_combinations.append(new_node)

    for i in range(0, len(dict_of_phenotypes.keys())):
        for comb in all_combinations:
            print(comb.node)
            for col in range(arr_.shape[1]):
                if col in comb.drugs_seen:
                    pass
                else:
                    # need if statement to say you can't leave total R
                    node_ = list(arr_[:, col])
                    if comb.node == total_R:
                        #comb.children.append(list(arr_[:, col]))
                        comb.drugs_seen.append(col)
                    else:
                        comb.children.append(list(arr_[:, col]))
                        comb.drugs_seen.append(col)

    df_ = pd.DataFrame()
    from_ = []
    to_ = []
    by_ = []
    for i in all_combinations:
        from_.append(i.parent)
        to_.append(i.node)
        by_.append(i.drugs_seen[0])
        for ind, c in enumerate(i.children):
            from_.append(i.node)
            to_.append(c)
            by_.append(i.drugs_seen[ind + 1])

    df_['From'] = from_
    df_['To'] = to_
    df_['by'] = by_

    fig_ = make_sankey(df_)
    return fig_

def make_sankey(df):
    from_ = [f'from {x}' for x in list(df.loc[:, 'From'])]
    to_ = [f'to {x}' for x in list(df.loc[:, 'To'])]
    newdf = pd.DataFrame()
    newdf['From'] = from_
    newdf['To'] = to_
    newdf['by'] = df.loc[:, 'by']
    print(newdf)
    counts = newdf.groupby(['From', 'To', 'by']).size().reset_index(name='counts')
    print(counts)
    past = counts['From'].tolist()
    current = counts['To'].tolist()
    unique_values = pd.unique(pd.concat((counts['From'], counts['To']), axis=0)).tolist()
    past_indices = [unique_values.index(i) for i in past]
    current_indices = [unique_values.index(i) for i in current]
    data_trace = dict(
        type='sankey',
        # width = 1118,
        # height = 772,
        node=dict(
            pad = 15,
            thickness = 15,
            line = dict(
            color = "black",
            width = 0.5
               ),
            label=unique_values,
        ),
        link=dict(
            source=past_indices,
            target=current_indices,
            value=counts['counts']
        ))

    layout = dict(
        height=1500,
        # width=80,
        title="Past and current job titles of data professionals",
        font=dict(
            size=20
        )
    )

    fig = go.Figure(data=[go.Sankey(data_trace)], layout=layout)
    fig.write_html('sankey.html')
    return fig


