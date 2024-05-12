import pandas as pd


def matrix_profile(final_graphs, drug1, drug2, drug3, drug_4, mic_cutoff):
    run = 0
    dict_of_phenotypes = {}
    df = final_graphs['DataFrame']
    cols = list(df.columns)
    print(cols)
    drugs = [drug1, drug2, drug3, drug_4]
    drugs_final = []
    [drugs_final.append(x) for x in drugs if x is not None]
    drug_to_alphabet = ['A', 'B', 'C', 'D']
    done_keys = []
    start_pheno = [f'{drugs_final[i]}_S' for i, e in enumerate(drugs_final)]
    total_R = [f'{drugs_final[i]}_R' for i, e in enumerate(drugs_final)]
    print(start_pheno)
    for i, drug in enumerate(drugs_final):
        if drug is None:
            pass
        else:
            print(drug)
            drug_i_index = cols.index(drug)
            list_dr_profiles = []
            for drug_second in drugs_final:
                drug_second_index = cols.index(drug_second)
                profile = df.iloc[drug_i_index, drug_second_index]
                if profile < 0:
                    list_dr_profiles.append(f'{drug_second}_S')
                elif profile == 0:
                    list_dr_profiles.append(f'{drug_second}_0')
                else:
                    list_dr_profiles.append(f'{drug_second}_R')
            dict_of_phenotypes[drug_to_alphabet[i]] = list_dr_profiles

    print(dict_of_phenotypes)

    return dict_of_phenotypes, start_pheno, total_R

