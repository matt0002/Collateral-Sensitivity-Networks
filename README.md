# Collateral Sensitivity Network Tool

# Introduction to CSNA:

---

CSNA is a python program designed to take in collateral sensitivity data and return intuitive and helpful graphs representing inter-drug relationships in a bacterial population. This application’s development was motivated by a lack of tools available that leverage collateral sensitivity i.e. resistance to one drug renders susceptibility to another drug. The outcome of this application leads to the optimal drug cycling therapy to effectively treat a bacterial infection. 

To utilize this tool, a matrix must be provided that represents the change in a bacterial population’s susceptibility profile to one drug when exposed to another drug. There are 3 possible inter-drug relationships and they are not always reciprocal. 

- Neutral:
    - No change in Drug A susceptibility when treated with Drug B.
- Collateral Sensitivity:
    - A population is more susceptible to Drug A after treatment with Drug B.
- Collateral Resistance:
    - A population is less susceptible to Drug A after treatment with Drug B.

This program consists of four major functionalities, a principal component analysis (PCA) to show drug profile similarities, a ternary diagram, an optimization function, and a sankey diagram to represent the flow of antibiotic resistance profiles 

## Example Data format:

---

|  | MER | TOB |
| --- | --- | --- |
| MER | 64 | -6 |
| TOB | 0 | 32 |

# Dependencies:

---

To run this program you will need to have python version ≥ 3.11. Python is a free to use programming language found **[here](https://www.python.org/).** 

### On Windows (tested in Windows 11)

1. 
    
    ```powershell
    mkdir CSNA
    cd CSNA
    ```
    

1. 
    
    ```
    py -3 -m pip install numpy pandas openpyxl xlsxwriter plotly dash 
    ```
    
2. Download the repository as a zip file and save the unzipped file in the new directory CSNA
3. To run the code:
    
    ```powershell
    py dash_main.py
    ```
