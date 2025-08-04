# Interactive Collateral Sensitivity Platform (ICSP)

# Dependencies:

---

To run this program you will need to have python version ≥ 3.11. Python is a free to use programming language found **[here](https://www.python.org/).** 

### On Windows (tested in Windows 11)

1. Navigate to PowerShell or terminal and make a new directory to store all files:
    
    ```powershell
    mkdir ICSP
    cd ICSP
    ```
    

1. Navigate to GitHub and download the package as a .zip file. Once in your downloads move the folder within the zip file to the new directoryjust made, ICSP.
2. Now back in poweshell, change directory to file directory.
    
    ```powershell
    cd ICSP
    cd Collateral-Sensitivity-Networks-main
    ```
    
3. Install all dependencies using pip:
    
    ```
    py -3 -m pip install numpy scikit-learn pandas openpyxl xlsxwriter plotly dash seaborn
    ```
    
4. To run the code:
    
    ```powershell
    py dash_main.py
    ```
    
5. ctrl+click on the https link that shows in terminal.
   
6. Once done with the code, re-enter powershell and press ctrl+C to end the code.

# Introduction to ICSP:

---

ICSP is a python program designed to take in collateral sensitivity data and return intuitive and helpful graphs representing inter-drug relationships in a bacterial population. This application’s development was motivated by a lack of tools available that leverage collateral sensitivity i.e. resistance to one drug renders susceptibility to another drug. The outcome of this application leads to optimal drug cycling therapy to effectively treat a bacterial infection. 

To utilize this tool, a matrix must be provided that represents the change in a bacterial population’s susceptibility profile to one drug when exposed to another drug. There are 3 possible inter-drug relationships and they are not always reciprocal:

- Insensitive:
    - No change in Drug A susceptibility when treated with Drug B.
- Collateral Sensitivity:
    - A population is more susceptible to Drug A after treatment with Drug B.
- Cross Resistance:
    - A population is less susceptible to Drug A after treatment with Drug B.

## Major Functionalities:

1. Data Visualization
    1. A principal component analysis (PCA), ternary diagram, and heatmap are utilized to show drug susceptibility profiles.  A sankey diagram is generated to represent the flow of antibiotic resistance profiles.
    
2. Optimization of Drugs
    1. An optimization function to show which n-drugs exhibit the strongest collateral sensitivities to each other.
    
3. Simulation
    1. A tunable differential equation with user-adjustable parameters. The user is able to choose the antibiotics which they want to cycle with and the simulation returns the bacterial populations along with all of the potential sub-populations. The results are the product of periodic cycling between drugs in the order they are chosen.
    
4. Antibiotic Cycling Optimization
    1. Using the parameters from the simulation differential evolution will be used to optimize cycling with the given drugs. Through a stochastic process, an optimal antibiotic cycling protocol will be found with an objective function to minimize the population size. Major differences between this function and the simulation is that in the simulation, the order of antibiotics is fixed, and the time for treatment with one antibiotic is fixed. Differential evolution varies the order and time of treatment of each antibiotic to minimize the population.

## Example Data format:

---

Attached is an example data file for reference, however a brief example is shown below. This work considered MIC fold-change as an indicator for drug interactions.

|  | Antibiotic A | Antibiotic B | Antibiotic C | Antibiotic D |
| --- | --- | --- | --- | --- |
| Antibiotic A | 64 | -6 | … | … |
| Antibiotic B | 0 | 32 | … | … |
| Antibiotic C | … | … | … | … |
| Antibiotic D | … | … | … | … |

Strains treated with a previous antibiotic or that are resistant to an antibiotic should be positioned on the vertical axis. On the horizontal axis should be drugs exposed to the respective vertical axis. For example, this table would read for a bacteria resistant to antibiotic A, antibiotic B is more effective.

**Please use a .csv file**

## General Workflow:

---

![alt text](https://github.com/matt0002/Collateral-Sensitivity-Networks/tree/main/assets/Figure1z1.pdf)
