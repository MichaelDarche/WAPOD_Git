*************************************************
On Windows/VS Code :

        py -m venv venv # Create an environnement
        
        venv\Scripts\Activate.ps1 # Activate the environnement*
        
        py -m pip install --upgrade pip setuptools wheel # Update pip
        
        pip install -r requirements.txt # Install the needed dependencies
        
    * If there is an issue with the rights of execution, run
        Set-ExecutionPolicy -Scope CurrentUser -ExecutionPolicy RemoteSigned
*************************************************
For every plots, go into the folder .\Fig#\ an run "py Fig#.py"
*************************************************
# Input parameters

The simulation is configured using the following files:

* Demarrer.txt : spatial domain and discretisation
* Milieu.txt : material properties (mass density and wave speed) of the media
* Frontiere.txt : interface properties (position, modulated compliance and inertia)
* Source.txt : source parameters

# Parameters for the file name

fs: source frequency

fm: modulation frequencies

# Raw datas

* tp+fs+fm : Time vector

Microstructured medium (velocity–stress formulation)

* Vfilm+fs+fm : velocity fields obtained with velocity-stress formulation on the space domain for every computed time steps 
* Sfilm+fs+fm : stress fields obtained with velocity-stress formulation on the space domain for every computed time steps 

Homogenized medium

* Uhomo+fs+fm : displacement fields obtained with displacement formulation on homogenised media
* VfilmSSH+fs+fm : velocity fields obtained with velocity-stress formulation on homogenised media (Strangs Splitting) on the space domain for every computed time steps 
* SfilmSSH+fs+fm : stress fields obtained with velocity-stress formulation on homogenised media (Strangs Splitting) on the space domain for every computed time steps 

# Plots

- Fig2 - Leading order homogenisation

        CS : Jump values and mean values at the interface
        
        0, Low, Mid, High : Cases of dissipative interfaces

- Fig3 - Leading order homogenisation - Impedance matching

- Fig4 - Leading order homogenisation - Case Ni=2

- Fig5 - Order-2 homogenization without modulation

        Uz_K-1e9_+fs+fm : Leading order displacement
        
        U_K-1e9_+fs+fm : Order-2 displacement
        
        P1, P2 : Order-2 correctors 

- Fig6 - Order-2 homogenization with increasing modulation frequency

        Uz_K-1e9_+fs+fm : Leading order displacement
        
        U_K-1e9_+fs+fm : Order-2 displacement
        
        P1, P2 : Order-2 correctors 

- Fig7 - Non reciprocal features

