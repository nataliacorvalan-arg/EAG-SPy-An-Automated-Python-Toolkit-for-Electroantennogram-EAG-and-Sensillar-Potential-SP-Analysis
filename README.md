# EAG-SPy-An-Automated-Python-Toolkit-for-Electroantennogram-EAG-and-Sensillar-Potential-SP-Analysis
EAG-SPy is a Python-based analysis toolkit for automated processing, quantification, and visualization of electroantennogram (EAG) recordings and sensillar potentials (SPs) stored in .abf (Axon Binary File) format.   

The script performs:
- Gaussian low-pass filtering in the frequency domain (FFT-based)
- Automatic SP detection using a relative threshold criterion
- Quantification of amplitude, latency, duration, and kinetic parameters
- Half-rise and half-decay time computation (absolute and normalized)
- Area under the curve (AUC) measurements
- Linear regression analysis of depolarization slopes (50% and 80%)
- First derivative (dV/dt) analysis for kinetic characterization
- Automatic generation of publication-ready figures
- Structured export of results in .csv and .txt formats

The software is designed for:

- robust analysis under noisy experimental conditions
- non-technical end users
- direct execution from a standalone .exe (no Python installation required)

System Requirements

- Operating system: Windows
- Input data: one or more files with extension .abf
- No Python, libraries, or additional software required

How to Use

1. Copy all .abf files you want to analyze into the EAG_SPy_v1_0_0.dist/ folder before running the program.  
2. Do not move the .exe out of its distribution folder. 
4. Double-click EAG_SPy_v1_0_0.exe.
5. The program will prompt you to enter:
	- Gaussian low-pass filter cutoff frequency
	- Analysis time interval (Important: the start time (in seconds) of the analysis interval must match the onset of the olfactory stimulus to ensure accurate latency calculations.
	- Time points t1 and t2 for voltage evaluation
6. The analysis runs automatically.
7. Results and figures are saved automatically in the same folder.
8. When finished, the console window remains open. Press Enter to exit

Output Files

Numerical results:
- sp_results.csv
- sp_results.txt
  
Each row corresponds to one sweep and includes:
- absolute SP minimum
- ΔSP (latency-based amplitude)
- latency
- duration
- area under the curve (AUC)
- rise and decay half-times
- linear regression slopes (50% and 80%)
- minimum derivative (dV/dt)
- temporal offsets between kinetic markers

Automatically generated folders:

sp_figures/
Per-sweep figures showing:
- raw and filtered signals
- latency point
- SP minimum
- regression lines

fitted_curves/
- .txt files containing fitted regression curves

derivative_curves/
- dV/dt curves (numerical data)

derivative_figures/
- Figures of voltage derivatives

histogram_figures/
- Summary histograms of key parameters

Error Handling:

If an unexpected error occurs: 
- a message “Unexpected error occurred” is displayed
- the full traceback is printed
- the console window remains open
- press Enter to exit
This allows users to read and report errors easily.

Notes

- Sweeps without a valid response are kept in the output with NaN values
- The program does not stop due to individual invalid sweeps
- The analysis is robust to:
	baseline drift
	electrical noise
	mechanical artifacts

Usage Recommendations

- Ensure .abf files contain valid sweep data (channel 0 used by default).
- Gaussian low-pass cutoff, interval start/end, and t1/t2 must be chosen according to your experimental protocol.
- The start time (in seconds) of the analysis interval must match the onset of the olfactory stimulus to ensure accurate latency calculations.
- Inspect generated figures for quality control

Software Information

The software is openly available via Zenodo at: https://doi.org/10.5281/zenodo.18760115

- Type: scientific electrophysiology analysis tool
- Distribution: standalone executable
- License: according to academic distribution terms
- DOI: automatically embedded in exported result files

If you use this software, please cite:
 Corvalan, N. A. (2026). "Automated Python Toolkit for Electroantennogram (EAG) and
 Sensillar Potential (SP) Analysis". Zenodo. DOI: 10.5281/zenodo.18760115

