# =========================================================================================================================================
# EAG-SPy: An Automated Python Toolkit for Electroantennogram (EAG) and Sensillar Potential (SP) Analysis".
#
# Author: Natalia Andrea Corvalan
# Affiliation: Universidad Nacional de Córdoba (UNC) – CONICET
# Contact: natalia.corvalan@gmail.com | natalia.corvalan@unc.edu.ar
# Year: 2026
# =========================================================================================================================================
# The latency detection approach is conceptually based on:
#   Pézier et al. (2007)
#   "Ca2+ Stabilizes the Membrane Potential of Moth Olfactory Receptor Neurons
#   at Rest and Is Essential for Their Fast Repolarization."
#   Chemical Senses 32: 305–317.
#   doi:10.1093/chemse/bjl059
#-----------------------------------------------------------------------------------------------------------------------------------------
# What does this script do?

#       - Filtering, automatic detection, and quantification of sensillar potentials (SPs) stored in .abf (Axon Binary File) format;
#       - Individual visualization and statistical summary.

# Summary:
#   ├── Processes all sweeps of each .abf file in the directory.
#   ├── Applies a Gaussian lowpass filter to each signal for frequency filtering: gaussian_lowpass_filter function.
#   ├── Computes:
#           ├── Sensillar Potential (SP), latency, duration, area under the curve, half-rise time (t1/2 rise),
#           ├── Half-decay time (t1/2 decay), potential at two user-defined timepoints t1 (Vm_t1) and t2 (Vm_t2),
#           ├── Percentage of SP drop at those times (drop_t1 and drop_t2, %).
#   ├── Saves results in .csv and .txt files named resultados_SP.csv and resultados_SP.txt, respectively.
#   ├── Displays plots in /figuras_SP/ and /figuras_histogramas/

#-----------------------------------------------------------------------------------------------------------------------------------------
# Required libraries (install via Windows PowerShell): pip install pyabf numpy matplotlib pandas openpyxl

#       pyabf: loads and manages .abf files → pip install pyabf
#       numpy, fft: numerical operations and frequency filtering → numpy.fft.fft, numpy.fft.ifft, numpy.fft.fftfreq
#       matplotlib.pyplot: for plotting
#       pandas: for organizing and exporting results → pip install pandas openpyxl
#-----------------------------------------------------------------------------------------------------------------------------------------
# Includes:

# 1. Gaussian_lowpass_filter function
#       Applies a Gaussian lowpass filter to a signal. Filtering is applied to each full sweep before any temporal cropping or analysis.
#           1.1. Transforms the signal into the frequency domain (FFT).
#           1.2. Applies a Gaussian profile centered at 0 Hz with cutoff defined by cutoff_hz.
#                That is, it applies a Gaussian filter to suppress high-frequency components above the cutoff frequency (e.g= 5).
#           1.3. Returns the smoothed signal via inverse FFT (iFFT): Transforms the filtered signal back to the time domain.

# 2. User input
#       Prompts for: 
#           cutoff: cutoff frequency for the filter (in Hz)
#           start and end: time window for analysis within each sweep (only the segment in [start, end] is used)
#           t1 and t2: specific timepoints to evaluate voltage and compute % drop (within [SP, end]) or rise (within [start, SP])
#           # Note: The selected start time for the analysis interval must coincide with the start of the olfactory stimulus 
#           to ensure correct latency calculation.

# 3. Parameters computed per sweep:
#       - SP (mV): minimum voltage within the analysis window [start, end]
#       - Latency (s): time from the start of the interval to the relative threshold of SP detection.
#             Calculated using time_start_drop and drop_threshold (adjustable by the user in the "SP Latency Calculation..." section of the script, e.g.:
#             drop_threshold = 0.05 → SP detection requires at least a 5% drop from stimulus onset (default value = 0.5);
#             drop_threshold = 0.10 → requires a 10% drop, and so on)
#       - Delta SP (mV): difference between SP (minimum voltage) and latency point 
#       - Duration (s): time span during which the signal remains below 50% of the SP (width of the event measured between the points
#         where the signal drops below and rises back above 50% of the SP value)
#       - Area Under the Curve (AUC, in mV·s): numerical integration of the signal during the duration window, proportional to event energy. 
#       - Area_half: Restricted area only between t1/2 rise and t1/2 decay.
#       - t1/2 rise (s): absolute (s) and relative (s/mV) time to reach 50% of SP in the rising phase.
#       - t1/2 decay (s): absolute (s) and relative (s/mV) time from SP to when the signal returns to 50% in the falling phase
#       - Vm at t1 and t2: voltages at user-defined timepoints t1 and t2.
#       - % drop at t1 and t2: relative to SP
#       - Fits a straight line (y=a+bx) between:
#           └── Latency point and 50% rise.
#           └── Latency point and 80% rise
#           # Stores slope values (mV/s).
#           # Saves fitted curves as .txt files in a fitted_curves folder
#           # Plots fitted lines on the figure.
#	    - Derivative calculation (dV/dt):
#           # Computes the 1st derivative of the filtered signal within the selected interval (start, end)
#           # Identifies the minimum derivative value (most negative slope), corresponding to the maximum depolarization rate, 
#             and the time elapsed between the onset of the SP (latency) and the point of maximum depolarization rate 
#            (Δt = time_min_derivative − time_drop_start).
#           # Saves 1st derivative curves as .txt files in a derivative_curves folder
#           # Generates derivative plots highlighting the minimum slope point, annotated with its value in mV/s.
#           # Creation of output folders
#           	└──derivative_curves: stores the numerical derivative data as .txt files.
#           	└──derivative_figures: stores the derivative plots as .jpg images.

# 4. Plots for each sweep (folder: figuras_SP/)
#       Overlays raw and filtered signal.
#       Visually marks:
#           SP
#           Latency
#           50% SP threshold
#           Slope 50% rise
#           Slope 80 % rise
#           t1 and t2 points with labels
#           Draws fitted slope lines
#       Saves the figures in the figuras_SP folder.

# 5. Export of results: saves all sweep parameters in a list of dictionaries "sp_results"
#       As .csv (semicolon-separated), with selected columns rounded.
#       As .txt in a readable tabular format.

# 6. Summary visualizations: Extended histogram summary (folder: histogram_figures)
#       Histograms of SP values, Delta SP, latency, duration, relative t1/2 rise, % drop at t1 and t2.
#       Saved in histogram_figures/extended_summary_histograms.png.

#-----------------------------------------------------------------------------------------------------------------------------------------
    # Output structure:

            #├── SP_results.csv
            #├── SP_results.txt
            #├── SP_figures/
            #│   ├── file1.abf_sweep.png
            #│   ├── file2_sweep.png 
            #│   ├── file3_sweep.png ...
            #├── histogram_figures/
            #│   └── extended_summary_histograms.png
            #├── fitted_curves/
            #│   ├── file1.abf_sweep0_fit50
            #│   ├── file1.abf_sweep0_fit80
            #│   ├── file2.abf_sweep0_fit50...
            #├── derivative_curves/
            #├── derivative_figures/
    
#-----------------------------------------------------------------------------------------------------------------------------------------
# If you use this software, please cite:
# Corvalan, N. A. (2026). "EAG-SPy: An Automated Python Toolkit for Electroantennogram (EAG) and Sensillar Potential (SP) Analysis".
# Zenodo. DOI: 10.5281/zenodo.18760115
#-----------------------------------------------------------------------------------------------------------------------------------------

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyabf
from numpy.fft import fft, ifft, fftfreq

# ==============================================================================
# SOFTWARE METADATA
# ==============================================================================

SOFTWARE_NAME = "Automated Python Toolkit for Electroantennogram (EAG) and Sensillar Potential (SP) Analysis"
VERSION = "1.0.0"
AUTHOR = "Natalia Andrea Corvalán"
AFFILIATION = "Consejo Nacional de Investigaciones Científicas y Técnicas (CONICET) - Argentina"
ORCID = "https://orcid.org/0000-0003-2028-9609"
DOI = "10.5281/zenodo.18760115"
RELEASE_YEAR = "2026"
LICENSE = "GNU General Public License v3.0 or later"

CITATION_TEXT = (
    f"{AUTHOR} ({RELEASE_YEAR}). "
    f"{SOFTWARE_NAME} (Version {VERSION}). "
    f"Zenodo. https://doi.org/{DOI}"
)

def print_citation():
    print("\n" + "=" * 70)
    print(f"{SOFTWARE_NAME} - Version {VERSION}")
    print(f"Author: {AUTHOR}")
    print(f"Affiliation: {AFFILIATION}")
    print(f"ORCID: {ORCID}")
    print(f"License: {LICENSE}")
    print("\nIf you use this software in academic work, please cite:")
    print(CITATION_TEXT)
    print("=" * 70 + "\n")

def metadata_header():
    return f"# Generated with {SOFTWARE_NAME} v{VERSION} | DOI: https://doi.org/{DOI}\n"
print_citation()

# ------------------------------
# Gaussian low-pass filter
# ------------------------------
def gaussian_lowpass_filter(signal, fs, cutoff_hz):
    n = len(signal)
    freqs = fftfreq(n, d=1/fs)
    spectrum = fft(signal)
    gauss = np.exp(-0.5 * (freqs / cutoff_hz) ** 2)
    return np.real(ifft(spectrum * gauss))


# ------------------------------
# SP latency and kinetics 
# ------------------------------
def calculate_sp_latency(time, voltage, stimulus_start, drop_threshold=0.05):
    """
    Latency-based SP analysis.

    Definitions:
    - V_latency: voltage at first threshold crossing (kinetic onset)
    - V_min: absolute minimum voltage after latency
    - delta_sp: |V_min - V_latency|
    """

    idx_stimulus = np.argmin(np.abs(time - stimulus_start))

    time_post = time[idx_stimulus:]
    voltage_post = voltage[idx_stimulus:]

    V_start = voltage[idx_stimulus]

    # Preliminary minimum (robustness check)
    V_min_pre = np.min(voltage_post)
    preliminary_delta = abs(V_start - V_min_pre)
    if preliminary_delta <= 0:
        return (np.nan,) * 8

    # Threshold
    V_threshold = V_start - drop_threshold * preliminary_delta

    drop_indices = np.where(voltage_post <= V_threshold)[0]
    if len(drop_indices) == 0:
        return (np.nan,) * 8

    idx_drop_start = idx_stimulus + drop_indices[0]
    t_drop_start = time[idx_drop_start]
    V_latency = voltage[idx_drop_start]
    latency = t_drop_start - stimulus_start

    # Minimum AFTER latency
    voltage_after_latency = voltage[idx_drop_start:]
    if len(voltage_after_latency) == 0:
        return (np.nan,) * 8

    V_min = np.min(voltage_after_latency)
    idx_min = np.argmin(voltage_after_latency) + idx_drop_start

    delta_sp = abs(V_min - V_latency)
    if delta_sp == 0:
        return (np.nan,) * 8

    # ---- t1/2 rise ----
    V_half = V_latency + 0.5 * (V_min - V_latency)

    idx_half_rise = np.where(voltage[idx_drop_start:] <= V_half)[0]
    if len(idx_half_rise):
        t_half_rise = time[idx_drop_start + idx_half_rise[0]]
        t_half_rise_abs = t_half_rise - t_drop_start
        t_half_rise_rel = t_half_rise_abs / delta_sp
    else:
        t_half_rise_abs = np.nan
        t_half_rise_rel = np.nan

    # ---- t1/2 decay ----
    idx_half_decay = np.where(voltage[idx_min:] >= V_half)[0]
    if len(idx_half_decay):
        t_half_decay = time[idx_min + idx_half_decay[0]]
        t_half_decay_abs = t_half_decay - t_drop_start
        t_half_decay_rel = t_half_decay_abs / delta_sp
    else:
        t_half_decay_abs = np.nan
        t_half_decay_rel = np.nan

    return (
        latency,
        t_drop_start,
        delta_sp,   
        delta_sp,   
        t_half_rise_abs,
        t_half_rise_rel,
        t_half_decay_abs,
        t_half_decay_rel
    )


# ------------------------------
# User-defined parameters
# ------------------------------
if __name__ == "__main__":
    try:
        cutoff = float(input("Enter Gaussian low-pass filter cutoff frequency (-3 dB) in Hz: "))
        start = float(input("Enter start time of analysis interval (s): "))
        end = float(input("Enter end time of analysis interval (s): "))
        t1 = float(input("Enter first evaluation time t1 (s): "))
        t2 = float(input("Enter second evaluation time t2 (s): "))


        # ------------------------------
        # File preparation (standalone)
        # ------------------------------
        import os

        # Directory where the script or executable is located
        current_dir = os.path.dirname(os.path.abspath(__file__))

        # Force working directory
        os.chdir(current_dir)

        # Find .abf files in working directory (non-recursive)
        abf_files = [
            os.path.join(current_dir, f)
            for f in os.listdir(current_dir)
            if f.lower().endswith(".abf")
        ]

        if len(abf_files) == 0:
            print("\n⚠ No .abf files found in the working directory.")
        else:
            print(f"\nFound {len(abf_files)} .abf files")

        results = []

        # ------------------------------
        # Main processing
        # ------------------------------
        for file_path in abf_files:
            file_name = os.path.basename(file_path)
            abf = pyabf.ABF(file_path)
            fs = abf.dataRate

            print(f"\nAnalyzing file: {file_name}")

            for sweep in abf.sweepList:
                abf.setSweep(sweepNumber=sweep, channel=0)

                # --- RESET sweep-local variables (CRITICAL)
                fit_line_50 = fit_line_80 = None
                x_50 = y_50 = x_80 = y_80 = None

                time = abf.sweepX
                voltage = abf.sweepY

                filtered_voltage = gaussian_lowpass_filter(voltage, fs, cutoff)

                interval_idx = np.where((time >= start) & (time <= end))[0]
                if len(interval_idx) == 0:
                    continue

                time_i = time[interval_idx]
                volt_i = filtered_voltage[interval_idx]

                # Absolute SP (visual reference)
                V_min = np.min(volt_i)
                idx_sp = np.argmin(volt_i)
                time_sp = time_i[idx_sp]

                (
                    latency,
                    time_drop_start,
                    sp_amplitude,
                    delta_sp,
                    t1_2_rise_abs,
                    t1_2_rise_rel,
                    t1_2_decay_abs,
                    t1_2_decay_rel
                ) = calculate_sp_latency(time_i, volt_i, start)

                # Voltage at latency
                if not np.isnan(time_drop_start):
                    idx_lat = np.argmin(np.abs(time_i - time_drop_start))
                    V_latency = volt_i[idx_lat]
                else:
                    V_latency = np.nan

                # ---- Duration (50% ΔSP)
                if not np.isnan(delta_sp):
                    V_half = V_latency + 0.5 * (V_min - V_latency)
                    crossings = np.where(volt_i <= V_half)[0]
                    if len(crossings) >= 2:
                        duration = time_i[crossings[-1]] - time_i[crossings[0]]
                    else:
                        duration = np.nan
                else:
                    duration = np.nan

                # ---- Areas
                area = np.trapezoid(volt_i, time_i)

                if not np.isnan(t1_2_rise_abs) and not np.isnan(t1_2_decay_abs):
                    t_rise = time_drop_start + t1_2_rise_abs
                    t_decay = time_drop_start + t1_2_decay_abs
                    idx_r = np.argmin(np.abs(time_i - t_rise))
                    idx_d = np.argmin(np.abs(time_i - t_decay))
                    area_half = np.trapezoid(volt_i[idx_r:idx_d+1], time_i[idx_r:idx_d+1]) if idx_d > idx_r else np.nan
                else:
                    area_half = np.nan

                # ---- Voltages at t1 / t2
                idx_t1 = np.argmin(np.abs(time_i - t1))
                idx_t2 = np.argmin(np.abs(time_i - t2))
                volt_t1 = volt_i[idx_t1]
                volt_t2 = volt_i[idx_t2]

                drop_t1 = ((volt_t1 - V_latency) / delta_sp) * 100 if delta_sp != 0 else np.nan
                drop_t2 = ((volt_t2 - V_latency) / delta_sp) * 100 if delta_sp != 0 else np.nan

                results.append({
                    "file": file_name,
                    "sweep": sweep,
                    "SP_mV": V_min,
                    "delta_SP_mV": delta_sp,
                    "time_SP_s": time_sp,
                    "latency_s": latency,
                    "time_drop_start_s": time_drop_start,
                    "AUC_mV·s": area,
                    "Area t1/2 rise-decay_mV·s": area_half,
                    "t1/2_rise_abs_s": t1_2_rise_abs,
                    "t1/2_rise_rel_s/mV": t1_2_rise_rel,
                    "t1/2_decay_abs_s": t1_2_decay_abs,
                    "t1/2_decay_rel_s/mV": t1_2_decay_rel,
                    "duration_s": duration,
                    f"Vm_t1_{t1}_s": volt_t1,
                    f"Vm_t2_{t2}_s": volt_t2,
                    "drop_t1_%": drop_t1,
                    "drop_t2_%": drop_t2
                })

                print(f"- Sweep {sweep}: ΔSP={delta_sp:.2f} mV | Lat={latency:.3f}s | Dur={duration:.3f}s | AUC = {area:.3f} mV·s")

                # ------------------------------
                # Fitted curves (linear regression) between latency and 50% / 80% ΔSP
                # ------------------------------
                slope_50_fit = np.nan
                slope_80_fit = np.nan

                if not np.isnan(delta_sp) and not np.isnan(latency):

                    # Reference voltages (coherent with latency-based model)
                    V_50 = V_latency + 0.5 * (V_min - V_latency)
                    V_80 = V_latency + 0.8 * (V_min - V_latency)

                    # ---- 50% regression ----
                    idx_50 = np.where(volt_i <= V_50)[0]
                    if len(idx_50) > 0:
                        t_50 = time_i[idx_50[0]]

                        mask_50 = (time_i >= time_drop_start) & (time_i <= t_50)
                        x_50 = time_i[mask_50]
                        y_50 = volt_i[mask_50]

                        if len(x_50) > 1:
                            coeffs_50 = np.polyfit(x_50, y_50, 1)
                            slope_50_fit = coeffs_50[0]
                            fit_line_50 = np.polyval(coeffs_50, x_50)

                            # Save fitted curve
                            slopes_folder = os.path.join(current_dir, "fitted_curves")
                            os.makedirs(slopes_folder, exist_ok=True)
                            np.savetxt(
                                os.path.join(slopes_folder, f"{file_name}_sweep{sweep}_fit50.txt"),
                                np.column_stack((x_50, fit_line_50)),
                                fmt="%.6f",
                                header="Time(s)\tVoltage(mV)"
                            )

                    # ---- 80% regression ----
                    idx_80 = np.where(volt_i <= V_80)[0]
                    if len(idx_80) > 0:
                        t_80 = time_i[idx_80[0]]

                        mask_80 = (time_i >= time_drop_start) & (time_i <= t_80)
                        x_80 = time_i[mask_80]
                        y_80 = volt_i[mask_80]

                        if len(x_80) > 1:
                            coeffs_80 = np.polyfit(x_80, y_80, 1)
                            slope_80_fit = coeffs_80[0]
                            fit_line_80 = np.polyval(coeffs_80, x_80)

                            # Save fitted curve
                            slopes_folder = os.path.join(current_dir, "fitted_curves")
                            os.makedirs(slopes_folder, exist_ok=True)
                            np.savetxt(
                                os.path.join(slopes_folder, f"{file_name}_sweep{sweep}_fit80.txt"),
                                np.column_stack((x_80, fit_line_80)),
                                fmt="%.6f",
                                header="Time(s)\tVoltage(mV)"
                            )

                # Store slopes
                results[-1].update({
                    "slope_fit_50_mV/s": slope_50_fit,
                    "slope_fit_80_mV/s": slope_80_fit
                })

                # ------------------------------
                # Plot SP analysis
                # ------------------------------
                plt.figure()

                # Original vs filtered
                plt.plot(time_i, voltage[interval_idx], label="Original", alpha=0.4, color="gray")
                plt.plot(time_i, volt_i, label="Filtered", linewidth=1)

                # Absolute minimum (SP visual)
                plt.plot(time_sp, V_min, "ro", label=f"SP = {V_min:.2f} mV")

                # 50% ΔSP reference
                if not np.isnan(delta_sp):
                    V_half = V_latency + 0.5 * (V_min - V_latency)
                    plt.axhline(V_half, color="gray", linestyle="--", label="50% ΔSP")

                # Latency marker
                if not np.isnan(latency):
                    plt.axvline(time_drop_start, color="red", linestyle="--", alpha=0.6,
                                label=f"Latency = {latency:.3f}s")
                    plt.plot(time_drop_start, V_latency, "o", color="red")

                # Regression lines 
                if fit_line_50 is not None:
                    plt.plot(x_50, fit_line_50, "b--", linewidth=0.8, label="Fit 50% rise")

                if fit_line_80 is not None:
                    plt.plot(x_80, fit_line_80, "g--", linewidth=0.8, label="Fit 80% rise")

                # t1 / t2 markers
                plt.plot(time_i[idx_t1], volt_t1, "o", color="purple", label=f"t1 = {t1}s")
                plt.plot(time_i[idx_t2], volt_t2, "o", color="magenta", label=f"t2 = {t2}s")

                plt.xlabel("Time (s)")
                plt.ylabel("Voltage (mV)")
                plt.title(f"{file_name} - Sweep {sweep}")
                plt.legend()
                plt.grid(True)
                plt.tight_layout()

                figures_folder = os.path.join(current_dir, "sp_figures")
                os.makedirs(figures_folder, exist_ok=True)
                plt.savefig(os.path.join(figures_folder, f"{file_name}_sweep{sweep}.png"))
                plt.close()

        # ------------------------------
        # Derivative analysis: compute and save derivative curves
        # ------------------------------
        # Kinetic characterization: dV/dt reflects the rate of membrane potential change

        print("\nStarting derivative analysis...")

        # Output folders
        derivative_folder_txt = os.path.join(current_dir, "derivative_curves")
        derivative_folder_fig = os.path.join(current_dir, "derivative_figures")
        os.makedirs(derivative_folder_txt, exist_ok=True)
        os.makedirs(derivative_folder_fig, exist_ok=True)

        for file_path in abf_files:
            file_name = os.path.basename(file_path)
            abf = pyabf.ABF(file_path)
            fs = abf.dataRate

            for sweep in abf.sweepList:
                abf.setSweep(sweepNumber=sweep, channel=0)

                time = abf.sweepX
                voltage = abf.sweepY

                # Apply same Gaussian filter used for SP analysis
                filtered_voltage = gaussian_lowpass_filter(voltage, fs, cutoff)

                # Restrict to analysis interval
                interval_indices = np.where((time >= start) & (time <= end))[0]
                if len(interval_indices) == 0:
                    continue

                time_interval = time[interval_indices]
                voltage_interval = filtered_voltage[interval_indices]

                # First derivative dV/dt
                derivative = np.gradient(voltage_interval, time_interval)

                # Most negative derivative = maximal depolarization rate
                min_derivative = np.min(derivative)
                idx_min_derivative = np.argmin(derivative)
                time_min_derivative = time_interval[idx_min_derivative]

                # Update corresponding SP result entry
                for r in results:
                    if r["file"] == file_name and r["sweep"] == sweep:
                        r["min_derivative_mV/s"] = min_derivative
                        r["time_min_derivative_s"] = time_min_derivative

                        if not np.isnan(r.get("time_drop_start_s", np.nan)):
                            r["delta_t_latency_min_derivative_s"] = (
                                time_min_derivative - r["time_drop_start_s"]
                            )
                        else:
                            r["delta_t_latency_min_derivative_s"] = np.nan
                        break

                # Save derivative curve as TXT
                txt_out = os.path.join(
                    derivative_folder_txt,
                    f"{file_name}_sweep{sweep}_derivative.txt"
                )
                np.savetxt(
                    txt_out,
                    np.column_stack((time_interval, derivative)),
                    fmt="%.6f",
                    header="Time(s)\tdV/dt(mV/s)"
                )

                # Plot derivative curve
                plt.figure(figsize=(8, 4))
                plt.plot(time_interval, derivative, linewidth=1)
                plt.plot(
                    time_min_derivative,
                    min_derivative,
                    marker="x",
                    color="black",
                    markersize=8,
                    linestyle="None",
                    label=f"Min dV/dt = {min_derivative:.2f} mV/s"
                )

                plt.title(f"Derivative – {file_name} – Sweep {sweep}")
                plt.xlabel("Time (s)")
                plt.ylabel("dV/dt (mV/s)")
                plt.legend()
                plt.grid(True)
                plt.tight_layout()

                fig_out = os.path.join(
                    derivative_folder_fig,
                    f"{file_name}_sweep{sweep}_derivative.jpg"
                )
                plt.savefig(fig_out, dpi=300)
                plt.close()

        print(
            "Derivative analysis completed.\n"
            f"- TXT curves saved in: {derivative_folder_txt}\n"
            f"- Figures saved in: {derivative_folder_fig}"
        )

        # ------------------------------
        # Export results
        # ------------------------------
        df = pd.DataFrame(results)

        # Ensure all expected columns exist (avoid KeyError)
        expected_columns = [
            # Core SP metrics
            "SP_mV",
            "delta_SP_mV",
            "time_SP_s",
            "latency_s",
            "duration_s",

            # Areas
            "AUC_mV·s",
            "Area t1/2 rise-decay_mV·s",

            # Kinetics
            "t1/2_rise_abs_s",
            "t1/2_rise_rel_s/mV",
            "t1/2_decay_abs_s",
            "t1/2_decay_rel_s/mV",

            # Regressions
            "slope_fit_50_mV/s",
            "slope_fit_80_mV/s",

            # Voltage sampling
            f"Vm_t1_{t1}_s",
            f"Vm_t2_{t2}_s",

            # Relative drops
            "drop_t1_%",
            "drop_t2_%",

            # Derivative
            "min_derivative_mV/s",
            "time_min_derivative_s",
            "delta_t_latency_min_derivative_s",
        ]

        for col in expected_columns:
            if col not in df.columns:
                df[col] = np.nan

        # Metadata
        df.attrs["software"] = SOFTWARE_NAME
        df.attrs["version"] = VERSION
        df.attrs["doi"] = f"https://doi.org/{DOI}"

        # Rounding
        round_columns = {
            "SP_mV": 3,
            "delta_SP_mV": 3,
            "time_SP_s": 3,
            "latency_s": 3,
            "duration_s": 3,
            "AUC_mV·s": 3,
            "Area t1/2 rise-decay_mV·s": 3,
            "t1/2_rise_abs_s": 3,
            "t1/2_rise_rel_s/mV": 3,
            "t1/2_decay_abs_s": 3,
            "t1/2_decay_rel_s/mV": 3,
            f"Vm_t1_{t1}_s": 2,
            f"Vm_t2_{t2}_s": 2,
            "drop_t1_%": 2,
            "drop_t2_%": 2,
            "slope_fit_50_mV/s": 3,
            "slope_fit_80_mV/s": 3,
            "min_derivative_mV/s": 2,
            "time_min_derivative_s": 4,
            "delta_t_latency_min_derivative_s": 4
        }

        for col, dec in round_columns.items():
            if col in df.columns:
                df[col] = df[col].round(dec)

        # CSV export
        csv_path = os.path.join(current_dir, "sp_results.csv")
        with open(csv_path, "w", encoding="utf-8") as f:
            f.write(metadata_header())
            df.to_csv(f, sep=";", index=False)

        print(f"\nResults saved to: {csv_path}")

        # TXT export
        txt_path = os.path.join(current_dir, "sp_results.txt")
        with open(txt_path, "w", encoding="utf-8") as f:
            f.write(metadata_header())
            df.to_csv(f, sep="\t", index=False)

        print(f"Results also saved as .txt at: {txt_path}")

        # ------------------------------
        # Summary histograms
        # ------------------------------
        hist_folder = os.path.join(current_dir, "histogram_figures")
        os.makedirs(hist_folder, exist_ok=True)

        hist_path = os.path.join(hist_folder, "extended_summary_histograms.png")

        plt.figure(figsize=(18, 10))

        vars_plot = [
            ("SP_mV", "SP (mV)", "skyblue"),
            ("delta_SP_mV", "ΔSP (mV)", "lightcoral"),
            ("latency_s", "Latency (s)", "orange"),
            ("duration_s", "Duration (s)", "green"),
            ("t1/2_rise_rel_s/mV", "t1/2 rise (s/mV)", "purple"),
            ("drop_t1_%", f"Drop (% ΔSP) at t1 = {t1}s", "teal"),
            ("drop_t2_%", f"Drop (% ΔSP) at t2 = {t2}s", "brown"),
            ("min_derivative_mV/s", "Min dV/dt (mV/s)", "black"),
        ]

        for i, (col, title, color) in enumerate(vars_plot, start=1):
            plt.subplot(2, 4, i)
            plt.hist(df[col].dropna(), bins=15, edgecolor="black")
            plt.title(title)
            plt.xlabel(title)
            plt.ylabel("Frequency")

        plt.tight_layout()

        plt.figtext(
            0.99, 0.01,
            f"{SOFTWARE_NAME} v{VERSION} | DOI: https://doi.org/{DOI}",
            ha="right",
            fontsize=6
        )

        plt.savefig(hist_path, dpi=300)
        plt.close()

        print(f"Histograms saved in: {hist_folder}")

        print("\nAnalysis completed.")

        pass
    
    except Exception as e:
        print("\nUnexpected error occurred:")
        print(str(e))

        import traceback
        traceback.print_exc()

    finally:
        input("\nPress Enter to exit...")