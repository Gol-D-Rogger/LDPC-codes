#!/usr/bin/python3
import glob
import os
import re

import matplotlib.pyplot as plt
from binomial_distrib import *

folder_prefix = "4112_752"
str_filt_keep = "user_data"
str_filt_redact = "XXXXXXXXXXXXXX"
mode = "fail_rate"
add_binomial_fail_prob_curves = False
y_axis_text = "UBER" if mode == "uber" else "Codeword Failure Rate"
binom_curve_count = 0


def extract_error_rates(directory):
    results = {}
    for filename in sorted(os.listdir(directory)):
        if filename.startswith("rber_") and filename.endswith(".txt"):
            rber_match = re.search(r"rber_(\d+p\d+)\.txt", filename)
            if rber_match:
                rber = rber_match.group(1).replace("p", ".")
                filepath = os.path.join(directory, filename)
                try:
                    with open(filepath, "r") as f:
                        for line in f:
                            if "Overall Decoding Failure Rate:" in line:
                                failure_rate_match = re.search(r"Overall Decoding Failure Rate: (\d+\.\d+)", line)
                                if failure_rate_match:
                                    failure_rate = float(failure_rate_match.group(1))
                                    num_bytes_cw = sum([int(elem) for elem in folder_prefix.split("_")[:2]])
                                    if mode == "uber":
                                        results[float(rber)] = failure_rate / (num_bytes_cw * 8)
                                    else:
                                        results[float(rber)] = failure_rate
                except FileNotFoundError:
                    print(f"File not found: {filepath}")
                except Exception as e:
                    print(f"An error occurred while processing {filepath}: {e}")
    return results


def plot_multiple_error_rates(all_rates, labels, colors=None):
    plt.figure(figsize=(12, 8))
    ax = plt.gca()
    global binom_curve_count

    for i, rber_rates in enumerate(all_rates):
        if rber_rates:
            rbers = sorted(rber_rates.keys())
            rates = [rber_rates[r] for r in rbers]
            linestyle = "--" if labels[i].startswith("Failure probability") else "-"
            binom_curve_count += 1 if labels[i].startswith("Failure probability") else 0

            if binom_curve_count == 1:
                plt.gca().set_prop_cycle(None)

            plot_args = {"marker": "o", "linestyle": linestyle, "label": labels[i], "linewidth": 2}
            if colors and i < len(colors) and colors[i]:
                plot_args["color"] = colors[i]
            plt.plot(rbers, rates, **plot_args)

    plt.yscale("log")
    plt.xlabel("RBER", fontsize=14)
    plt.ylabel(f"{y_axis_text}", fontsize=14)
    plt.title("ScaleFlux, Inc. Proprietary and Confidential\nLDPC Decoder", fontsize=15)
    ax.tick_params(labelsize=11)
    plt.grid(True, which="major", linestyle="-", color="blue")
    plt.grid(True, which="minor", linestyle=":", axis="x")
    ax.yaxis.set_major_locator(plt.LogLocator(base=10.0, numticks=15))
    ax.xaxis.set_major_locator(plt.MultipleLocator(0.005))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.001))
    plt.legend()
    plt.savefig("customer_fig.png", dpi=300)
    plt.show()


def get_failure_prob_binomial(bit_error_probability, total_bits, max_error_correction_capability):
    sum_fail_prob = 0
    for k_i in range(max_error_correction_capability + 1, total_bits + 1):
        res = binom_pmf(p=bit_error_probability, k=k_i, n=total_bits)
        sum_fail_prob += res
    return sum_fail_prob


def get_binomial_fail_prob_curve(bit_error_range, total_bits, max_error_correction_capability):
    curve = {}
    for bit_error_prob in bit_error_range:
        res = get_failure_prob_binomial(bit_error_prob, total_bits, max_error_correction_capability)
        curve[bit_error_prob] = res
    return curve


if __name__ == "__main__":
    base_dir = "./_curve_data/"
    all_error_rates = []
    legend_labels = []

    for folder_prefix_str in ["4112_496", "4112_752"]:
        folder_prefix = folder_prefix_str
        for folder in sorted(glob.glob(os.path.join(base_dir, f"{folder_prefix}*"))):
            if os.path.isdir(folder):
                if str_filt_keep in folder and str_filt_redact not in folder:
                    print(f"Processing folder: {folder}")
                    error_rates = extract_error_rates(folder)
                    if error_rates:
                        all_error_rates.append(error_rates)
                        legend_label = os.path.basename(folder).replace(f"{folder_prefix}_", "")
                        legend_labels.append(legend_label)

    if add_binomial_fail_prob_curves:
        num_bytes_cw = sum([int(elem) for elem in folder_prefix.split("_")[:2]])
        for t, start, count in [(465, 0.009, 18), (860, 0.018, 20), (900, 0.0195, 20), (920, 0.02, 20)]:
            curve = get_binomial_fail_prob_curve([start + 0.0005 * x for x in range(count)], num_bytes_cw * 8, t)
            all_error_rates.append(curve)
            legend_labels.append(f"Failure probability for t = {t}")

    colors = ["#1f77b4", "gold", "green"]
    if all_error_rates:
        plot_multiple_error_rates(all_error_rates, legend_labels, colors)
    else:
        print(f"No data found in folders starting with '{folder_prefix}'")
