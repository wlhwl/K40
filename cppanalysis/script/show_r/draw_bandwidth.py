#! /usr/bin/env python3
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

if __name__ == "__main__":
    wf_hz2MB = 31 * 500 * 16 / 1024 / 1024 / 8
    tdc_hz2MB = 24 * 100 / 8 / 1024 / 1024
    one_time_sim = 1.56524 * 10 ** 7  # ns
    coin_num = [1, 2, 3, 4, 5, 6]
    time_window = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
    csvfile = "trigger_num.csv"
    trigger_list = pd.read_csv(csvfile, header=None)
    run_num = len(trigger_list[0])
    print("run %d times" % run_num)
    coin2 = []
    for i in range(10, 20):
        coin2.append(sum(trigger_list[i][:]))

    window10 = []
    for j in range(0, 6):
        window10.append(sum(trigger_list[10*j + 4][:]))

    coin2_trigger_rate = np.divide(coin2, one_time_sim * run_num) * 10 ** 9
    window10_trigger_rate = np.divide(window10, one_time_sim * run_num) * 10 ** 9

    coin22MB = coin2_trigger_rate * (wf_hz2MB + tdc_hz2MB)
    window102MB = window10_trigger_rate * (wf_hz2MB + tdc_hz2MB)

    plt.figure(dpi=800)
    plt.plot(time_window, coin2_trigger_rate, "--", color='blue')
    plt.xlabel("time window [ns]")
    plt.ylabel("MB/s per hDOM")
    plt.savefig("coin2.png")

    plt.figure(dpi=800)
    plt.plot(coin_num, window10_trigger_rate, "--", color='blue')
    plt.xlabel("coincidence photon number")
    plt.ylabel("MB/s per hDOM")
    plt.savefig("window10.png")

