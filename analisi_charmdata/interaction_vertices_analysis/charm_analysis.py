import numpy as np
import math
import functools as ft
from datetime import date
from itertools import repeat
from ROOT import TFile, TCanvas, TH1F, TF1, TPaveText, TAxis


def average_std_dev(x, x_mean):
    sum = 0
    for i, x_i in enumerate(x):
        sum += math.fabs(x_i - x_mean)
    return (1 / (i + 1)) * sum


def bin_content(hist):
    x = []
    for num_bins in range(hist.GetNbinsX()):
        x.append(hist.GetBinContent(num_bins + 1))
    return x


def bin_avg_dev_std_error(hist, avg_dev_std):
    for num_bins in range(hist.GetNbinsX()):
        error = avg_dev_std
        hist.SetBinError(num_bins + 1, error)
    return hist


def bin_res_error(hist, data_diff):
    for num_bins, residuals in enumerate(data_diff):
        error = residuals
        hist.SetBinError(num_bins + 1, error)
    return hist


def bin_poisson_error(hist):
    for num_bins in range(hist.GetNbinsX()):
        error = math.sqrt(hist.GetBinContent(num_bins + 1))
        hist.SetBinError(num_bins + 1, error)
    return hist


def bin_poisson_error_MC_scale_4(hist, err_MC, MC_tot_quarter_over_nquarter, MC_norm_const):
    #hist.Scale(1/MC_tot_quarter_over_nquarter)
    for num_bins in range(hist.GetNbinsX()):
        error = MC_norm_const * MC_tot_quarter_over_nquarter * math.sqrt(hist.GetBinContent(num_bins + 1))
        err_MC.append(error)
        hist.SetBinError(num_bins + 1, error)
        hist.SetBinContent(num_bins + 1, MC_tot_quarter_over_nquarter * hist.GetBinContent(num_bins + 1))
    return hist

def make_hist_conf(iconf, irun, file_dir):
    file_name = "CH" + str(iconf + 1) + "_vz_mc_data_R" + str(irun) + ".root"
    file_in = TFile(file_dir + file_name)
    canvas = file_in.Get("c3")
    hDT = canvas.GetPrimitive("DT" + str(irun) + "_vz")
    if iconf < 2:
        hDT.Rebin(3)
    else:
        hDT.Rebin(6)
    return hDT.GetNbinsX(), hDT.GetXaxis().GetXmin(), hDT.GetXaxis().GetXmax()


today = date.today()
print("Today's date:", today)

# work directory
file_dir = "/home/vale/Desktop/SHiP/analysis/20200610_proton_study/cut_0.05_with_tracks/vz_not_scaled/"

out_file = open("log_charm_errors.txt", "w")
out_file.write(f"Log file - Stima degli errori (V.G. {today})")

# RUN available for each configuration
CH_file = [[1, 2, 4, 5], [2, 3, 4], [1], [1], [1], [1]]
CH1_file = [1, 2, 4, 5]
CH2_file = [2, 3, 4]
CH3_file = [1]
CH4_file = [1]
CH5_file = [1]
CH6_file = [1]
conf_file = []
conf_file.append(CH1_file)

conf_file.append(CH2_file)
conf_file.append(CH3_file)
conf_file.append(CH4_file)
conf_file.append(CH5_file)
conf_file.append(CH6_file)

nconf = 6
run_scale_factor = [[] for x in range(nconf)]
# MC POT PER CONF
nPOT_MC = [135000, 135000, 35000, 27500, 40000, 22500]
# CONSTANT
MC_norm_const = [(nPOT_MC[0] / iPOT)  for iPOT, iFile in zip(nPOT_MC, CH_file)]
print(MC_norm_const)
norm_const = [(nPOT_MC[0] / iPOT) / len(iFile) for iPOT, iFile in zip(nPOT_MC, CH_file)]
print(norm_const)
MC_tot_quarter_over_nquarter = [4, 4, 4, 4, 4, 4]
# DT POT PER CONF
DT_POT = [[136270, 110352, 76307, 73575], [35763, 38853, 43131], [29586], [24578], [25593], [20443]]
for ich, ipot_mc in enumerate(nPOT_MC):
    for ipot_run in DT_POT[ich]:
        run_scale_factor[ich].append(ipot_mc/ipot_run)
print(run_scale_factor)

out_file.write(f"\n\nGeneral info")
out_file.write(f"\nNormalization constant for each configuration '(mcpot1/mcpotx)/nScan': {norm_const}")
out_file.write(f"\nScale factor for each run '(mcpot_conf/data_pot_run)': {run_scale_factor}")
#

# conf_avg_std_dev = []
all_data_point = []
err_eff = []
err_fit = []
err_tot = []
err_MC = []
err_MC_pr = []
err_MC_hd = []
tot_points = 0

hCH_DT = []
hCH_MC = []

MC_file_name = "full_vz_1quarter.root"
MC_file_in = TFile(file_dir + MC_file_name)

for index, conf_number in enumerate(conf_file):

    hMC_pr = MC_file_in.Get("vz_ch" + str(index+1) + "_pr")
    hMC_hd = MC_file_in.Get("vz_ch" + str(index + 1) + "_hd")

    if index < 2 :
        hMC_pr.Rebin(3)
        hMC_hd.Rebin(3)
    else:
        hMC_pr.Rebin(6)
        hMC_hd.Rebin(6)

    bin_poisson_error_MC_scale_4(hMC_pr, err_MC_pr, MC_tot_quarter_over_nquarter[index], MC_norm_const[index])
    bin_poisson_error_MC_scale_4(hMC_hd, err_MC_hd, MC_tot_quarter_over_nquarter[index], MC_norm_const[index])
    print("ERRORI_MC: ",err_MC_pr, err_MC_hd)

    list_avg_dev_std = []
    conf_counts = []
    eff_corr = []
    out_file.write("\n\nAnalyzing configuration CH" + str(index + 1))

    #nbins, xmin, xmax = make_hist_conf(index, conf_number[0], file_dir)
    #hCONF_DT = TH1F("conf_dt" + str(index + 1), "conf_dt" + str(index + 1), nbins, xmin, xmax)
    #hCONF_MC = TH1F("conf_mc" + str(index + 1), "conf_mc" + str(index + 1), nbins, xmin, xmax)

    #hCH_MC.append(TH1F("MC_CH" + str(index + 1), "MC_CH" + str(index + 1), nbins, xmin, xmax))
    #hCH_DT.append(TH1F("DT_CH" + str(index + 1), "DT_CH" + str(index + 1), nbins, xmin, xmax))
    #print("thifffffffff ", hCH_MC)

    for irun, run_number in enumerate(conf_number):
        file_name = "CH" + str(index + 1) + "_vz_mc_data_R" + str(run_number) + ".root"
        print(f"Reading file {file_name}")
        out_file.write(f"\nReading file {file_name}")

        file_in = TFile(file_dir + file_name)
        canvas = file_in.Get("c3")
        hMC = canvas.GetPrimitive("MC" + str(index + 1) + "_vz")
        hDT = canvas.GetPrimitive("DT" + str(run_number) + "_vz")
        hDT_corr = hDT.Clone()
        if index < 2 :
            hDT.Rebin(3)
            hMC.Rebin(3)
        else:
            hDT.Rebin(6)
            hMC.Rebin(6)
        hDT.Scale(run_scale_factor[index][irun])
        hMC.Fit("pol0", "Q")
        hDT.Fit("pol0", "Q")
        print("run_scale_factor ", run_scale_factor[index][irun])
        const_MC = hMC.GetFunction("pol0").GetParameter(0)
        const_DT = hDT.GetFunction("pol0").GetParameter(0)  # err_eff mi serve riscalare ai pot
        hMC.GetFunction("pol0").SetLineColor(1)
        hDT.GetFunction("pol0").SetLineColor(3)
        eff_corr.append(const_MC / const_DT)
        hDT.Scale(1/run_scale_factor[index][irun])
        hDT.Fit("pol0", "Q")
        const_DT_fit = hDT.GetFunction("pol0").GetParameter(0)  # per err_fit
        print("constanti ",const_DT, const_DT_fit)
        # print(f"MC constant is {const_MC}, DATA constant is {const_DT}")
        # print(f"Avg std dev correction is {eff_corr}")
        out_file.write("\nPol0 fit results:")
        out_file.write(f"\nMC constant is {const_MC}, DATA constant is {const_DT}")
        data_point = bin_content(hDT)  # counts of the run
        conf_counts.append(data_point)  # counts for all runs
        #data_diff = [abs(x - const_DT) for x in data_point]
        # print(data_point, data_diff)
        # print(f"Bin residuals: {data_diff}")
        average_dev_std = average_std_dev(data_point, const_DT_fit)
        list_avg_dev_std.append(average_dev_std)
        # print(f"The average std dev is {average_dev_std}")
        out_file.write(f"\nThe average std dev is {average_dev_std}")
        # print(list_avg_dev_std)

    #print("bincontent ", index, hCH_MC.GetBinContent(1))
    # Lista con tutti i bin content per ciasun run
    all_data_point.append(data_point)
    out_file.write(f"\n\nFinal results for CH" + str(index + 1))
    out_file.write(f"\nIl numero di run acquisiti è: {irun + 1}")
    conf_bins = len(all_data_point[index])
    out_file.write(f"\nIl numero di bin della configurazione è {conf_bins}")
    out_file.write(f"\nFit pol0 efficiencies correction for each run are {eff_corr[tot_points:]}")

    # CALCOLO ERRORE DI EFFICIENZA
    # Trasformo la lista in un array numpy
    arr = np.array(conf_counts)
    # Mi calcolo la trasposta, ogni sottolista ha i get bin content dell'i-esimo bin di ciascun run
    tr_arr = np.transpose(arr)
    print("transposta :",tr_arr)

    # Bincontent sommato su tutti i run di una configurazione (N_iTot)
    tot_counts = [sum(icounts) for icounts in zip(*conf_counts)]
    print(f"tot bin counts {tot_counts}")

    double_product = [a*b for a, b in zip(eff_corr, run_scale_factor[index])]
    print("double_product: ", double_product, eff_corr, run_scale_factor[index])

    for itot_counts, icounts in zip(tot_counts, tr_arr):
        print("iscale_check: ", itot_counts, norm_const[index], eff_corr, icounts)
        err_eff.append(math.sqrt(itot_counts) * (norm_const[index] / itot_counts) * (np.dot(icounts, double_product)))
        #np.linalg.multi_dot([icounts, eff_corr, run_scale_factor[index]])
        #scl_prod = np.dot(icounts, run_scale_factor[index])
        print("scalar_product: ", np.dot(icounts, double_product))
    # print(f"Bin errors from efficiency {err_eff}")

    conf_scale_factor = sum(run_scale_factor[index])
    weights_scale_factor = [iscale/conf_scale_factor for iscale in(run_scale_factor[index])]
    print(weights_scale_factor)
    out_file.write(f"\nBin errors from efficiency {err_eff[tot_points:]}")
    out_file.write(f"\nThe average std deviation for each run is {list_avg_dev_std}")
    conf_avg_std_dev = (norm_const[index]) * math.sqrt(sum(iw * iavg * iavg for iw, iavg in zip(weights_scale_factor,list_avg_dev_std)))
    out_file.write(f"\nThe mean error from the average std dev is {conf_avg_std_dev}")
    err_fit.extend(repeat(conf_avg_std_dev, conf_bins))

    # ERRORE SUL MC
    bin_poisson_error_MC_scale_4(hMC, err_MC, MC_tot_quarter_over_nquarter[index], MC_norm_const[index])
    print(err_MC)
    out_file.write(f"\nMC bin errors from Poisson {err_MC[tot_points:]}")
    # hMC.Scale(4)  # riporto al confronto coi dati

    tot_points = (ft.reduce(lambda count, l: count + len(l), all_data_point, 0))
    print(tot_points)

out_file.write(f"\n\nFinal results")
out_file.write(f"\nIl numero di configurazioni acquisite è: {index + 1}")
# out_file.write(f"\n\nErrors list: 'err_eff - err_fit - err_tot'")
out_file.write(f"\n\nErrors list: 'ibin - err_DATA - err_MC - err_MC_pr - err_MC_hd'")

err_tot = ([math.sqrt((a * a + b * b)) for a, b in zip(err_eff, err_fit)])
#print(err_tot)

for points, (i, j, k, l, m, n) in enumerate(zip(err_eff, err_fit, err_tot, err_MC, err_MC_pr, err_MC_hd)):
    print(points, i, j, k, l, m, n)
    # out_file.write(f"\n{points} {i} {j} {k} {l}")
    out_file.write(f"\n{points} {k} {l} {m} {n}")

#print(hCH_MC)
#input()
"""
for i, ich in enumerate(hCH_MC):
    c1.cd(i + 1)
    hCH_MC[0].Draw("")

input()

for ibins in range(0,hMC.GetNbinsX()):
    print(ibins, hMC.GetBinError(ibins+1))


hMC.SetTitle("")
#hDT.Draw("sames")
t1 = TPaveText(0.32,0.92,0.63,0.98,"brNDC")
t1.AddText("Original comparison")
t1.Draw("same")
c1.cd(2)
hMC.Draw("hist")
hDT_corr2.Draw("sames")
t2 = TPaveText(0.32,0.92,0.63,0.98,"brNDC")
t2.AddText("Correction from integrals ratio")
t2.Draw("same")
c1.cd(3)
hMC.Draw("hist")
hDT_corr3.Draw("sames")
t3 = TPaveText(0.32,0.92,0.63,0.98,"brNDC")
t3.AddText("Correction from residuals")
t3.Draw("same")
c1.cd(4)
hMC.Draw("hist")
hDT_corr.Draw("sames")
t4 = TPaveText(0.32,0.92,0.63,0.98,"brNDC")
t4.AddText("Correction from average standard deviation")
t4.Draw("same")

print("Premere qualsiasi tasto per uscire")
input()
"""
