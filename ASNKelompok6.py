import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import plotly.graph_objs as go
from streamlit_option_menu import option_menu
import math
import streamlit as st 




########
column_names = ['ECG']
data=pd.read_csv('ECG5minutes.txt',delimiter="\t", names=column_names)
data["sample interval"] = np.arange(len(data))
data["elapsed time"] = (data["sample interval"])*(1/200)
x=data["elapsed time"]
y=data["ECG" ] - (sum(data["ECG" ]/len(data["ECG"]))) #agar turun ke baseline

fs=int(round(1/(data.iloc[1,2]-data.iloc[0,2])))
jumlahdata = int(np.size(x))

#LPF
fc_lpf = 11
fc_lpf=float(fc_lpf)

lpf_ecg = np.zeros(jumlahdata) 
for n in range(3):
    lpf_ecg[-n] = lpf_ecg[0]
    y[-n] = y[0]

# Coefficients of LPF
T = 1 / fs
w = 2 * math.pi * fc_lpf
a0 = w**2
a1 = 2*(w**2)
b1 = (8/(T**2)) - (2*(w**2))
C0 = ((4/(T**2)) - ((2 * (math.sqrt(2))*w)/T) + (w**2))
C1 = ((4/(T**2)) + ((2 * (math.sqrt(2))*w)/T) + (w**2))

# BUTTERWORTH LOWPASS FILTER EQUATION
for n in range(jumlahdata):
    lpf_ecg[n] = ((b1 * lpf_ecg[n-1]) - (C0 * lpf_ecg[n-2]) + (a0 * y[n]) + (a1 * y[n-1]) + (a0 * y[n-2])) / C1
#fc_hpf =input("FREQUENCY CUT-OFF FOR HIGHPASS FILTER :")
fc_hpf = 5
fc_hpf=float(fc_hpf)

#HPF
hpf_ecg = np.zeros(np.size(lpf_ecg))
for n in range(3):
    hpf_ecg[-n] = hpf_ecg[0]

# Coefficients of LPF
T = 1 / fs
w = 2 * math.pi * fc_hpf
e0 = 4*T
e1 = 8*T
e2 = 4*T
d0 = ((2*(w**2)*(T**2))-8)
d1 = (((w**2)*(T**2)) - (2 * (math.sqrt(2)) * T * w)+4)
d2 = (((w**2)*(T**2)) + (2 * (math.sqrt(2)) * T* w)+4)
# BUTTERWORTH LOWPASS FILTER EQUATION
for n in range(np.size(lpf_ecg)):
    hpf_ecg[n] = ((e0 * lpf_ecg[n]) - (e1 * lpf_ecg[n-1]) +(e2 * lpf_ecg[n-2])-(d0 * hpf_ecg[n-1])- (d1 * hpf_ecg[n-2]))/ d2

#DERIVATIVE
drv=np.zeros(np.size(hpf_ecg))
for n in range (np.size(hpf_ecg)-2):
   drv[n]= (1/8)*(-(hpf_ecg[n-2]) - (2*hpf_ecg[n-1]) + (2*hpf_ecg[n+1]) + (hpf_ecg[n+2]))

# SQUARING PROCEDURE METHOD
sqr=np. zeros(np.size(drv) )
for n in range (np.size(drv)):
  sqr[n]=(drv[n])**2
# MAV PROCEDURE METHOD
w = 10 
mav = np.zeros(np.size(sqr))
for n in range(np.size(sqr)):
    for i in range(w):
        mav[n] = mav[n] + sqr[n - i]
    mav[n] = mav[n] / w
tinggi=0
tinggi=np.zeros(np.size(mav))
for n in range (np. size(mav) -1): 
    if (tinggi < mav[n]) .all():
      tinggi [n]=mav[n]

thr=tinggi*0.5
thrqrs=np.zeros(np.size(mav))
for n in range (np. size(mav)-1):
  if (mav[n] >= thr).all():
     thrqrs [n]=1
  elif (mav[n]<thr).all():
    thrqrs [n]=0
      
#NUMBERS OF R TO R CALCULATIONS
ptp = 0
waktu = np.zeros(np.size(thrqrs))
selisih = np.zeros(np.size(thrqrs))

for n in range(np.size(thrqrs) - 1):
    if thrqrs[n] < thrqrs[n + 1]:
        waktu[ptp] = n / fs;
        selisih[ptp] = waktu[ptp] - waktu[ptp - 1]
        ptp += 1
ptp = ptp - 1

#CALCULATION OF THE AMOUNT OF R
j = 0
peak = np.zeros(np.size(thrqrs))
for n in range(np.size(thrqrs)-1):
    if thrqrs[n] == 1 and thrqrs[n-1] == 0:
        peak[j] = n
        j += 1



#BPM CALCULATIONS':
temp = 0
interval = np.zeros(np.size(thrqrs))
BPM = np.zeros(np.size(thrqrs))

for n in range(ptp):
    interval[n] = (peak[n] - peak[n-1]) * (1/fs)
    BPM[n] = 60 / interval[n]
    temp = temp+BPM[n]
    rata = temp / (n - 1)
    
#TIME DOMAIN
RR_SDNN=0
for n in range (ptp):
   RR_SDNN += (((selisih[n])-(60/rata))**2)

SDNN = math.sqrt (RR_SDNN/ (ptp-1))

RR_RMSSD=0
for n in range (ptp):
   RR_RMSSD += ((selisih[n+1]-selisih[n])**2)
RMSSD =  math. sqrt (RR_RMSSD/(ptp-1))

# FIND NN50 ALGORITHM
NN50 = 0

for n in range (ptp): 
    if (abs(selisih[n+1]-selisih[n])>0.05):
      NN50 +=1
pNN50 = (NN50/ (ptp-1)) *100 

dif = 0
for n in range (ptp):
  dif += abs(selisih[n]-selisih[n+1])
RRdif = dif/(ptp-1)

RR_SDSD = 0
for n in range (ptp):
  RR_SDSD += (((abs(selisih[n]-selisih[n+1]))-RRdif)**2)
SDSD = math.sqrt(RR_SDSD/(ptp-2))

bpm_rr = np.zeros(ptp)
for n in range (ptp):
  bpm_rr[n] = 60/selisih[n]
  if bpm_rr [n]>100:
    bpm_rr[n]=rata

n = np. arange(0,ptp,1,dtype=int)

bpm_rr_baseline = bpm_rr - 70

# Fungsi Fourier Transform dan perhitungan frekuensi
def fourier_transform(signal):
    N = len(signal)
    fft_result = np.zeros(N, dtype=complex)
    for k in range(N):
        for n in range(N):
            fft_result[k] += signal[n] * np.exp(-2j * np.pi * k * n / N)
    return fft_result

def calculate_frequency(N, sampling_rate):
    return np.arange(N) * sampling_rate / N

sampling_rate = 1

# Subset 0-50
n_subset = n[0:50]
bpm_rr_baseline_subset = bpm_rr_baseline[0:50]
M = len(bpm_rr_baseline_subset) - 1
hamming_window = np.zeros(M+1)
for i in range(M+1):
    hamming_window[i] = 0.54 - 0.46 * np.cos(2 * np.pi * i / M)
bpm_rr_baseline_windowed = bpm_rr_baseline_subset * hamming_window
fft_result = fourier_transform(bpm_rr_baseline_windowed) # Compute Fourier Transform
sampling_rate = 1
fft_freq = calculate_frequency(len(bpm_rr_baseline_windowed), sampling_rate)
half_point = len(fft_freq) // 2
fft_freq_half = fft_freq[:half_point]
fft_result_half = fft_result[:half_point]

# Subset 50-100
n_subset1 = n[50:100]
bpm_rr_baseline_subset1 = bpm_rr_baseline[50:100]
M1 = len(bpm_rr_baseline_subset1) -1
hamming_window1 = np.zeros(M1+1)
for i in range(M1+1):
    hamming_window1[i] = 0.54 - 0.46 * np.cos(2 * np.pi * i /M1 )
bpm_rr_baseline_windowed1 = bpm_rr_baseline_subset1 * hamming_window1
fft_result1 = fourier_transform(bpm_rr_baseline_windowed1)
fft_freq1 = calculate_frequency(len(bpm_rr_baseline_windowed1), sampling_rate)
half_point1 = len(fft_freq1) // 2
fft_freq_half1 = fft_freq1[:half_point1]
fft_result_half1 = fft_result1[:half_point1]

# Subset 100-150
n_subset2 = n[100:150]
bpm_rr_baseline_subset2 = bpm_rr_baseline[100:150]

M2 = len(bpm_rr_baseline_subset2) - 1
hamming_window2 = np.zeros(M2+1)
for i in range(M2+1):
    hamming_window2[i] = 0.54 - 0.46 * np.cos(2 * np.pi * i / M2)
bpm_rr_baseline_windowed2 = bpm_rr_baseline_subset2 * hamming_window2

fft_result2 = fourier_transform(bpm_rr_baseline_windowed2)
fft_freq2 = calculate_frequency(len(bpm_rr_baseline_windowed2), sampling_rate)
half_point2 = len(fft_freq2) // 2
fft_freq_half2 = fft_freq2[:half_point2]
fft_result_half2 = fft_result2[:half_point2]

# Subset 150-200
n_subset3 = n[150:200]
bpm_rr_baseline_subset3 = bpm_rr_baseline[150:200]

M3 = len(bpm_rr_baseline_subset3) - 1
hamming_window3 = np.zeros(M3+1)
for i in range(M3+1):
    hamming_window3[i] = 0.54 - 0.46 * np.cos(2 * np.pi * i / M3)
bpm_rr_baseline_windowed3 = bpm_rr_baseline_subset3 * hamming_window3

fft_result3 = fourier_transform(bpm_rr_baseline_windowed3)
fft_freq3 = calculate_frequency(len(bpm_rr_baseline_windowed3), sampling_rate)
half_point3 = len(fft_freq3) // 2
fft_freq_half3 = fft_freq3[:half_point3]
fft_result_half3 = fft_result3[:half_point3]

# Subset 200-250
n_subset4 = n[201:251]
bpm_rr_baseline_subset4 = bpm_rr_baseline[201:251]

M4 = len(bpm_rr_baseline_subset4) - 1
hamming_window4 = np.zeros(M4+1)
for i in range(M4+1):
    hamming_window4[i] = 0.54 - 0.46 * np.cos(2 * np.pi * i / M4)
bpm_rr_baseline_windowed4 = bpm_rr_baseline_subset4 * hamming_window4

fft_result4 = fourier_transform(bpm_rr_baseline_windowed4)
fft_freq4 = calculate_frequency(len(bpm_rr_baseline_windowed4), sampling_rate)
half_point4 = len(fft_freq4) // 2
fft_freq_half4 = fft_freq4[:half_point4]
fft_result_half4 = fft_result4[:half_point4]

# Subset 250-300
n_subset5 = n[251:301]
bpm_rr_baseline_subset5 = bpm_rr_baseline[251:301]

M5 = len(bpm_rr_baseline_subset5) - 1
hamming_window5 = np.zeros(M5+1)
for i in range(M5+1):
    hamming_window5[i] = 0.54 - 0.46 * np.cos(2 * np.pi * i / M5)
bpm_rr_baseline_windowed5 = bpm_rr_baseline_subset5 * hamming_window5

fft_result5 = fourier_transform(bpm_rr_baseline_windowed5)
fft_freq5 = calculate_frequency(len(bpm_rr_baseline_windowed5), sampling_rate)
half_point5 = len(fft_freq5) // 2
fft_freq_half5 = fft_freq5[:half_point5]
fft_result_half5 = fft_result5[:half_point5]

# Subset 300-350
n_subset6 = n[300:350]
bpm_rr_baseline_subset6 = bpm_rr_baseline[300:350]

M6 = len(bpm_rr_baseline_subset6) - 1
hamming_window6 = np.zeros(M6+1)
for i in range(M6+1):
    hamming_window6[i] = 0.54 - 0.46 * np.cos(2 * np.pi * i / M6)
bpm_rr_baseline_windowed6 = bpm_rr_baseline_subset6 * hamming_window6

fft_result6 = fourier_transform(bpm_rr_baseline_windowed6)
fft_freq6 = calculate_frequency(len(bpm_rr_baseline_windowed6), sampling_rate)
half_point6 = len(fft_freq6) // 2
fft_freq_half6 = fft_freq6[:half_point6]
fft_result_half6 = fft_result6[:half_point6]

#FFT TOTAL
FFT_TOTAL = (fft_result + fft_result1 +  fft_result2 + fft_result3 + fft_result4 + fft_result5 + fft_result6)/7
FFT_FREQ_TOTAL = (fft_freq + fft_freq1 + fft_freq2 +  fft_freq3 +  fft_freq4 +  fft_freq5 +  fft_freq6)/7

half_point_total = len(FFT_FREQ_TOTAL) // 2
fft_freq_total = FFT_FREQ_TOTAL[:half_point_total]
fft_result_total = FFT_TOTAL[:half_point_total]

#FFT SPECTRUM
def manual_interpolation(x, xp, fp):
    return np.interp(x, xp, fp)

# Frequency bands
#x_ulf = np.linspace(0.0001, 0.003, 99)
x_vlf = np.linspace(0.003, 0.04, 99)
x_lf = np.linspace(0.04, 0.15, 99)
x_hf = np.linspace(0.15, 0.4, 99)

# Compute the interpolated values manually
#y_ulf = manual_interpolation(x_ulf, fft_freq_total, np.abs(fft_result_total))
y_vlf = manual_interpolation(x_vlf, fft_freq_total, np.abs(fft_result_total))
y_lf = manual_interpolation(x_lf, fft_freq_total, np.abs(fft_result_total))
y_hf = manual_interpolation(x_hf, fft_freq_total, np.abs(fft_result_total))






#DISPLAY STREAMLIT
st.set_page_config(
    page_title="Tugas 1 ASN",
    page_icon="📊",
)

with st.sidebar:
    selected = option_menu("TUGAS 1", ["Home", "Signal Processing","HRV Analysis","DWT"], default_index=0)

if selected == "Home":
   st.title('Project ASN Kelompok 6')
   st.subheader("Anggota kelompok")
   new_title = '<p style="font-family:Georgia; color: black; font-size: 20px;">Afifah Hasnia Nur Rosita - 5023211007</p>'
   st.markdown(new_title, unsafe_allow_html=True)
   new_title = '<p style="font-family:Georgia; color: black; font-size: 20px;">Syahdifa Aisyah Qurrata Ayun - 5023211032</p>'
   st.markdown(new_title, unsafe_allow_html=True)
   new_title = '<p style="font-family:Georgia; color: black; font-size: 20px;">Sharfina Nabila Larasati - 5023211055</p>'
   st.markdown(new_title, unsafe_allow_html=True)

   

if selected == "Signal Processing":
    selected1 = option_menu(None, ["Data & Graphic", "Filter","Method & Calculation"], 
    menu_icon="cast", default_index=0, orientation="horizontal")
    
    if selected1 == 'Data & Graphic':
        st.title('Data & Graphic Input')
        st.header("Data Input")
        st.write(data)

        # Create the figure with Plotly
        fig = go.Figure(data=go.Scatter(x=x[0:2000], y=y[0:2000], mode='lines'))
        fig.update_layout(
            title="Original Signal",
            xaxis_title="Elapsed Time",
            yaxis_title="Amplitude (mV)",
            xaxis=dict(showline=True, showgrid=True),
            yaxis=dict(showline=True, showgrid=True),
            
        )
        
        # Display the figure in Streamlit
        st.header("Graphic Input")
        st.plotly_chart(fig) 
        
        
        new_title = '<p style="font-family:Georgia; color: black; font-size: 20px;">Nilai FS</p>'
        st.markdown(new_title, unsafe_allow_html=True)
        st.write(fs)
        new_title = '<p style="font-family:Georgia; color: black; font-size: 20px;">Jumlah Semua Data</p>'
        st.markdown(new_title, unsafe_allow_html=True)
        st.write(jumlahdata)
    elif selected1 == 'Filter':
        st.header("LPF")

        fig_LPF = go.Figure(data=go.Scatter(x=x[0:2000], y=lpf_ecg[0:1000], mode='lines'))
        fig_LPF.update_layout(
            title="LPF",
            xaxis_title="Sequence (n)",
            yaxis_title="Amplitude",
            xaxis=dict(showline=True, showgrid=True),
            yaxis=dict(showline=True, showgrid=True),
            template='plotly_dark'

         )
        st.plotly_chart(fig_LPF)
    
        st.header("HPF")
        fig_HPF = go.Figure(data=go.Scatter(x=x[0:2000], y=hpf_ecg[0:1000], mode='lines'))
        fig_HPF.update_layout(
            title="HPF",
            xaxis_title="Sequence (n)",
            yaxis_title="Amplitude",
            xaxis=dict(showline=True, showgrid=True),
            yaxis=dict(showline=True, showgrid=True),
            template='plotly_dark'
         )
        st.plotly_chart(fig_HPF)
    elif selected1 == 'Method & Calculation':
     optimizer_options = ['', 'Derivative', 'Squaring', 'Moving Average', 'Thresholding','Calculation']
     selected_optimizer = st.selectbox('Method & Calculation', optimizer_options)

     if selected_optimizer == 'Derivative':
        fig_DRV = go.Figure(data=go.Scatter(x=x[9:1000], y=drv[0:1000], mode='lines'))
        fig_DRV.update_layout(
            title="DERIVATIVE",
            xaxis_title="Sequence (n)",
            yaxis_title="Amplitude",
            xaxis=dict(showline=True, showgrid=True),
            yaxis=dict(showline=True, showgrid=True),
            template='plotly_dark'
        )
        st.header("DERIVATIVE")
        st.plotly_chart(fig_DRV)
     elif selected_optimizer == 'Squaring':
        fig_sqr = go.Figure(data=go.Scatter(x=x[0:1000], y=sqr[0:1000], mode='lines'))
        fig_sqr.update_layout(
            title="SQUARING",
            xaxis_title="Sequence (n)",
            yaxis_title="Amplitude",
            xaxis=dict(showline=True, showgrid=True),
            yaxis=dict(showline=True, showgrid=True),
            template='plotly_dark'
        )
        st.header("SQUARING")
        st.plotly_chart(fig_sqr)
     elif selected_optimizer == 'Moving Average':
        fig_mav = go.Figure(data=go.Scatter(x=x[0:1000], y=mav[0:1000], mode='lines'))
        fig_mav.update_layout(
            title="MAV",
            xaxis_title="Time",
            yaxis_title="Amplitude",
            xaxis=dict(showline=True, showgrid=True),
            yaxis=dict(showline=True, showgrid=True),
            template='plotly_dark'
        )
        st.header("MAV")
        st.plotly_chart(fig_mav)
     elif selected_optimizer == 'Thresholding':
        fig = go.Figure(data=go.Scatter(x=x[0:4000], y=y[0:4000], mode='lines'))
        fig.update_layout(
            title="RAW SIGNAL",
            xaxis_title="Elapsed Time",
            yaxis_title="Amplitude (mV)",
            xaxis=dict(showline=True, showgrid=True),
            yaxis=dict(showline=True, showgrid=True),
            template='plotly_dark'
        )
        st.subheader("THRESHOLDING")
        st.plotly_chart(fig)

        fig = go.Figure(data=go.Scatter(x=x[0:4000], y=thrqrs[0:4000], mode='lines'))
        fig.update_layout(
            title="SIGNAL THRESHOLD",
            xaxis_title="Time",
            yaxis_title="Amplitude",
            xaxis=dict(showline=True, showgrid=True),
            yaxis=dict(showline=True, showgrid=True),
            template='plotly_dark'
        )
        st.plotly_chart(fig)
     elif selected_optimizer == 'Calculation':
            new_title = '<p style="font-family:Georgia; color: black; font-size: 15px;">NUMBERS OF R TO R CALCULATIONS</p>'
            st.markdown(new_title, unsafe_allow_html=True)          
            st.write(ptp)
            new_title = '<p style="font-family:Georgia; color: black; font-size: 15px;">CALCULATION OF THE AMOUNT OF R</p>'
            st.markdown(new_title, unsafe_allow_html=True)
            st.write(j)
            new_title = '<p style="font-family:Georgia; color: black; font-size: 15px;">BPM CALCULATIONS</p>'
            st.markdown(new_title, unsafe_allow_html=True)
            st.write(rata)
if selected == "HRV Analysis":
    sub_selected = st.sidebar.radio(
        "Pilih Metode HRV Analysis",
        ["Time Domain Analysis", "Frequency Domain Analysis", "Non Liniear Analysis"],
        index=0
    )

    if sub_selected == 'Time Domain Analysis':
        new_title = '<p style="font-family:Georgia; color:black; font-size: 25px; text-align: center;">Time Domain Analysis</p>'
        st.markdown(new_title, unsafe_allow_html=True)
        optimizer_options1 = ['','SDNN', 'RMSSD', "pNN50", "SDSD"]
        selected_optimizer1 = st.selectbox('Time-domain analysis', optimizer_options1)

        if selected_optimizer1 == 'SDNN':
            st.write(SDNN)
        elif selected_optimizer1 == 'RMSSD':
            st.write(RMSSD)
        elif selected_optimizer1 == 'pNN50':
            st.write(pNN50)
        elif selected_optimizer1 == 'SDSD':
            st.write(SDSD)
        ########
        fig_Tachogram = go.Figure(data=go.Scatter(x=n, y=bpm_rr, mode='lines'))
        fig_Tachogram.update_layout(
            title="TACHOGRAM",
            xaxis_title="n",
            yaxis_title="BPM",
            xaxis=dict(showline=True, showgrid=True),
            yaxis=dict(showline=True, showgrid=True)
        )
        st.plotly_chart(fig_Tachogram)

        fig_histogram = go.Figure(data=go.Histogram(x=bpm_rr, nbinsx=ptp))

        fig_histogram.update_layout(
            title="Histogram Interval RR",
            xaxis_title="Interval RR",
            yaxis_title="Banyak Data",
            xaxis=dict(showline=True, showgrid=True),
            yaxis=dict(showline=True, showgrid=True),
            bargap=0.2,  # Optional: Adjusts the gap between bars
            bargroupgap=0.1,  # Optional: Adjusts the gap between groups
        )

        st.plotly_chart(fig_histogram)
    
    elif sub_selected == 'Frequency Domain Analysis':
       new_title = '<p style="font-family:Georgia; color:black; font-size: 25px; text-align: center;">Frequency Domain Analysis</p>'
       st.markdown(new_title, unsafe_allow_html=True)
       selected2 = option_menu(None, ["Tachogram", "Segmentation","Spectrum"], 
       menu_icon="cast", default_index=0, orientation="horizontal")

       if selected2 == 'Tachogram':
          st.title('Tachogram BPM Baseline')
          
          # Plotting dengan Plotly
          n = np.arange(0, ptp, 1, dtype=int)
          fig = go.Figure(data=go.Scatter(x=n, y=bpm_rr_baseline, mode='lines'))
          fig.update_layout(
            title="TACHOGRAM",
            xaxis_title="n",
            yaxis_title="BPM",
            xaxis=dict(showline=True, showgrid=True),
            yaxis=dict(showline=True, showgrid=True)
          )
          # Display the figure in Streamlit
          st.plotly_chart(fig) 
        
       elif selected2 == 'Segmentation':
         optimizer_options = ['Subset 0-50', 'Subset 50-100', 'Subset 100-150', 'Subset 150-200', 'Subset 200-250','Subset 250-300', 'Subset 300-350', 'FFT TOTAL']
         selected_optimizer = st.selectbox('Segmentation', optimizer_options)

         if selected_optimizer == 'Subset 0-50':
          st.title('Subset 0-50')
          #segmentation data
          st.header("Tachogram BPM Baseline (Segment 0-50)")
          fig_segment = go.Figure(data=go.Scatter(x=n_subset, y=bpm_rr_baseline_subset, mode='lines'))
          fig_segment.update_layout(
             title="Segmented TACHOGRAM (0-50)",
             xaxis_title="n",
             yaxis_title="BPM",
             xaxis=dict(showline=True, showgrid=True),
             yaxis=dict(showline=True, showgrid=True)
             )
          st.plotly_chart(fig_segment)

          #hamming window
          st.header("Tachogram (Subset 0-50) with Hamming Window")
          fig_segment = go.Figure(data=go.Scatter(x=n_subset, y=bpm_rr_baseline_windowed, mode='lines'))
          fig_segment.update_layout(
             title="WINDOWED TACHOGRAM (0-50)",
             xaxis_title="n",
             yaxis_title="BPM",
             xaxis=dict(showline=True, showgrid=True),
             yaxis=dict(showline=True, showgrid=True)
             )
          st.plotly_chart(fig_segment)

          #FFT 
          st.header("FFT of TACHOGRAM (Subset 0-50)")
          fig_segment = go.Figure(data=go.Scatter(x=fft_freq_half, y=np.abs(fft_result_half), mode='lines'))
          fig_segment.update_layout(
             title="FFT 0-50",
             xaxis_title="Frequency (Hz)",
             yaxis_title="Magnitude",
             xaxis=dict(showline=True, showgrid=True),
             yaxis=dict(showline=True, showgrid=True)
             )
          st.plotly_chart(fig_segment)

         elif selected_optimizer == 'Subset 50-100':
          st.title('Subset 50-100')
          #segmentation data
          st.header("Tachogram BPM Baseline (Segment 50-100)")
          fig_segment = go.Figure(data=go.Scatter(x=n_subset1, y=bpm_rr_baseline_subset1, mode='lines'))
          fig_segment.update_layout(
             title="Segmented TACHOGRAM (50-100)",
             xaxis_title="n",
             yaxis_title="BPM",
             xaxis=dict(showline=True, showgrid=True),
             yaxis=dict(showline=True, showgrid=True)
             )
          st.plotly_chart(fig_segment)

          #hamming window
          st.header("Tachogram (Subset 50-100) with Hamming Window")
          fig_segment = go.Figure(data=go.Scatter(x=n_subset1, y=bpm_rr_baseline_windowed1, mode='lines'))
          fig_segment.update_layout(
             title="WINDOWED TACHOGRAM (50-100)",
             xaxis_title="n",
             yaxis_title="BPM",
             xaxis=dict(showline=True, showgrid=True),
             yaxis=dict(showline=True, showgrid=True)
             )
          st.plotly_chart(fig_segment)

          #FFT 
          st.header("FFT of TACHOGRAM (Subset 50-100)")
          fig_segment = go.Figure(data=go.Scatter(x=fft_freq_half1, y=np.abs(fft_result_half1), mode='lines'))
          fig_segment.update_layout(
             title="FFT 50-100",
             xaxis_title="Frequency (Hz)",
             yaxis_title="Magnitude",
             xaxis=dict(showline=True, showgrid=True),
             yaxis=dict(showline=True, showgrid=True)
             )
          st.plotly_chart(fig_segment)

         elif selected_optimizer == 'Subset 100-150':
          st.title('Subset 100-150')
          #segmentation data
          st.header("Tachogram BPM Baseline (Segment 100-150)")
          fig_segment = go.Figure(data=go.Scatter(x=n_subset2, y=bpm_rr_baseline_subset2, mode='lines'))
          fig_segment.update_layout(
             title="Segmented TACHOGRAM (100-150)",
             xaxis_title="n",
             yaxis_title="BPM",
             xaxis=dict(showline=True, showgrid=True),
             yaxis=dict(showline=True, showgrid=True)
             )
          st.plotly_chart(fig_segment)

          #hamming window
          st.header("Tachogram (Subset 100-150) with Hamming Window")
          fig_segment = go.Figure(data=go.Scatter(x=n_subset2, y=bpm_rr_baseline_windowed2, mode='lines'))
          fig_segment.update_layout(
             title="WINDOWED TACHOGRAM (100-150)",
             xaxis_title="n",
             yaxis_title="BPM",
             xaxis=dict(showline=True, showgrid=True),
             yaxis=dict(showline=True, showgrid=True)
             )
          st.plotly_chart(fig_segment)

          #FFT 
          st.header("FFT of TACHOGRAM (Subset 100-150)")
          fig_segment = go.Figure(data=go.Scatter(x=fft_freq_half2, y=np.abs(fft_result_half2), mode='lines'))
          fig_segment.update_layout(
             title="FFT 100-150",
             xaxis_title="Frequency (Hz)",
             yaxis_title="Magnitude",
             xaxis=dict(showline=True, showgrid=True),
             yaxis=dict(showline=True, showgrid=True)
             )
          st.plotly_chart(fig_segment)
        
         elif selected_optimizer == 'Subset 150-200':
          st.title('Subset 150-200')
          #segmentation data
          st.header("Tachogram BPM Baseline (Segment 150-200)")
          fig_segment = go.Figure(data=go.Scatter(x=n_subset3, y=bpm_rr_baseline_subset3, mode='lines'))
          fig_segment.update_layout(
             title="Segmented TACHOGRAM (150-200)",
             xaxis_title="n",
             yaxis_title="BPM",
             xaxis=dict(showline=True, showgrid=True),
             yaxis=dict(showline=True, showgrid=True)
             )
          st.plotly_chart(fig_segment)

          #hamming window
          st.header("Tachogram (Subset 150-200) with Hamming Window")
          fig_segment = go.Figure(data=go.Scatter(x=n_subset3, y=bpm_rr_baseline_windowed3, mode='lines'))
          fig_segment.update_layout(
             title="WINDOWED TACHOGRAM (150-200)",
             xaxis_title="n",
             yaxis_title="BPM",
             xaxis=dict(showline=True, showgrid=True),
             yaxis=dict(showline=True, showgrid=True)
             )
          st.plotly_chart(fig_segment)

          #FFT 
          st.header("FFT of TACHOGRAM (Subset 150-200)")
          fig_segment = go.Figure(data=go.Scatter(x=fft_freq_half3, y=np.abs(fft_result_half3), mode='lines'))
          fig_segment.update_layout(
             title="FFT 150-200",
             xaxis_title="Frequency (Hz)",
             yaxis_title="Magnitude",
             xaxis=dict(showline=True, showgrid=True),
             yaxis=dict(showline=True, showgrid=True)
             )
          st.plotly_chart(fig_segment)
    
         elif selected_optimizer == 'Subset 200-250':
          st.title('Subset 200-250')
          #segmentation data
          st.header("Tachogram BPM Baseline (Segment 200-250)")
          fig_segment = go.Figure(data=go.Scatter(x=n_subset4, y=bpm_rr_baseline_subset4, mode='lines'))
          fig_segment.update_layout(
             title="Segmented TACHOGRAM (200-250)",
             xaxis_title="n",
             yaxis_title="BPM",
             xaxis=dict(showline=True, showgrid=True),
             yaxis=dict(showline=True, showgrid=True)
             )
          st.plotly_chart(fig_segment)

          #hamming window
          st.header("Tachogram (Subset 200-250) with Hamming Window")
          fig_segment = go.Figure(data=go.Scatter(x=n_subset4, y=bpm_rr_baseline_windowed4, mode='lines'))
          fig_segment.update_layout(
             title="WINDOWED TACHOGRAM (200-250)",
             xaxis_title="n",
             yaxis_title="BPM",
             xaxis=dict(showline=True, showgrid=True),
             yaxis=dict(showline=True, showgrid=True)
             )
          st.plotly_chart(fig_segment)

          #FFT 
          st.header("FFT of TACHOGRAM (Subset 200-250)")
          fig_segment = go.Figure(data=go.Scatter(x=fft_freq_half4, y=np.abs(fft_result_half4, mode='lines')))
          fig_segment.update_layout(
             title="FFT 200-250",
             xaxis_title="Frequency (Hz)",
             yaxis_title="Magnitude",
             xaxis=dict(showline=True, showgrid=True),
             yaxis=dict(showline=True, showgrid=True)
             )
          st.plotly_chart(fig_segment)
        
         elif selected_optimizer == 'Subset 250-300':
          st.title('Subset 250-300')
          #segmentation data
          st.header("Tachogram BPM Baseline (Segment 250-300)")
          fig_segment = go.Figure(data=go.Scatter(x=n_subset5, y=bpm_rr_baseline_subset5, mode='lines'))
          fig_segment.update_layout(
             title="Segmented TACHOGRAM (250-300)",
             xaxis_title="n",
             yaxis_title="BPM",
             xaxis=dict(showline=True, showgrid=True),
             yaxis=dict(showline=True, showgrid=True)
             )
          st.plotly_chart(fig_segment)

          #hamming window
          st.header("Tachogram (Subset 250-300) with Hamming Window")
          fig_segment = go.Figure(data=go.Scatter(x=n_subset5, y=bpm_rr_baseline_windowed5, mode='lines'))
          fig_segment.update_layout(
             title="WINDOWED TACHOGRAM (250-300)",
             xaxis_title="n",
             yaxis_title="BPM",
             xaxis=dict(showline=True, showgrid=True),
             yaxis=dict(showline=True, showgrid=True)
             )
          st.plotly_chart(fig_segment)

          #FFT 
          st.header("FFT of TACHOGRAM (Subset 250-300)")
          fig_segment = go.Figure(data=go.Scatter(x=fft_freq_half5, y=np.abs(fft_result_half5), mode='lines'))
          fig_segment.update_layout(
             title="FFT 250-300",
             xaxis_title="Frequency (Hz)",
             yaxis_title="Magnitude",
             xaxis=dict(showline=True, showgrid=True),
             yaxis=dict(showline=True, showgrid=True)
             )
          st.plotly_chart(fig_segment)

         elif selected_optimizer == 'Subset 300-350':
          st.title('Subset 300-350')
          #segmentation data
          st.header("Tachogram BPM Baseline (Segment 300-350)")
          fig_segment = go.Figure(data=go.Scatter(x=n_subset6, y=bpm_rr_baseline_subset6, mode='lines'))
          fig_segment.update_layout(
             title="Segmented TACHOGRAM (300-350)",
             xaxis_title="n",
             yaxis_title="BPM",
             xaxis=dict(showline=True, showgrid=True),
             yaxis=dict(showline=True, showgrid=True)
             )
          st.plotly_chart(fig_segment)

          #hamming window
          st.header("Tachogram (Subset 300-350) with Hamming Window")
          fig_segment = go.Figure(data=go.Scatter(x=n_subset6, y=bpm_rr_baseline_windowed6, mode='lines'))
          fig_segment.update_layout(
             title="WINDOWED TACHOGRAM (300-350)",
             xaxis_title="n",
             yaxis_title="BPM",
             xaxis=dict(showline=True, showgrid=True),
             yaxis=dict(showline=True, showgrid=True)
             )
          st.plotly_chart(fig_segment)

          #FFT 
          st.header("FFT of TACHOGRAM (Subset 300-350)")
          fig_segment = go.Figure(data=go.Scatter(x=fft_freq_half6, y=np.abs(fft_result_half6), mode='lines'))
          fig_segment.update_layout(
             title="FFT 300-350",
             xaxis_title="Frequency (Hz)",
             yaxis_title="Magnitude",
             xaxis=dict(showline=True, showgrid=True),
             yaxis=dict(showline=True, showgrid=True)
             )
          st.plotly_chart(fig_segment)

        ##FFT TOTAL
         elif selected_optimizer == 'FFT TOTAL':
          st.title('FFT TOTAL')
          #segmentation data
          st.header("Tachogram BPM Baseline (Segment 300-350)")
          fig_segment = go.Figure(data=go.Scatter(x=fft_freq_total, y=np.abs(fft_result_total), mode='lines'))
          fig_segment.update_layout(
             title="FFT TOTAL of TACHOGRAM",
             xaxis_title="Frequency (Hz)",
             yaxis_title="Magnitude",
             xaxis=dict(showline=True, showgrid=True),
             yaxis=dict(showline=True, showgrid=True)
             )
          st.plotly_chart(fig_segment)

       if selected2 == 'Spectrum':
        st.title('FFT Spectrum (Welchs periodogram)')
          
        #segmentation data
        st.header("FFT Spectrum")
        fig_segment = go.Figure(data=go.Scatter(x=fft_freq_total, y=np.abs(fft_result_total), mode='lines', line=dict(color='black', width=0.3)))

        # Add the filled regions for different frequency bands
        fig_segment.add_trace(go.Scatter(x=x_vlf, y=y_vlf, fill='tozeroy', mode='none', fillcolor='rgba(166, 81, 216, 0.2)', name='VLF'))
        fig_segment.add_trace(go.Scatter(x=x_lf, y=y_lf, fill='tozeroy', mode='none', fillcolor='rgba(81, 166, 216, 0.2)', name='LF'))
        fig_segment.add_trace(go.Scatter(x=x_hf, y=y_hf, fill='tozeroy', mode='none', fillcolor='rgba(216, 166, 81, 0.2)', name='HF'))

        # Update the layout
        fig_segment.update_layout(
        title="FFT Spectrum (Welch's periodogram)",
        xaxis_title="Frequency (Hz)",
        yaxis_title="Density",
        xaxis=dict(showline=True, showgrid=True, range=[0, 0.5]),
        yaxis=dict(showline=True, showgrid=True, range=[0, max(np.abs(fft_result_total))]),
        legend=dict(title="Frequency Bands")
        )

        # Plot the figure using Streamlit
        st.plotly_chart(fig_segment)

    elif sub_selected == 'Non Liniear Analysis':
       new_title = '<p style="font-family:Georgia; color:black; font-size: 25px; text-align: center;">Non Liniear Analysis</p>'
       st.markdown(new_title, unsafe_allow_html=True) 
      
        

          





        

       

        


    
    



        
    






 


        
        






        





        
        
    
