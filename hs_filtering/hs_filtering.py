import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpld3
import seaborn as sns

import skimage.morphology as skimo
import skimage.measure as skim
import skimage.filters as skif
import skimage as ski
import scipy.ndimage as ndi
import scipy as sc
import scipy.signal as si
import math
import numpy as np
import hyperspy.api as hs

import urllib, base64
import io

plt.ioff()
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["figure.figsize"] = (5,5)
plt.rcParams['agg.path.chunksize'] = 10000


__all__ = [
    'GlobalData',
    'load_data_from_file',
    'init_empty_data',
    'init_load_data',
    'identify_energy_peaks_from_eds',
    'get_number_of_peaks',
    'add_new_energy_peak',
    'delete_energy_peak',
    'get_identified_energy_peaks_data',
    'get_energy_peaks_data_weighted_sum',
    'get_labeled_data_with_boundaries',
    'eds_summed_interactive_plot_to_html',
    'eds_summed_plot_to_html',
    'haadf_plot_to_html',
    'microstructure_plot_to_html',
    'plot_comparison'
]


"""Klasa zawierająca wszystkie dane programu"""
class GlobalData:
    DATADIR = "data\\"
    data_loaded = False
    peaks_identification_occured = False
    channels_separation_occured = False
    channels_summation_occured = False
    summary_occured = False


# __________________________________________________________________
# Funkcje wykorzystywane przez views.py


"""Funkcja odpowiedzialna za załadowanie danych z pliku bcf"""
def load_data_from_file(filename):
    data_file = hs.load(GlobalData.DATADIR + filename)
    data = np.array(data_file[1].data)
    return (data, data_file)


"""Funkcja inicjująca dane w GlobalData, wykorzystywana także do czyszczenia starych danych podczas ładowania nowego pliku bcf"""
def init_empty_data():
    # Dane do Identyfikacji pików
    GlobalData.isig_l = None
    GlobalData.isig_r = None
    GlobalData.signal_no_background = None
    GlobalData.peaks_identification_occured = False

    # Dane do Separacji kanałów
    GlobalData.peaks = None
    GlobalData.peaks_kev = None
    GlobalData.data_peaks = None
    GlobalData.data_peaks_filtered = None
    GlobalData.gauss_filter_coeff = None
    GlobalData.channels_separation_occured = False

    # Dane do Tworzenia ważonej sumy kanałów
    GlobalData.summed_data_peaks = None
    GlobalData.weighted_sum_coeffs = None
    GlobalData.channels_summation_occured = False

    # Dane do wyświetlania wyników
    GlobalData.area_threshold = None
    GlobalData.dilation_range = None
    GlobalData.closing_range = None
    GlobalData.labeled_data = None
    GlobalData.summed_data_peaks_with_boundaries = None
    GlobalData.boundaries = None
    GlobalData.summary_occured = False


"""Inicjalizacja danych w GlobalData na podstawie danych załadowanych z pliku"""
def init_load_data(data, data_file):
    GlobalData.data = data
    GlobalData.data_file = data_file
    GlobalData.haadf = data_file[0]
    GlobalData.eds = data_file[1]
    GlobalData.eds.axes_manager[0].name = 'x'
    GlobalData.eds.axes_manager[1].name = 'y'

    # Dane do Identyfikacji pików
    GlobalData.sum_all = GlobalData.eds.sum('x')
    GlobalData.sum_all = GlobalData.sum_all.sum('y')
    GlobalData.isig_l_default = GlobalData.sum_all.axes_manager[0].value2index(0.1)
    GlobalData.isig_r_default = GlobalData.sum_all.axes_manager[0].value2index(20.)
    GlobalData.isig_min = 0
    GlobalData.isig_max = GlobalData.sum_all.axes_manager[0].axis.size - 1

    # Dane do Separacji kanałów i Tworzenia ważonej sumy kanałów
    GlobalData.gauss_filter_coeff_default = 3

    # Dane do wyświetlania wyników
    GlobalData.area_threshold_default = 700
    GlobalData.dilation_range_default = 30
    GlobalData.closing_range_default = 30

    GlobalData.data_loaded = True


"""Tworzenie przetworzonego wykresu EDS oraz wywołanie funkcji identyfikującej piki"""
def identify_energy_peaks_from_eds(isig_l, isig_r):
    isig_l_keV = GlobalData.sum_all.axes_manager[0].axis[isig_l]
    isig_r_keV = GlobalData.sum_all.axes_manager[0].axis[isig_r]

    signal_cropped = GlobalData.sum_all.isig[isig_l_keV:isig_r_keV]
    signal_no_background = signal_cropped.remove_background(
        signal_range=(isig_l_keV, isig_r_keV), fast=False)
    peaks, peaks_kev = get_peaks(
        GlobalData.eds.axes_manager, signal_no_background)
    GlobalData.signal_no_background = signal_no_background
    GlobalData.peaks = peaks
    GlobalData.peaks_kev = peaks_kev
    GlobalData.isig_l = isig_l
    GlobalData.isig_r = isig_r

    GlobalData.peaks_identification_occured = True
    GlobalData.channels_separation_occured = False
    GlobalData.channels_summation_occured = False
    GlobalData.summary_occured = False
    init_weighted_sum_coeffs_list()


def get_number_of_peaks():
    return len(GlobalData.peaks)


"""Funkcja umożliwiająca ręczne dodawanie pików przez użytkownika"""
def add_new_energy_peak(peaks, pkev, axes_manager, pair_kev):
    v2i = axes_manager[2].value2index
    left = v2i(pair_kev[0])
    right = v2i(pair_kev[1])
    peaks.append((left, right))
    pkev.append(pair_kev)

    GlobalData.channels_separation_occured = False
    GlobalData.channels_summation_occured = False  
    GlobalData.summary_occured = False
    init_weighted_sum_coeffs_list()


"""Funkcja umożliwiająca ręczne usuwanie pików przez użytkownika"""
def delete_energy_peak(id):
    del GlobalData.peaks[id]
    del GlobalData.peaks_kev[id]

    GlobalData.channels_separation_occured = False
    GlobalData.channels_summation_occured = False  
    GlobalData.summary_occured = False
    init_weighted_sum_coeffs_list()


"""Funkcja tworząca przetworzone, przefiltrowane dane dla każdego z pików"""
def get_identified_energy_peaks_data(gauss_filter_coeff):
    data_peaks = np.array(get_energy_peaks_data(GlobalData.data, GlobalData.peaks))
    coeffs_vert = [x[0] for x in trend_coefficients(data_peaks, 0)]
    coeffs_hor = [x[0] for x in trend_coefficients(data_peaks, 1)]
    data_peaks_unwedge_v = unwedge(data_peaks, coeffs_vert, 0)
    data_peaks_unwedge = unwedge(data_peaks_unwedge_v, coeffs_hor, 1)
    data_peaks_unwedge_filtered = filter_all_peaks_gauss(data_peaks_unwedge, gauss_filter_coeff)

    GlobalData.data_peaks = data_peaks_unwedge
    GlobalData.data_peaks_filtered = data_peaks_unwedge_filtered
    GlobalData.channels_separation_occured = True
    GlobalData.channels_summation_occured = False
    GlobalData.summary_occured = False


"""Funkcja filtrująca sumę ważoną map rozkładu intensywności dla podanych pików"""
def get_energy_peaks_data_weighted_sum(data, coeff, gauss_filter_coeff):
    summed_data_peaks = energy_peaks_data_weighted_sum(data, coeff)
    summed_data_peaks_filtered = filter_one_peak_gauss(summed_data_peaks, gauss_filter_coeff)

    GlobalData.summed_data_peaks = summed_data_peaks_filtered
    GlobalData.channels_summation_occured = True 
    GlobalData.summary_occured = False

 
"""Funkcja identyfikująca granice ziaren oraz oznaczająca ziarna"""
def get_labeled_data_with_boundaries(area_threshold, dilation_range, closing_range):
    labeled_data, contours = filterset(
        GlobalData.summed_data_peaks, 0, area_threshold, dilation_range, closing_range, 1)

    GlobalData.labeled_data = labeled_data
    GlobalData.boundaries = contours
    GlobalData.summed_data_peaks_with_boundaries = GlobalData.summed_data_peaks + contours * 0.5
    GlobalData.summary_occured = True


 
# __________________________________________________________________
# Funkcje pomocnicze


"""Funkcja zwracająca zidentyfikowane piki oraz odpowiadające im przedziały w keV"""
def get_peaks(axes_manager, signal_no_background):
    import math
    # findPeaks returns 3D array, first dimension is neglected
    peaks = signal_no_background.find_peaks1D_ohaver(amp_thresh=1)[0]
    p = []
    pkev = []
    v2i = axes_manager[2].value2index
    for r in peaks:
        if(not math.isnan(r[2])):
            left = v2i(r[0]-r[2])
            right = v2i(r[0]+r[2])
            p.append((left, right))
            pkev.append((r[0] - r[2], r[0]+r[2]))
    return (p, pkev)


"""Funkcja zwracająca dane EDS dla każdego z pików energii"""
def get_energy_peaks_data(data, peaks):
    data_peaks = np.ndarray(shape=(len(peaks),data.shape[0],data.shape[1]))
    for i in range(len(peaks)):
        data_peaks[i] = data[:,:,peaks[i][0]:peaks[i][1]].sum(axis=2)
    return data_peaks


"""Funkcja zwracająca sumę ważoną map rozkładu intensywności dla podanych pików"""
def energy_peaks_data_weighted_sum(data_peaks, weighted_sum_coeffs):
    r = np.ndarray(shape=data_peaks[0].shape)
    for d, c in zip(data_peaks, weighted_sum_coeffs):
        r = np.add(r, (d*c)/d.max())
    return r


def trend_coefficients(data, axis):
    import statistics as st
    from sklearn.linear_model import LinearRegression
    coeffs = []
    if (axis == 0) or (axis == 1):
        for p in range (data.shape[0]):
            trendCoefs = []
            for i in range (0,data.shape[axis + 1],int(data.shape[axis + 1]/10)):
                if axis == 0:
                    cut = data[p][:,i]
                else:
                    cut = data[p][i,:]
                model = LinearRegression()
                x = np.array(range(len(cut))).reshape((-1, 1))
                y = np.array(cut)
                model.fit(x, y)
                pred = model.predict(x)
                trendCoefs.append(model.coef_[0])
                #print (p, i, trendCoefs)
            coeffs.append([st.mean(trendCoefs),trendCoefs])
    else:
        print("Axis can be 0 or 1")
    return coeffs


def unwedge(data, coefficients, axis):
    data_unwedged = np.ndarray(shape=data.shape)
    for datalayer, dataunwedged, coeff in zip(data, data_unwedged, coefficients):
        for i in range(datalayer.shape[axis]):
            if axis == 0:
                dataunwedged[i,:] = datalayer[i,:]-(i*coeff)
            else:
                dataunwedged[:,i] = datalayer[:,i]-(i*coeff)
    return data_unwedged  


"""Filtrowanie obrazu dla pojedynczego piku metodą Gaussa"""
def filter_one_peak_gauss(matrices, size):
    filtered = np.ndarray(shape=matrices.shape)
    filtered = ski.filters.gaussian(matrices, size)
    return filtered


"""Filtrowanie obrazu dla listy pików metodą Gaussa"""
def filter_all_peaks_gauss(matrices, size):
    filtered = np.ndarray(shape=matrices.shape)
    r = range(0,matrices.shape[0])
    for i in r:
        filtered[i] = ski.filters.gaussian(matrices[i], size)
    return filtered


"""Inicjowanie pustego wektora współczynników dla sumy ważonej kanałów"""
def init_weighted_sum_coeffs_list():
    coeffs = [0.] * len(GlobalData.peaks)
    GlobalData.weighted_sum_coeffs = coeffs


def addBoundaries(array, width=1):
    array[0:width,:]=1
    array[:,0:width]=1
    array[array.shape[0]-width:array.shape[0],:]=1
    array[:,array.shape[1]-width:array.shape[1]]=1
    return array


def filterset(signal, pics = 1, sizesToRemove = 700, dilationSquare = 30, closingSquare = 30,boundary=1):
    grayscale2 =signal
    gaussed = skif.gaussian(grayscale2,3)
    if pics: plt.imshow(gaussed, cmap='gray')
    if pics: skif.thresholding.try_all_threshold(gaussed,figsize=(25, 25))
    thresh = skif.thresholding.threshold_li(gaussed)
    binary = gaussed > thresh
    if pics: plt.imshow(binary, cmap='gray')
    
    label_objects, nb_labels = ndi.label(binary)
    sizes = np.bincount(label_objects.ravel())
    mask_sizes = sizes > sizesToRemove
    mask_sizes[0] = 0
    binary_cleaned = mask_sizes[label_objects]

    if pics: plt.imshow(binary_cleaned, cmap='gray')
    diamclo = skimo.dilation(binary_cleaned, skimo.square(dilationSquare))
    if pics: plt.imshow(diamclo, cmap='gray')
    diamclo2 = skimo.binary_closing(diamclo, skimo.square(closingSquare))
    diamclo2 = addBoundaries(diamclo2, boundary)
    if pics: plt.imshow(diamclo2, cmap='gray', interpolation='none')
    eroded = skimo.thin(diamclo2).astype(int)
    if pics: plt.imshow(eroded, cmap='gray', interpolation='none')
    erodedWide = skimo.dilation(eroded, skimo.square(5))
    if pics: plt.imshow(erodedWide, cmap='gray', interpolation='none')
    diamclo4 =  np.where(erodedWide==0, 1, erodedWide)
    diamclo4 =  np.where(erodedWide==1, 0, diamclo4)
    all_labels = skim.label(diamclo4)
    blobs_labels = skim.label(diamclo4, background=0)
    return blobs_labels, erodedWide



# __________________________________________________________________
# Funkcje do wyświetlania danych


"""Wyświetlanie interaktywnego wykresu EDS po zsumowaniu wartości z wszystkich pikseli (x, y) obrazu"""
def eds_summed_interactive_plot_to_html():
    figsizeSwitch = plt.rcParams["figure.figsize"]
    plt.rcParams["figure.figsize"] = (8,3)

    fig, ax = plt.subplots()
    ax.plot(GlobalData.sum_all)
    ax.locator_params(axis='x', nbins=10)
    plt.close()
    html = mpld3.fig_to_html(fig)

    plt.rcParams["figure.figsize"] = figsizeSwitch
    return html


"""Wyświetlanie wykresu EDS po zsumowaniu wartości z wszystkich pikseli (x, y) obrazu z zachowaniem jednostek"""
def eds_summed_plot_to_html():
    figsizeSwitch = plt.rcParams["figure.figsize"]
    plt.rcParams["figure.figsize"] = (8,5)

    scale = GlobalData.sum_all.axes_manager[0].scale
    offset = GlobalData.sum_all.axes_manager[0].offset

    from matplotlib.ticker import FuncFormatter
    fig, ax = plt.subplots()
    ax.plot(GlobalData.sum_all)
    ax.locator_params(axis='x', nbins=10)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: "%.1f" % (x * scale + offset)))
    ax.xaxis.labelpad = 10
    ax.yaxis.labelpad = 10
    ax.set_xlabel("Energy axis (keV)")
    ax.set_ylabel("X-rays (Counts)")

    buf = io.BytesIO()
    fig.savefig(buf, format='png')
    buf.seek(0)
    string = base64.b64encode(buf.read())
    uri = 'data:image/png;base64,' + urllib.parse.quote(string)
    plt.close()

    plt.rcParams["figure.figsize"] = figsizeSwitch
    return uri


"""Wyświetlanie obrazów mikrostruktury"""
def microstructure_plot_to_html(sig, figsize_x, figsize_y, cmap):
    figsizeSwitch = plt.rcParams["figure.figsize"]
    plt.rcParams["figure.figsize"] = (figsize_x, figsize_y)

    fig, ax = plt.subplots()
    if cmap == 0:
        plt.imshow(sig, cmap='nipy_spectral')
    else:
        ax = sns.heatmap(sig, xticklabels=False, yticklabels=False)
    buf = io.BytesIO()
    fig.savefig(buf, format='png')
    buf.seek(0)
    string = base64.b64encode(buf.read())
    uri = 'data:image/png;base64,' + urllib.parse.quote(string)
    plt.close()

    plt.rcParams["figure.figsize"] = figsizeSwitch
    return uri


    """Wyświetlanie obrazów HAADF"""
def haadf_plot_to_html(data, figsize_x, figsize_y, cbar_enabled=True):
    figsizeSwitch = plt.rcParams["figure.figsize"]
    plt.rcParams["figure.figsize"] = (figsize_x, figsize_y)

    scale_x = GlobalData.haadf.axes_manager[0].scale
    offset_x = GlobalData.haadf.axes_manager[0].offset
    scale_y = GlobalData.haadf.axes_manager[1].scale
    offset_y = GlobalData.haadf.axes_manager[1].offset

    from matplotlib.ticker import FuncFormatter
    fig, ax = plt.subplots()
    ax = sns.heatmap(data, cmap=plt.cm.gray, cbar=cbar_enabled)
    ax.locator_params(axis='x', nbins=10)
    ax.locator_params(axis='y', nbins=10)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: "%.1f" % (x * scale_x + offset_x)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, pos: "%.1f" % (y * scale_y + offset_y)))
    ax.xaxis.labelpad = 10
    ax.yaxis.labelpad = 10
    ax.set_xlabel("width axis (μm)")
    ax.set_ylabel("height axis (μm)")
    
    buf = io.BytesIO()
    fig.savefig(buf, format='png')
    buf.seek(0)
    string = base64.b64encode(buf.read())
    uri = 'data:image/png;base64,' + urllib.parse.quote(string)
    plt.close()

    plt.rcParams["figure.figsize"] = figsizeSwitch
    return uri


"""Wyświetlanie porównania wykresów"""
def plot_comparison(original, filtered, filter_name):
    figsizeSwitch = plt.rcParams["figure.figsize"]

    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(20, 10), sharex=True,
                                   sharey=True)
    ax1.imshow(original, cmap=plt.cm.gray)
    ax1.set_title('original')
    ax1.axis('off')
    ax2.imshow(filtered, cmap=plt.cm.gray)
    ax2.set_title(filter_name)
    ax2.axis('off')

    buf = io.BytesIO()
    fig.savefig(buf, format='png')
    buf.seek(0)
    string = base64.b64encode(buf.read())
    uri = 'data:image/png;base64,' + urllib.parse.quote(string)
    
    plt.close()

    plt.rcParams["figure.figsize"] = figsizeSwitch
    return uri
