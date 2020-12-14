from django.shortcuts import render
from django.http import HttpResponseRedirect
from django.urls import reverse
from django.views.generic import TemplateView

from .hs_filtering import *
import os


# Ustalanie rozmiaru wykresów mikrostruktury
figsize_x = 10
figsize_y = 8


class LoadFile(TemplateView):
    template_name = 'load_file.html'
    context_object_name = 'bcf_files'

    def get_context_data(self):
        context = super(LoadFile, self).get_context_data()
        context['bcf_files'] = search_for_bcf_files_in_datadir(
            GlobalData.DATADIR)
        context['navbar_properties'] = get_navbar_properties()
        return context


def file_loaded(request):
    try:
        filename = request.POST['filename']
    except:
        return render(request, 'load_file.html', {
            'bcf_files': search_for_bcf_files_in_datadir(GlobalData.DATADIR),
            'navbar_properties': get_navbar_properties(),
            'error_message': "Nie wybrano żadnego pliku.",
        })
    else:
        try:
            data, data_file = load_data_from_file(filename)
        except:
            return render(request, 'load_file.html', {
                'bcf_files': search_for_bcf_files_in_datadir(GlobalData.DATADIR),
                'navbar_properties': get_navbar_properties(),
                'error_message': "Błąd; nie można otworzyć wskazanego pliku.",
            })
        else:
            init_empty_data()
            init_load_data(data, data_file)

        return HttpResponseRedirect(reverse('peaks_identification'))


def peaks_identification(request):
    if not GlobalData.data_loaded:
        return HttpResponseRedirect(reverse('load_file'))
    if request.method == 'GET':
        if GlobalData.peaks_identification_occured:
            return render(request, 'peaks.html', {
                'haadf': haadf_plot_to_html(GlobalData.haadf.data, figsize_x, figsize_y),
                'eds_summed': eds_summed_plot_to_html(),
                'eds_summed_interactive': eds_summed_interactive_plot_to_html(),
                'isig_l': GlobalData.isig_l,
                'isig_r': GlobalData.isig_r,
                'isig_l_default': GlobalData.isig_l_default,
                'isig_r_default': GlobalData.isig_r_default,
                'peaks_data': zip(GlobalData.peaks, GlobalData.peaks_kev),
                'navbar_properties': get_navbar_properties(),
            })
        else:
            return render(request, 'peaks.html', {
                'haadf': haadf_plot_to_html(GlobalData.haadf.data, figsize_x, figsize_y),
                'eds_summed': eds_summed_plot_to_html(),
                'eds_summed_interactive': eds_summed_interactive_plot_to_html(),
                'isig_l_default': GlobalData.isig_l_default,
                'isig_r_default': GlobalData.isig_r_default,
                'navbar_properties': get_navbar_properties(),
            })
    else:
        try:
            isig_l = request.POST['isig-l']
            isig_r = request.POST['isig-r']

            # Sprawdzanie poprawności podawanych danych
            if len(isig_l) == 0:
                isig_l = GlobalData.isig_l_default
            else:
                try:
                    isig_l = int(isig_l)
                except:
                    raise Exception("Podana wartość nie jest liczbą całkowitą")

            if len(isig_r) == 0:
                isig_r = GlobalData.isig_r_default
            else:
                try:
                    isig_r = int(isig_r)
                except:
                    raise Exception("Podana wartość nie jest liczbą całkowitą")

            if isig_l < GlobalData.isig_min:
                raise Exception(
                    "Podana wartość lewego przedziału mniejsza od minimalnej")

            if isig_r > GlobalData.isig_max:
                raise Exception(
                    "Podana wartość prawego przedziału większa od maksymalnej")

            if isig_l >= isig_r:
                raise Exception(
                    "Lewa wartość przedziału musi być mniejsza od prawej")

            # Identyfikacja pików
            try:
                identify_energy_peaks_from_eds(isig_l, isig_r)
            except Exception as e:
                raise Exception(e)

            return render(request, 'peaks.html', {
                'haadf': haadf_plot_to_html(GlobalData.haadf.data, figsize_x, figsize_y),
                'eds_summed': eds_summed_plot_to_html(),
                'eds_summed_interactive': eds_summed_interactive_plot_to_html(),
                'isig_l': GlobalData.isig_l,
                'isig_r': GlobalData.isig_r,
                'isig_l_default': GlobalData.isig_l_default,
                'isig_r_default': GlobalData.isig_r_default,
                'peaks_data': zip(GlobalData.peaks, GlobalData.peaks_kev),
                'navbar_properties': get_navbar_properties(),
            })
        except Exception as e:
            return render(request, 'peaks.html', {
                'haadf': haadf_plot_to_html(GlobalData.haadf.data, figsize_x, figsize_y),
                'eds_summed': eds_summed_plot_to_html(),
                'eds_summed_interactive': eds_summed_interactive_plot_to_html(),
                'isig_l_default': GlobalData.isig_l_default,
                'isig_r_default': GlobalData.isig_r_default,
                'error_message': e,
                'navbar_properties': get_navbar_properties(),
            })


def add_new_peak(request):
    if not GlobalData.data_loaded:
        return HttpResponseRedirect(reverse('load_file'))
    if not GlobalData.peaks_identification_occured:
        return HttpResponseRedirect(reverse('peaks_identification'))
    try:
        isig_l_peak = request.POST['isig-l']
        isig_r_peak = request.POST['isig-r']

        # Sprawdzanie poprawności podawanych danych
        try:
            isig_l_peak = int(isig_l_peak)
        except:
            raise Exception("Podana wartość nie jest liczbą całkowitą")

        try:
            isig_r_peak = int(isig_r_peak)
        except:
            raise Exception("Podana wartość nie jest liczbą całkowitą")

        if isig_l_peak < GlobalData.isig_min:
            raise Exception(
                "Podana wartość lewego przedziału piku mniejsza od minimalnej")

        if isig_r_peak > GlobalData.isig_max:
            raise Exception(
                "Podana wartość prawego przedziału piku większa od maksymalnej")

        if isig_l_peak >= isig_r_peak:
            raise Exception(
                "Lewa wartość przedziału musi być mniejsza od prawej")

        # Dodawanie peaku
        kev_l = GlobalData.sum_all.axes_manager[0].index2value(isig_l_peak)
        kev_r = GlobalData.sum_all.axes_manager[0].index2value(isig_r_peak)
        add_new_energy_peak(GlobalData.peaks, GlobalData.peaks_kev,
                            GlobalData.eds.axes_manager, (kev_l, kev_r))

        return HttpResponseRedirect(reverse('peaks_identification'))
    except Exception as e:
        return render(request, 'peaks.html', {
            'haadf': haadf_plot_to_html(GlobalData.haadf.data, figsize_x, figsize_y),
            'eds_summed': eds_summed_plot_to_html(),
            'eds_summed_interactive': eds_summed_interactive_plot_to_html(),
            'isig_l_default': GlobalData.isig_l_default,
            'isig_r_default': GlobalData.isig_r_default,
            'isig_l': GlobalData.isig_l,
            'isig_r': GlobalData.isig_r,
            'peaks_data': zip(GlobalData.peaks, GlobalData.peaks_kev),
            'error_message': e,
            'navbar_properties': get_navbar_properties(),
        })


def delete_peak(request):
    if not GlobalData.data_loaded:
        return HttpResponseRedirect(reverse('load_file'))
    if not GlobalData.peaks_identification_occured:
        return HttpResponseRedirect(reverse('peaks_identification'))
    try:
        peak_id = request.POST['peakId']
        try:
            delete_energy_peak(int(peak_id))
        except Exception as e:
            raise Exception(e)

        return HttpResponseRedirect(reverse('peaks_identification'))
    except Exception as e:
        return render(request, 'peaks.html', {
            'haadf': haadf_plot_to_html(GlobalData.haadf.data, figsize_x, figsize_y),
            'eds_summed': eds_summed_plot_to_html(),
            'eds_summed_interactive': eds_summed_interactive_plot_to_html(),
            'isig_l_default': GlobalData.isig_l_default,
            'isig_r_default': GlobalData.isig_r_default,
            'isig_l': GlobalData.isig_l,
            'isig_r': GlobalData.isig_r,
            'peaks_data': zip(GlobalData.peaks, GlobalData.peaks_kev),
            'error_message': e,
            'navbar_properties': get_navbar_properties(),
        })


def channels_separation(request):
    if not GlobalData.data_loaded:
        return HttpResponseRedirect(reverse('load_file'))
    if not GlobalData.peaks_identification_occured:
        return HttpResponseRedirect(reverse('peaks_identification'))
    if GlobalData.peaks_identification_occured and get_number_of_peaks() == 0:
        error_message = "Ilość pików energii równa zero. W celu dalszego przetwarzania zmień zakres wyszukiwania pików lub dodaj piki ręcznie."
        return render(request, 'peaks.html', {
            'haadf': haadf_plot_to_html(GlobalData.haadf.data, figsize_x, figsize_y),
            'eds_summed': eds_summed_plot_to_html(),
            'eds_summed_interactive': eds_summed_interactive_plot_to_html(),
            'isig_l_default': GlobalData.isig_l_default,
            'isig_r_default': GlobalData.isig_r_default,
            'isig_l': GlobalData.isig_l,
            'isig_r': GlobalData.isig_r,
            'peaks_data': zip(GlobalData.peaks, GlobalData.peaks_kev),
            'error_message': error_message,
            'navbar_properties': get_navbar_properties(),
        })
    if request.method == 'GET':
        if not GlobalData.channels_separation_occured and not GlobalData.channels_summation_occured:
            return render(request, 'channels.html', {
                'gauss_filter_coeff_default': GlobalData.gauss_filter_coeff_default,
                'navbar_properties': get_navbar_properties(),
            })
        elif GlobalData.channels_separation_occured and not GlobalData.channels_summation_occured:
            data_peaks_plots = get_data_peaks_plots(
                GlobalData.data_peaks_filtered)
            return render(request, 'channels.html', {
                'peaks_data_with_plots': zip(data_peaks_plots, GlobalData.peaks_kev, GlobalData.weighted_sum_coeffs),
                'gauss_filter_coeff': GlobalData.gauss_filter_coeff,
                'gauss_filter_coeff_default': GlobalData.gauss_filter_coeff_default,
                'navbar_properties': get_navbar_properties(),
            })
        else:
            data_peaks_plots = get_data_peaks_plots(
                GlobalData.data_peaks_filtered)
            weighted_sum_data_peaks_plot = microstructure_plot_to_html(
                GlobalData.summed_data_peaks, figsize_x, figsize_y, 1)
            return render(request, 'channels.html', {
                'peaks_data_with_plots': zip(data_peaks_plots, GlobalData.peaks_kev, GlobalData.weighted_sum_coeffs),
                'weighted_sum_data_peaks_plot': weighted_sum_data_peaks_plot,
                'gauss_filter_coeff': GlobalData.gauss_filter_coeff,
                'gauss_filter_coeff_default': GlobalData.gauss_filter_coeff_default,
                'navbar_properties': get_navbar_properties(),
            })
    else:
        try:
            gauss_filter_coeff = request.POST['gaussFiltCoeff']

            # Sprawdzanie poprawności podawanych danych
            if len(gauss_filter_coeff) == 0:
                gauss_filter_coeff = GlobalData.gauss_filter_coeff_default
            else:
                try:
                    gauss_filter_coeff = int(gauss_filter_coeff)
                except:
                    raise Exception("Podana wartość nie jest liczbą całkowitą")
            GlobalData.gauss_filter_coeff = gauss_filter_coeff

            # Pobieranie danych dla wszystkich zidentyfikowanych pików energii
            get_identified_energy_peaks_data(
                GlobalData.gauss_filter_coeff)
            data_peaks_plots = get_data_peaks_plots(
                GlobalData.data_peaks_filtered)

            return render(request, 'channels.html', {
                'peaks_data_with_plots': zip(data_peaks_plots, GlobalData.peaks_kev, GlobalData.weighted_sum_coeffs),
                'gauss_filter_coeff': GlobalData.gauss_filter_coeff,
                'gauss_filter_coeff_default': GlobalData.gauss_filter_coeff_default,
                'navbar_properties': get_navbar_properties(),
            })
        except Exception as e:
            return render(request, 'channels.html', {
                'error_message': e,
                'gauss_filter_coeff_default': GlobalData.gauss_filter_coeff_default,
                'navbar_properties': get_navbar_properties(),
            })


def weighted_sum_of_multiplied_channels(request):
    if not GlobalData.data_loaded:
        return HttpResponseRedirect(reverse('load_file'))
    if not GlobalData.peaks_identification_occured:
        return HttpResponseRedirect(reverse('peaks_identification'))
    if not GlobalData.channels_separation_occured:
        return HttpResponseRedirect(reverse('channels_separation'))
    try:
        # przetwarzanie wag sumy
        coeffs = []
        for i in range(0, len(GlobalData.weighted_sum_coeffs)):
            coeff = request.POST['peak{}_weight'.format(i)]
            if len(coeff) == 0:
                coeff = GlobalData.weighted_sum_coeffs[i]
            else:
                try:
                    coeff = float(coeff)
                except:
                    raise Exception(
                        "Podana wartość wagi nie jest liczbą zmiennoprzecinkową")
            coeffs.append(coeff)
        GlobalData.weighted_sum_coeffs = coeffs

        # Tworzenie sumy ważonej z wybranych zidentyfikowanych pików energii
        get_energy_peaks_data_weighted_sum(
            GlobalData.data_peaks, GlobalData.weighted_sum_coeffs, GlobalData.gauss_filter_coeff)

        return HttpResponseRedirect(reverse('channels_separation'))
    except Exception as e:
        data_peaks_plots = get_data_peaks_plots(GlobalData.data_peaks_filtered)
        return render(request, 'channels.html', {
            'peaks_data_with_plots': zip(data_peaks_plots, GlobalData.peaks_kev, GlobalData.weighted_sum_coeffs),
            'gauss_filter_coeff': GlobalData.gauss_filter_coeff,
            'gauss_filter_coeff_default': GlobalData.gauss_filter_coeff_default,
            'error_message': e,
            'navbar_properties': get_navbar_properties(),
        })


def summary(request):
    if not GlobalData.data_loaded:
        return HttpResponseRedirect(reverse('load_file'))
    if not GlobalData.peaks_identification_occured:
        return HttpResponseRedirect(reverse('peaks_identification'))
    if not GlobalData.channels_separation_occured or not GlobalData.channels_summation_occured:
        return HttpResponseRedirect(reverse('channels_separation'))
    if request.method == 'GET':
        if GlobalData.summary_occured:
            return render(request, 'summary.html', {
                'labeled_data_plot': microstructure_plot_to_html(GlobalData.labeled_data, figsize_x, figsize_y, 0),
                'summed_data_plot': haadf_plot_to_html(GlobalData.summed_data_peaks, figsize_x, figsize_y),
                'summed_data_boundaries_plot': haadf_plot_to_html(GlobalData.summed_data_peaks_with_boundaries, figsize_x, figsize_y),
                'parameters': [GlobalData.area_threshold, GlobalData.dilation_range, GlobalData.closing_range],
                'default_parameters': [GlobalData.area_threshold_default, GlobalData.dilation_range_default, GlobalData.closing_range_default],
                'navbar_properties': get_navbar_properties(),
            })
        else:
            return render(request, 'summary.html', {
                'default_parameters': [GlobalData.area_threshold_default, GlobalData.dilation_range_default, GlobalData.closing_range_default],
                'navbar_properties': get_navbar_properties(),
            })
    else:
        try:
            area_threshold = request.POST['areaThreshold']
            dilation_range = request.POST['dilationRange']
            closing_range = request.POST['closingRange']

            # Sprawdzanie poprawności podawanych danych
            if len(area_threshold) == 0:
                area_threshold = GlobalData.area_threshold_default
            else:
                try:
                    area_threshold = int(area_threshold)
                except:
                    raise Exception(
                        "Podana wartość progu powierzchni obrazu nie jest liczbą całkowitą")

            if len(dilation_range) == 0:
                dilation_range = GlobalData.dilation_range_default
            else:
                try:
                    dilation_range = int(dilation_range)
                except:
                    raise Exception(
                        "Podana wartość zakresu rozszerzania jasnych obszarów nie jest liczbą całkowitą")

            if len(closing_range) == 0:
                closing_range = GlobalData.closing_range_default
            else:
                try:
                    closing_range = int(closing_range)
                except:
                    raise Exception(
                        "Podana wartość zakresu domykania jasnych obrazów nie jest liczbą całkowitą")

            get_labeled_data_with_boundaries(
                area_threshold, dilation_range, closing_range)

            GlobalData.area_threshold = area_threshold
            GlobalData.dilation_range = dilation_range
            GlobalData.closing_range = closing_range

            return render(request, 'summary.html', {
                'labeled_data_plot': microstructure_plot_to_html(GlobalData.labeled_data, figsize_x, figsize_y, 0),
                'summed_data_plot': haadf_plot_to_html(GlobalData.summed_data_peaks, figsize_x, figsize_y),
                'summed_data_boundaries_plot': haadf_plot_to_html(GlobalData.summed_data_peaks_with_boundaries, figsize_y, figsize_y, False),
                'parameters': [GlobalData.area_threshold, GlobalData.dilation_range, GlobalData.closing_range],
                'default_parameters': [GlobalData.area_threshold_default, GlobalData.dilation_range_default, GlobalData.closing_range_default],
                'navbar_properties': get_navbar_properties(),
            })
        except Exception as e:
            return render(request, 'summary.html', {
                'default_parameters': [GlobalData.area_threshold_default, GlobalData.dilation_range_default, GlobalData.closing_range_default],
                'navbar_properties': get_navbar_properties(),
                'error_message': e,
            })


# _______________________________________________________
# Funkcje pomocnicze
def search_for_bcf_files_in_datadir(datadir):
    for root, dirs, files in os.walk(datadir):
        bcf_files = []
        for file in files:
            if file.endswith('.bcf'):
                bcf_files.append(file)
    return bcf_files


def get_data_peaks_plots(data_peaks):
    data_peaks_plots = []
    for data_peak in data_peaks:
        data_peaks_plots.append(microstructure_plot_to_html(
            data_peak, figsize_x, figsize_y, 1))
    return data_peaks_plots


def get_navbar_properties():
    return [GlobalData.data_loaded, GlobalData.peaks_identification_occured, GlobalData.channels_summation_occured]
