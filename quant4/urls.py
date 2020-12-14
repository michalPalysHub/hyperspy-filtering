from django.contrib import admin
from django.urls import path, include

from django.conf import settings
from django.conf.urls.static import static

from hs_filtering import views

urlpatterns = [
    path('admin/', admin.site.urls),
    path('', views.LoadFile.as_view(), name='load_file'),
    path('loaded/', views.file_loaded, name='file_loaded'),
    path('peaks/', views.peaks_identification, name='peaks_identification'),
    path('add/', views.add_new_peak, name='add_peak'),
    path('delete/', views.delete_peak, name='delete_peak'),
    path('channels', views.channels_separation, name="channels_separation"),
    path('sum', views.weighted_sum_of_multiplied_channels, name="sum_channels"),
    path('summary', views.summary, name="summary"),
]
