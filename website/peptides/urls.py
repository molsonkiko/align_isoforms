from django.urls import path

from . import views

app_name = 'peptides'
urlpatterns = [
    path('', views.index_view, name='index'),
    path('about', views.about_view, name='about'),
    path('alignments/<str:acc_nums>/', views.alignments_view, name='alignments'),
    path('download_alignment/<str:prots>/', views.download_alignment, name='download_alignment'),
    path('get_protein/', views.get_protein, name='get_protein'),
    path('interaction_plot/<str:acc_num>', views.interaction_plot_histograms, name='interaction_plot'),
    path('download_interaction_plot_data/<str:acc_num>', views.download_interaction_plot_data, name='download_interaction_plot_data'),
    path('peptides/', views.peptides_csv, name='peptides'),
    path('proteins/<str:acc_num>/', views.protein_view, name='proteins'),
    path('proteins/<str:acc_num>/json', views.protein_json, name='proteins_json'),
    path('proteins/json_schema', views.protein_json_schema, name = 'protein_json_schema'),
    path('request_alignment/', views.request_alignment, name='request_alignment'),
    path('site_map', views.site_map, name='site_map'),
]