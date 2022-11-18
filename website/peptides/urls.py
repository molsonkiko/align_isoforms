from django.urls import path

from . import views

app_name = 'peptides'
urlpatterns = [
    path('', views.index_view, name='index'),
    path('about', views.about, name='about'),
    path('proteins/<str:acc_num>/', views.protein_view, name='proteins'),
    # path('alignments/<int:pk>/', views.AlignmentView.as_view(), name='alignments'),
    path('get_protein/', views.get_protein, name='get_protein'),
    path('alignments/<str:acc_nums>/', views.alignments_view, name='alignments'),
    path('download_alignment/<str:prots>/', views.download_alignment, name='download_alignment'),
    path('proteins/<str:acc_num>/json', views.protein_json, name='proteins_json'),
    path('proteins/json_schema', views.protein_json_schema, name = 'protein_json_schema'),
    path('peptides/', views.peptides_csv, name='peptides'),
    path('site_map', views.site_map, name='site_map'),
]