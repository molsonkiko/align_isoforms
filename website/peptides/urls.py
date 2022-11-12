from django.urls import path

from . import views

app_name = 'peptides'
urlpatterns = [
    path('', views.IndexView.as_view(), name='index'),
    path('proteins/<str:acc_num>/', views.protein_view, name='proteins'),
    # path('alignments/<int:pk>/', views.AlignmentView.as_view(), name='alignments'),
    path('get_protein/', views.get_protein, name='get_protein'),
    path('alignments/<str:acc_nums>/', views.alignments_view, name='alignments'),
]