# Generated by Django 4.1.3 on 2022-11-15 02:43

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("peptides", "0001_initial"),
    ]

    operations = [
        migrations.AlterField(
            model_name="peptide",
            name="prot",
            field=models.CharField(max_length=10),
        ),
        migrations.AlterField(
            model_name="protein",
            name="acc_num",
            field=models.CharField(max_length=10, unique=True),
        ),
    ]
