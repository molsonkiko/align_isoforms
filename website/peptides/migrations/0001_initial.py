# Generated by Django 4.1.3 on 2022-11-14 23:47

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = []

    operations = [
        migrations.CreateModel(
            name="Alignment",
            fields=[
                (
                    "prots",
                    models.CharField(max_length=250, primary_key=True, serialize=False),
                ),
                ("alignment", models.CharField(max_length=720000)),
            ],
            options={
                "abstract": False,
            },
        ),
        migrations.CreateModel(
            name="Isoform",
            fields=[
                ("pkey", models.AutoField(primary_key=True, serialize=False)),
            ],
            options={
                "abstract": False,
            },
        ),
        migrations.CreateModel(
            name="Peptide",
            fields=[
                ("pkey", models.AutoField(primary_key=True, serialize=False)),
                ("prot", models.CharField(max_length=15)),
                ("peptide", models.CharField(max_length=10000)),
                ("location", models.IntegerField(default=-1)),
            ],
        ),
        migrations.CreateModel(
            name="Protein",
            fields=[
                ("prot_id", models.AutoField(primary_key=True, serialize=False)),
                ("acc_num", models.CharField(max_length=15, unique=True)),
                ("sequence", models.CharField(max_length=45000)),
            ],
        ),
        migrations.AddIndex(
            model_name="protein",
            index=models.Index(fields=["acc_num"], name="acc_num_idx"),
        ),
        migrations.AddIndex(
            model_name="peptide",
            index=models.Index(fields=["prot"], name="prot_idx"),
        ),
        migrations.AddField(
            model_name="isoform",
            name="prot_1",
            field=models.ForeignKey(
                on_delete=django.db.models.deletion.CASCADE,
                related_name="prot_1",
                to="peptides.protein",
            ),
        ),
        migrations.AddField(
            model_name="isoform",
            name="prot_2",
            field=models.ForeignKey(
                on_delete=django.db.models.deletion.CASCADE,
                related_name="prot_2",
                to="peptides.protein",
            ),
        ),
    ]
