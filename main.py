#!/usr/bin/env python

from limaps import Project


def create_a_project():
    proj = (
        Project(
            foqthreshold=0.2,
            colnum=8,
            rownum=6,
            grouporder="h",
        )
        .set_targetfile(
            targetfile=None  # or "/home/hayashi/Desktop/220426_remi2_N2_subtractor/220426_remi2_1_2.csv"
        )
        .set_groupnames(["N2", "empty"], [(1, 6), (7, 8)])
        .read_dataframe()
        .process_samplegroups()
    )
    proj = (
        proj.saveafig(
            proj.plot_samplegroups_grid("foq", 60, "fq_duration"),
            "gridfoq.png",
            dpi=150,
        )
        .saveafig(
            proj.plot_samplegroups_grid("area", 120, "fq_duration"),
            f"gridarea.png",
            dpi=150,
        )
        .saveafig(
            proj.dotplots("fqlt_duration"),
            "fq_duration_dotplot.png",
            dpi=150,
        )
    )

    proj.to_pickle(compression=False)


if __name__ == "__main__":
    create_a_project()
