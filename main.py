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
        .set_targetfile(targetfile="./220426_remi2_1_2.csv")
        .set_groupnames(
            uniquegroupnames=["N2"],
            gindex=[(1, 6)],
        )
        .read_dataframe()
        .process_samplegroups()
    )
    proj = (
        proj.saveafig(
            figure=proj.plot_samplegroups_grid("foq", 60, "fq_duration"),
            filename="gridfoq.png",
            dpi=150,
        )
        .saveafig(
            figure=proj.plot_samplegroups_grid("area", 120, "fq_duration"),
            filename="gridarea.png",
            dpi=150,
        )
        .saveafig(
            figure=proj.dotplots("fqlt_duration"),
            filename="fq_duration_dotplot.png",
            dpi=150,
        )
    )

    proj.to_pickle().create_summary_slide()


if __name__ == "__main__":
    create_a_project()
