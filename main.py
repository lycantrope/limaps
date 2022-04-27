#!/usr/bin/env python

from limaps import Project
import pandas as pd
import numpy as np

if __name__ == "__main__":
    proj = (
        Project()
        .set_targetfile("./test/210913_remi2_1_2.csv")
        .set_groupnames(["N2", "rem8"], [(1, 2), (3, 8)])
        .read_dataframe()
        .process_samplegroups()
    )
    proj = (
        proj
        .saveafig(
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

    proj.to_pickle(False)
    #proj = Project.from_pickle("./test/210913_remi2_1_2.pkl")
    # print("loaded")
    # pd.testing.assert_frame_equal(proj.data, new_proj.data)
    # pd.testing.assert_index_equal(proj.groups[0], new_proj.groups[0])
    # for org, save in zip(proj.samplegroups, new_proj.samplegroups):
    #     for ind, ind_ in zip(org.fullindlist, save.fullindlist):
    #         pd.testing.assert_series_equal(ind.rawdata, ind_.rawdata)
