# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 10:17:25 2017

@author: taizo kawano
"""

from typing import List, Tuple

import matplotlib.pyplot as plt
import matplotlib.axes as axes
import numpy as np
import pandas as pd


def get_jitter(data, bins=10):
    histdata = np.histogram(data, bins=bins)
    # array for how far shift from the middle line
    shift_val = []
    for num in histdata[0]:
        jitter = (
            np.arange(1, num + 1, dtype=int) // 2 * np.power(-1, np.arange(num))
            + (1 - num % 2) / 2
        )
        np.random.shuffle(jitter)
        shift_val.extend(jitter)
    shift_arr = np.asarray(shift_val, dtype=float)
    shift_arr += (np.random.rand(len(shift_arr)) - 0.5) / 10
    shift_arr = shift_arr / shift_arr.max()
    return shift_arr[np.argsort(np.argsort(data))]


def add_jitter_plot(
    data: np.ndarray,
    *,
    ax: axes.Axes,
    center_pos: int,
    thickness: float = 1,
    color: str = "black",
    outlier: str = "oc",
    size: float = 1,
    coeff: float = 0.3,
    jitter_bin: int = 10,
):
    q25, q50, q75 = np.percentile(data, [25, 50, 75])
    ax.hlines(
        y=q50,
        xmin=center_pos - 0.3,
        xmax=center_pos + 0.3,
        colors="black",
        linewidths=thickness,
    )
    ax.hlines(
        y=[q25, q75],
        xmin=center_pos - 0.1,
        xmax=center_pos + 0.1,
        colors="black",
        linewidths=thickness,
    )
    ax.vlines(x=center_pos, ymin=q25, ymax=q75, colors="black", linewidths=thickness)
    upperlimit = q75 + 1.5 * (q75 - q25)
    lowerlimit = q25 - 1.5 * (q75 - q25)
    if outlier == "line":
        ax.vlines(
            x=center_pos,
            ymin=lowerlimit,
            ymax=upperlimit,
            colors="gray",
            linewidths=thickness / 2,
        )

    jitter = get_jitter(data, jitter_bin)
    if outlier != "oc":
        ax.scatter(
            jitter * coeff + center_pos,
            data,
            s=sizes,
            color="black",
        )
        return

    outliner_mask = np.logical_or(
        np.array(data) > upperlimit, np.array(data) < lowerlimit
    )

    sizes = np.ones(len(data)) * size
    facecolor = np.full(len(data), fill_value=color)

    sizes[outliner_mask] *= 2
    facecolor[outliner_mask] = "none"
    ax.scatter(
        jitter * coeff + center_pos,
        data,
        s=sizes,
        color=color,
        facecolors=facecolor,
        marker="o",
    )
    return None


########################################
# dot plot function
# median with inter quatile range are shown as bars
# _thedata: list of list: [[],[]...]. not numpy array.
# ylim: (min, max)
# or melted form of pandas.Dataframe.
# if melted form, label must come first column
# **kwargs could be labels, groupnames ylim, size,thickness,
# binnum, coeff, sort, figsize
#####################################
def dotplots(
    groups: List[pd.Series],
    figsize: Tuple[int, int] = None,
    labels: List[str] = None,
    ylim: Tuple[int, int] = None,
    size: int = 1,
    thickness: float = 1,
    outlier: str = "none",
    binnum: int = 10,
    coeff: float = 0.5,
):

    # basic figure setup.
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(1, 1, 1)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.yaxis.set_ticks_position("left")
    ax.xaxis.set_ticks_position("bottom")

    if labels is not None and len(labels) != len(groups):
        raise ValueError("The number of labels is not equal to groups number")

    labelswithsize = [f"{g.name}\n{g.count()}" for g in groups]
    if labels:
        labelswithsize = [f"{l}\n{g.count()}" for l, g in zip(labels, groups)]

    ymin = min((g.min() for g in groups))
    ymax = max((g.max() for g in groups))
    if ylim is not None and len(ylim) == 2:
        ymin, ymax = ylim

    y_margin = (ymax - ymin) * 0.05  # add 5% margin on both min and max

    ax.set_ylim(ymin - y_margin, ymax + y_margin)
    ax.set_xlim(-0.5, len(groups) - 0.5)

    ax.set_xticks(np.arange(0, len(groups), dtype=int))
    ax.set_xticklabels(labels)
    ax.set_xticklabels(labelswithsize)

    for i in range(groupnumber):
        q25, q50, q75 = np.percentile(adata, [25, 50, 75])
        ax.hlines(
            y=q50,
            xmin=i - 0.3,
            xmax=i + 0.3,
            colors="black",
            linewidths=thickness,
        )
        ax.hlines(
            y=[q25, q75],
            xmin=i - 0.1,
            xmax=i + 0.1,
            colors="black",
            linewidths=thickness,
        )
        ax.vlines(x=i, ymin=q25, ymax=q75, colors="black", linewidths=thickness)
        upperlimit = q75 + 1.5 * (q75 - q25)
        lowerlimit = q25 - 1.5 * (q75 - q25)
        if outlier == "line":
            ax.vlines(
                x=i,
                ymin=lowerlimit,
                ymax=upperlimit,
                colors="gray",
                linewidths=thickness / 2,
            )

        outliner_mask = np.logical_or(
            np.array(adata) > upperlimit, np.array(adata) < lowerlimit
        )
        outdata = adata[outliner_mask]

        pdata = adata[~outliner_mask]

        outlier_params = dict(
            s=size * 2,
            color="black",
            facecolors="none",
            marker="o",
        )

        jitter = get_jitter(pdata, binnum)
        ax.scatter(
            jitter * coeff + i,
            pdata,
            s=size,
            color="black",
        )
    fig.tight_layout()

    return fig
