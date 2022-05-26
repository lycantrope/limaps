import logging
import random
from typing import Optional, Tuple

import matplotlib.colors as mcolors
import matplotlib.patches as patch
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from .individual import Individual
from .lethargus import Lethargus

logger = logging.getLogger(__name__)


class LethargusDetector:
    def __init__(self, foqthreshold, minduration=1.5, mininterval=10):
        self.foqthreshod = foqthreshold
        self.minduration = minduration
        self.mininterval = mininterval

    def setadata(self, ind: Individual) -> "LethargusDetector":
        self.ind = ind
        self.interval = self.ind.interval
        return self

    def prepfig(
        self, figsize: Tuple[float, float] = (8, 2), **kwargs
    ) -> Tuple[Figure, Axes, Axes]:
        fig = plt.figure(figsize=figsize)
        return self.prepplot(fig, **kwargs)

    def prepplot(
        self,
        fig: Figure,
        labeloff: bool = False,
        overlay=None,
        xlim: Optional[Tuple[float, float]] = None,
        ylim=(0, 1.1),
    ) -> Tuple[Figure, Axes, Axes]:
        xlim = xlim or (0, len(self.ind.rawdata))

        ax = fig.add_subplot(1, 1, 1)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.spines["right"].set_visible(False)
        # ax.spines['top'].set_visible(False)
        ax2 = ax.twiny()
        ax2.xaxis.tick_top()
        # originalxtick = ax.get_xticklabels()
        hrtick = np.arange(xlim[1] * self.ind.interval / 60 / 60).astype(int)
        ax.set_xticks(hrtick * 60 * 60 / self.ind.interval)
        ax.set_xticklabels([])
        if labeloff:
            ax.set_yticklabels([])
        label = "_".join(
            [
                self.ind.date,
                self.ind.groupname,
                str(self.ind.expnum),
                str(self.ind.samplenum),
            ]
        )
        if overlay is not None:
            if not isinstance(overlay[0], list):
                ollist = [overlay]
            else:
                ollist = overlay
            for aol in ollist:
                ostart = aol[0] * 60 * 60 / self.ind.interval
                oend = aol[1] * 60 * 60 / self.ind.interval
                rect = plt.Rectangle(
                    (ostart, 0), oend - ostart, 1, alpha=0.2, color="blue"
                )
                ax.add_patch(rect)

        ax.annotate(
            label,
            xy=(0.01, 0.9),
            xycoords="axes fraction",
            fontsize=8,
            horizontalalignment="left",
            verticalalignment="bottom",
        )

        # ax.set_xticklabels([int(x) for x in ax.get_xticks()])
        fig.tight_layout()
        # ax.xaxis.tick_top()
        ax.set_xticklabels(hrtick)
        return fig, ax, ax2

    def processdata(self) -> Tuple[Figure, Axes]:
        # here the method to detect lethargus. may be multiple.
        fig, ax, ax2 = self.prepfig()
        # ax2.plot(self.ind.foq, linewidth =0.5,
        #        color = "black", linestyle ="-")
        ax.plot(self.ind.foq, linewidth=0.5, color="black", linestyle="-")
        rawdata = self.ind.calc_raw_area(60)
        ax.plot(
            rawdata.rolling(window=60).median(),
            linewidth=0.5,
            color="gray",
            linestyle="-",
        )
        ax2.axhline(y=self.foqthreshod, linewidth=0.5, color="black", linestyle=":")
        if self.ind.letharguslist:
            return fig, ax

        # process individual data
        self.ind = (
            self.ind.screen_by_quiet_foq(self.foqthreshod)
            .screen_by_margin(self.minduration, prepostmargin_hr=1.0)
            .plot_candidates(ax)
        )

        while not self.ind.check_interval(self.mininterval):
            self.ind.filter_by_totalq()

        # plot confirmed graph
        self.ind.plot_candidates(ax, lt_checked=True).confirm_result()

        return fig, ax

    def manual_processdata(self, minimal_interval=20) -> Tuple[Figure, Axes]:
        fig, ax, ax2 = self.prepfig()
        # ax2.plot(self.ind.foq, linewidth =0.5,
        #        color = "black", linestyle ="-")
        ax.plot(self.ind.foq, linewidth=0.5, color="black", linestyle="-")
        rawdata = self.ind.calc_raw_area(60)
        area = rawdata.rolling(window=60).median()
        ax.plot(
            area,
            linewidth=0.5,
            color="gray",
            linestyle="-",
        )
        ax2.axhline(y=self.foqthreshod, linewidth=1, color="black", linestyle=":")

        # process individual data
        self.ind = self.ind.screen_by_quiet_foq(self.foqthreshod).plot_candidates(ax)

        onset = self.ind.fq_onsetcandidates.flatten()
        offset = self.ind.fq_exitcandidates.flatten()
        while len(offset) > 1 and offset[0] < onset[0]:
            offset = offset[1:]

        while len(onset) > 1 and onset[-1] > offset[-1]:
            onset = onset[:-1]

        if not len(onset) * len(offset):
            logger.error(
                f"No onset and offset were found when screening at {self.foqthreshod}"
            )
            return fig, ax
        get_color = lambda: random.choice([*mcolors.TABLEAU_COLORS.values()])
        color = get_color()
        # 15 min interval
        interval_frame = 60 * minimal_interval / self.interval
        for i, (on, off) in enumerate(zip(onset, offset)):
            if i != 0 and on - offset[i - 1] > interval_frame:
                pre_color = color
                while pre_color == color:
                    color = get_color()
            rect = patch.Rectangle(
                (on, 0),
                off - on,
                height=0.9,
                alpha=0.2,
                color=color,
            )
            ax.add_patch(rect)

        lt = {}
        while lt.get("start") is None:
            clickpos = plt.ginput(n=1)
            if clickpos and clickpos[0]:
                idx_x, _ = clickpos[0]
                idx_pos = idx_x * len(self.ind.foq)
                lt["start"] = onset[np.argmin(np.abs(onset - idx_pos))]

        ax.axvline(lt["start"], linewidth=1.5, color="red", linestyle=":")
        while lt.get("end") is None:
            clickpos = plt.ginput(n=1)
            if clickpos and clickpos[0]:
                idx_x, _ = clickpos[0]
                idx_pos = idx_x * len(self.ind.foq)
                lt["end"] = offset[np.argmin(np.abs(offset - idx_pos))]
        ax.axvline(lt["end"], linewidth=1.5, color="red", linestyle=":")

        plt.close(fig)
        fig, ax, ax2 = self.prepfig()
        # ax2.plot(self.ind.foq, linewidth =0.5,
        #        color = "black", linestyle ="-")
        ax.plot(self.ind.foq, linewidth=0.5, color="black", linestyle="-")
        ax.plot(
            area,
            linewidth=0.5,
            color="gray",
            linestyle="-",
        )
        ax2.axhline(y=self.foqthreshod, linewidth=1, color="black", linestyle=":")
        rect = patch.Rectangle(
            (lt["start"], 0),
            lt["end"] - lt["start"],
            height=0.9,
            alpha=0.4,
            color="red",
        )
        ax.add_patch(rect)
        clickpos = None
        while not clickpos:
            clickpos = plt.ginput(n=1)
        idx_x, _ = clickpos[0]
        idx_pos = idx_x * len(self.ind.foq)
        plt.close(fig)
        if idx_pos < lt["start"] or idx_pos > lt["end"]:
            return self.manual_processdata()

        self.ind.letharguslist = [
            Lethargus(
                self.ind.interval,
                lt["start"],
                lt["end"],
            ).calc_measurements(self.ind.qaboolean)
        ]
        self.ind.manual = True
        return fig, ax
