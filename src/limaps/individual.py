import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List

import matplotlib.axes as axes
import matplotlib.figure as figure
import matplotlib.patches as patch
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .lethargus import Lethargus

logger = logging.getLogger(__name__)


@dataclass
class Individual:
    date: str = field()
    groupname: str = field()
    expnum: int = field()
    interval: float = field()
    samplenum: int = field()
    rawdata: pd.Series = field(repr=False)

    # post init value
    raw_df: pd.DataFrame = field(init=False, repr=False, default_factory=pd.DataFrame)

    letharguslist: List[Lethargus] = field(init=False, repr=False, default_factory=list)

    def __post_init__(self):
        self.raw_df = self.raw_df.assign(
            raw=self.rawdata,
            qaboolean=self.calc_qa(self.rawdata),
            foq=self.calc_foq(self.rawdata),
            curmaxarea=self.calc_cummax_area(),
        )

        self.curmaxarea = None

        # if there are multiple candidate, human correction required
        self.fq_qbooleanvec = None
        self.fq_onsetcandidates = None
        self.fq_exitcandidates = None
        self.fq_oescreendmatrix: List[Lethargus] = []

        # 20180226 lethargusdetector implemented
        self.letharguslist = None

        # area rate depend lethargus
        self.normalized_area = self.calc_normalize_area()
        self.area_rate = self.calc_area_rate()  # 20 min running median

    @property
    def qaboolean(self) -> pd.Series:
        return self.raw_df.get("qaboolean", None)

    @property
    def foq(self) -> pd.Series:
        return self.raw_df.get("foq", None)

    @property
    def label_str(self) -> str:
        """

        Returns:
            str: "date_groupname_expnum_samplenum"
        """
        return "_".join(
            [
                self.date,
                self.groupname,
                str(self.expnum),
                str(self.samplenum),
            ]
        )

    @property
    def has_valid_lethargus(self) -> bool:
        """Return True if lethargus list is not empty

        Returns:
            bool: letharguslist is not None and len(letharguslist) > 0
        """
        return self.letharguslist is not None and len(self.letharguslist) > 0

    #################################################################
    def calc_qa(self, _dataseries: pd.Series) -> pd.Series:
        pixthreash = 1
        return (_dataseries < pixthreash).astype("B")

    #################################################################
    # calculate fraction of q
    # _dataseries is the rawdata which consist of pixel number above threshold
    def calc_foq(self, data: pd.Series) -> pd.Series:
        """return the 10-min rolling mean value that drop the nan value

        Args:
            data (pd.Series): a time-lapsed subtraction value

        Returns:
            pd.Series: 10-min rolling mean (np.float16) dropping the nan
        """
        return (
            self.calc_qa(data)
            .rolling(window=int(600 / self.interval), center=False)
            .mean()
            .dropna()
            .reset_index(drop=True)
        ).astype(np.float16)

    #################################################################
    # area normalized with max area. not curmax
    def calc_raw_area(self, window: int) -> pd.Series:
        maxraw: float = (
            self.rawdata.rolling(
                window=int(window / self.interval),
                center=True,
            )
            .median()
            .max()
        )
        return self.rawdata / maxraw

    def calc_cummax_area(self) -> pd.Series:
        """Calculate the accumulated max from 1-min rolling median

        Returns:
            pd.Series: 1-min rolling-median accumulated max
        """
        windowsize = int(60 / self.interval)
        return (
            self.rawdata.rolling(window=windowsize, center=False)
            .median()
            .fillna(method="bfill")
            .cummax()
            .astype(np.float)
        )

    def calc_normalize_area(self) -> pd.Series:
        """Calculate the normalized rawdata via rawdata / (1min-rolling-median-cummax)

        Returns:
            pd.Series: normalized area
        """
        self.curmaxarea = self.calc_cummax_area()
        return self.rawdata / self.curmaxarea

    def calc_normalized_area_with_rolling_median(self, window: int) -> pd.Series:
        """Calculate rolling median of normalized_area

        Args:
            window (int): Rolling window for calculate median (sec)

        Returns:
            pd.Series: rolling median of normalized_area
        """
        return self.normalized_area.rolling(
            window=round(window / self.interval),
            center=True,
        ).median()

    def calc_area_rate(self) -> pd.Series:
        """20 mins rolling median of normalized_area

        Returns:
            pd.Series: 20mins of rolling median of normalized_area
        """
        return self.calc_normalized_area_with_rolling_median(1200)

    #################################################################
    def preparefig(self):
        fig = plt.figure(figsize=(8, 2))
        ax = fig.add_subplot(111)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlim(0, len(self.rawdata))
        ax.set_ylim(0, 1.1)
        ax.yaxis.set_ticks_position("left")
        ax.xaxis.set_ticks_position("bottom")

        return fig, ax

    def plotnormalizedarea(self, _ax, _window):
        _ax.plot(
            self.calc_normalized_area_with_rolling_median(_window),
            color="gray",
            linewidth=0.5,
        )
        return _ax

    def plot_foq(self, _ax):
        _ax.plot(self.foq, linestyle="-", linewidth=1, color="black")
        return _ax

    def plotlowlevel(self, _ax):
        # _ax.plot(self.normalizedarea.rolling(window=int(600/interval), center= True).median(), color = "gray", linewidth =0.5)
        _ax.plot(self.area_rate, color="gray", linewidth=0.5)

        _ax.plot(self.foq, color="black")
        # foq threashold to determine lethargus
        # lthresh = 0.05
        # ax.axhline(y = lthresh, color = "gray")
        return _ax

    def plot_candidates(self, ax, lt_checked=False) -> "Individual":
        plot_kwargs = {
            False: dict(height=0.5, alpha=0.2, color="gray"),
            True: dict(height=0.9, alpha=0.2, color="blue"),
        }
        for lt in self.fq_oescreendmatrix:
            rect = patch.Rectangle(
                (lt.start, 0),
                lt.end - lt.start,
                **plot_kwargs[lt_checked],
            )
            ax.add_patch(rect)
        return self

    def plot_lethargus(
        self,
        ax: axes.Axes,
        lethargus: Lethargus,
        threshold: float = 0.05,
    ):
        rect = patch.Rectangle(
            (lethargus.start, 0),
            lethargus.end - lethargus.start,
            1,
            alpha=0.15,
            color="blue",
        )
        ax.add_patch(rect)
        ax.axhline(
            y=threshold,
            color="red",
            linestyle=":",
            linewidth=0.5,
            alpha=0.3,
        )

        return ax

    def saveafig(self, targetdir: str, _fig: figure.Figure, opt: str = "foq"):
        filename = f"{self.label_str}_{opt}.png"
        _fig.savefig(Path(targetdir).joinpath(filename), dpi=100)

    def get_summary_df(self) -> pd.DataFrame:
        columns = [
            "date",
            "groupname",
            "expnum",
            "samplenum",
            "ltnum",
            "fqlt_start",
            "fqlt_end",
            "fqlt_duration",
            "totalq",
            "meanquiescent",
            "meanquiescentout",
            "numberofbout",
            "qmean",
            "amean",
            "qmedian",
            "amedian",
            "qfreq",
        ]
        if self.letharguslist is None or not self.letharguslist:
            return pd.DataFrame(columns=columns)
        else:
            return (
                pd.DataFrame(
                    [lt.get_measures() for lt in self.letharguslist],
                )
                .assign(
                    date=self.date,
                    groupname=self.groupname,
                    expnum=str(self.expnum),
                    samplenum=str(self.samplenum),
                )
                .rename_axis("ltnum")
                .reset_index()[columns]
            )

    def savefoqincsvfile(self, targetdir: Path) -> "Individual":
        # save the lethargus period and foq with the date_groupname_expnum_# nam)
        filename = f"{self.label_str}_foq.csv"
        filepath = Path(targetdir).joinpath(filename)
        fq_period = max(self.fq_oescreendmatrix, key=lambda lt: lt.totalq, default=None)
        fq_period_str = "{}\n{}\n"
        if fq_period is not None:
            fq_period_str = fq_period_str.format(fq_period.start, fq_period.end)
        with open(filepath, "w") as file:
            file.write(f"{self.label_str}_{self.interval}\n")
            file.write(fq_period_str)
        self.foq.dropna().astype(np.float16).to_csv(
            filepath,
            mode="a",
            header=False,
            index=False,
        )
        return self

    def screen_by_quiet_foq(self, foq_threshold: float) -> "Individual":
        # Hayashi lab critera
        # quiescence onset time was defined as a time point
        # after which the fractional quiescence remains > 0.05 for at least 60 minutes.
        # The quiescence exit time was defined as the time point
        # after which the fractional quiescence had reached < 0.1. (actually 0.05)
        self.fq_qbooleanvec = self.foq > foq_threshold
        qbooleandiff = self.fq_qbooleanvec.astype(int).diff()
        self.fq_onsetcandidates = np.where(qbooleandiff == 1)[0]
        self.fq_exitcandidates = np.where(qbooleandiff == -1)[0]
        return self

    def screen_by_margin(
        self,
        min_duration_hr: float,
        prepostmargin_hr: float = 1.0,
    ) -> "Individual":

        if self.fq_qbooleanvec is None or not len(self.fq_onsetcandidates):
            return self

        continuouslength = int(min_duration_hr * 60 * 60 / self.interval)
        # pre and post requirments could be adjustable
        prepostmargin_hr = int(prepostmargin_hr * 60 * 60 / self.interval)

        onset_candidates = []
        for oi in self.fq_onsetcandidates:
            # logger.info("oi in onsetcandidates "+str(oi))
            if oi < prepostmargin_hr:
                msg = "pre-period imaging duraion is short"
                logger.warn(msg)
                continue
            elif len(self.fq_qbooleanvec) - oi < prepostmargin_hr:
                msg = "post-period imaging duraion is short"
                logger.warn(msg)
                continue

            sumofq = np.sum(self.fq_qbooleanvec[oi : (oi + continuouslength)])
            logger.info(f"sumofq {sumofq}")
            if sumofq == continuouslength:
                msg = "suit the criteria"
                logger.info(msg)
                onset_candidates.append(oi)

        logger.info(f"temponsetcandidates {onset_candidates}")
        # exit when goes under threshold
        exit_candidates = self.fq_exitcandidates
        logger.info(f"tempexitcandidates {exit_candidates}")

        onsetexitmatrix: List[Lethargus] = []
        for onset in onset_candidates:
            # open is no foud exit yet
            logger.info(f"slice {onset}")
            if any(filter(lambda lt: lt.start < onset < lt.end, onsetexitmatrix)):
                continue

            offset = min(filter(lambda off: off > onset, exit_candidates), default=None)
            if offset is not None:
                onsetexitmatrix.append(Lethargus(self.interval, onset, offset))

        logger.info(f"onsetexitmatrix {onsetexitmatrix}")
        self.fq_oescreendmatrix = sorted(
            (lt.calc_measurements(self.qaboolean) for lt in onsetexitmatrix),
            key=lambda lt: lt.start,
        )
        if self.letharguslist is None or not self.letharguslist:
            logger.warn(f"{self.label_str}: No lethargus detected")
        return self

    def check_interval(self, min_interval_hr: float) -> bool:
        # lethargus list
        ll = self.fq_oescreendmatrix
        if len(ll) <= 1:
            return True
        min_interval_frame = min_interval_hr * 60 * 60 / self.interval
        for pre, post in zip(ll[:-1], ll[1:]):
            pre.distanceproblem = (
                abs(post.end - pre.end) < min_interval_frame * 60 * 60 / self.interval
            )
        return any((lt.distanceproblem for lt in ll))

    def filter_by_totalq(self) -> "Individual":
        filter_matrix: List[Lethargus] = []
        while len(self.fq_oescreendmatrix) > 1 and any(
            (lt.distanceproblem for lt in self.fq_oescreendmatrix)
        ):
            pre, post = self.fq_oescreendmatrix[0], self.fq_oescreendmatrix[1]
            if not pre.distanceproblem:
                filter_matrix.append(pre)
            elif pre.totalq > post.totalq:
                self.fq_oescreendmatrix[1] = pre
                self.fq_oescreendmatrix[1].distanceproblem = False
                filter_matrix.append(self.fq_oescreendmatrix[1])
            self.fq_oescreendmatrix = self.fq_oescreendmatrix[1:]

        filter_matrix.extend(self.fq_oescreendmatrix)
        self.fq_oescreendmatrix = filter_matrix
        return self

    def screen_by_activity(self, act_threshold: float = 0.75) -> "Individual":
        self.ar_qbooleanvec = self.area_rate < act_threshold
        arbooleandiff = self.ar_qbooleanvec.astype(int).diff()
        self.ar_onsetcandidates = np.where(arbooleandiff == 1)[0]
        self.ar_exitcandidates = np.where(arbooleandiff == -1)[0]
        return self

    def confirm_result(self) -> "Individual":

        # shallow copy and sorted with total quiescence
        if self.fq_oescreendmatrix is None or not self.fq_oescreendmatrix:
            self.fq_oescreendmatrix = []
        self.letharguslist = sorted(
            self.fq_oescreendmatrix,
            key=lambda lt: lt.fq_duration,
            reverse=True,
        )
        return self

    def __setstate__(self, kws: Dict[str, Any]) -> None:
        letharguslist = kws.get("letharguslist")
        self.__dict__.update(kws)
        self.__dict__["raw_df"] = pd.DataFrame()
        self.__post_init__()
        self.__dict__["letharguslist"] = letharguslist

    def __getstate__(self) -> Dict[str, Any]:
        return dict(
            date=self.date,
            groupname=self.groupname,
            expnum=self.expnum,
            interval=self.interval,
            samplenum=self.samplenum,
            rawdata=self.rawdata.astype(np.int16),
            letharguslist=self.letharguslist,
        )
