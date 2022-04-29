import functools
import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Set, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from .error import IndividualIntervalNotConsistent
from .individual import Individual
from .lethargus import Lethargus
from .lethargusdetector import LethargusDetector

logger = logging.getLogger(__name__)


@dataclass
class Samplegroup:
    date: str
    groupname: str
    expnum: str
    threshold: float = 0.05
    homepath: Path = field(init=False)
    fullindlist: List[Individual] = field(init=False, default_factory=list, repr=False)
    indset: Set[int] = field(init=False, repr=False, default_factory=set)
    is_dummy: bool = field(repr=False, default=False)

    def __post_init__(self):
        self.summarydf = None
        self.fq_maxduration = None
        self.ltfoqmatrix = None
        self.ltfoqmatrixtail = None

    def __setstate__(self, kws: Dict[str, Any]) -> None:
        kws["homepath"] = Path(kws.get("homepath"))
        self.__dict__.update(kws)

    def __getstate__(self) -> Dict[str, Any]:
        kwargs = self.__dict__.copy()
        kwargs["homepath"] = os.fspath(kwargs["homepath"])
        return kwargs

    @property
    def is_valid(self) -> bool:
        return any(ind.letharguslist for ind in self.fullindlist)

    @property
    def interval(self) -> Optional[float]:
        if self.fullindlist is None or not self.fullindlist:
            return None
        return self.fullindlist[0].interval

    @property
    def label_str(self) -> str:
        return "_".join([self.groupname])

    @property
    def indlist(self) -> List[Individual]:
        return sorted(
            (ind for ind in self.fullindlist if ind.samplenum in self.indset),
            key=lambda ind: ind.samplenum,
        )

    def set_directory(self, path: str) -> "Samplegroup":
        self.homepath = Path(path)
        return self

    def processindividuals(
        self,
        *,
        interval: float,
        dataframe: pd.DataFrame,
        min_duration: float = 1.5,
        min_interval: int = 10,
    ) -> "Samplegroup":

        ld = LethargusDetector(self.threshold, min_duration, min_interval)
        dataframe.rename(
            columns={n: i + 1 for i, n in enumerate(dataframe.columns)}
        ).agg(
            self.__aggfunc_for_ind,
            date=self.date,
            groupname=self.groupname,
            expnum=self.expnum,
            interval=interval,
            ld=ld,
        )
        return self

    def __aggfunc_for_ind(
        self,
        data: pd.Series,
        *,
        date: str,
        groupname: str,
        expnum: int,
        interval: float,
        ld: LethargusDetector,
    ):
        ind = Individual(
            date,
            groupname,
            expnum,
            interval,
            samplenum=data.name,
            rawdata=data,
        )
        if self.is_dummy:
            self.append_individual(ind)
            return
        logger.info("---------------------")
        logger.info(ind.label_str)

        fig, ax = ld.setadata(ind).processdata()

        ind.savefoqincsvfile(self.homepath)

        self.append_individual(ind)

        ax.annotate(
            ind.label_str,
            xy=(0.01, 0.9),
            xycoords="axes fraction",
            fontsize=8,
            horizontalalignment="left",
            verticalalignment="bottom",
        )

        ind.saveafig(self.homepath, fig)

    def append_individual(self, _ind: Individual) -> "Samplegroup":
        self.fullindlist.append(_ind)
        if _ind.has_valid_lethargus:
            self.indset.add(_ind.samplenum)
        return self

    # assign by samplenum
    def delete_individual(self, samplenum: int) -> "Samplegroup":
        if samplenum in self.indset:
            self.indset.remove(samplenum)
        return self

    def delete_lethargus(self, _samplenum: int, lethargusnum: int) -> "Samplegroup":
        # del self.indlist[index-1]
        ind_dict = {ind.samplenum: ind for ind in self.fullindlist}
        ind = ind_dict.get(_samplenum)
        if ind is not None:
            del ind.letharguslist[lethargusnum]

        return self

    def makedf(self, *, ltlist: List[Individual] = None) -> "Samplegroup":
        self.summarydf = None

        if ltlist is None:
            ltlist = self.indlist

        if ltlist:
            self.summarydf = pd.concat(
                (ind.get_summary_df() for ind in ltlist), ignore_index=True
            )
        return self

    def calcmaxduration(self, ordinalnum=0) -> "Samplegroup":
        max_durations = [
            ind.letharguslist[0].end - ind.letharguslist[0].start + 1
            for ind in self.fullindlist
            if ind.has_valid_lethargus and len(ind.letharguslist) > ordinalnum
        ]
        self.fq_maxduration = max(max_durations, default=0)
        self.numofadequatesample = len(max_durations)
        return self

    # adataseries could be foq, arearate etc. alignhead True or tail False
    # def makeltmatrix(self, adataseries, start, end, alignhead):
    # def makeltmatrix(self, alignhead, **kwargs):
    # **kwargs could have which lethargus uses (L3 or L4..)
    def makeltmatrix(
        self, align="head", padding: float = 1.0, ordinalnum=0
    ) -> "Samplegroup":
        if self.fq_maxduration is None:
            self.calcmaxduration(ordinalnum)

        err = f"Fail to make ltmatrix, no lethargus was detected in {self.groupname}"
        if not self.numofadequatesample > 0:
            logger.error(err)
            return self

        if len(set(ind.interval for ind in self.fullindlist)) > 1:
            # the data is not consistent
            raise IndividualIntervalNotConsistent()

        ltfoqmatrix = pd.concat(
            [
                (
                    ind.letharguslist[ordinalnum].crop_and_pad_foq(
                        ind.foq,
                        padding=padding,
                        mode="both",
                        flip=(align == "tail"),
                    )
                )
                for ind in self.fullindlist
                if ind.has_valid_lethargus and len(ind.letharguslist) > ordinalnum
            ],
            axis=1,
        ).values.T

        if align == "head":
            self.ltfoqmatrix = ltfoqmatrix
        else:
            self.ltfoqmatrixtail = np.fliplr(ltfoqmatrix)

        return self

    def makeamatrix(
        self,
        datatype,
        samplenumlist=None,
        xlim=None,
        autoelim=False,
    ):
        indlist = self.fullindlist
        if xlim is None:
            matrixcol = len(indlist[0].rawdata)
            xlim = (0, matrixcol)
        matrixcol = xlim[1] - xlim[0]

        if samplenumlist is None:
            matrixrow = len(self.fullindlist)
        else:
            matrixrow = len(samplenumlist)
            indlist = []

        ematrix = np.full(shape=(matrixrow, matrixcol), fill_value=np.nan)

        for i, ind in enumerate(indlist):
            if datatype == "foq":
                data = ind.foq[xlim[0] : (xlim[1] - 1)]

            elif datatype == "area":
                tempraw = ind.calc_raw_area(60)
                data = tempraw[xlim[0] : (xlim[1] - 1)]

            # ematrix[i][xlim[0]:(xlim[1]-1)] = data
            ematrix[i][0 : len(data)] = data
            if autoelim:
                if np.mean(data) > 0.9:  # arbitrary value....
                    ematrix[i][0 : len(data)] = np.nan

        rmatrix = ematrix

        return rmatrix

    def savesummarydf(self) -> "Samplegroup":
        if self.summarydf is None or self.summarydf.empty:
            return self
        dataname = "_".join([self.groupname])
        interval = int(self.fullindlist[0].interval)
        suffix = f"_{interval}_summary.csv"

        csvpath = os.path.join(self.homepath, dataname + suffix)
        self.summarydf.to_csv(
            csvpath,
            index=False,
        )
        return self

    def pipe(self, func: Callable, *args, **kwargs) -> "Samplegroup":
        func(self, *args, **kwargs)
        return self

    def makealignfigure(
        self,
        matrix: np.ndarray,
        alignhead: bool,
        representtype: str,
    ):
        fig = plt.figure(figsize=(8, 4))
        ax = fig.add_subplot(1, 1, 1)
        ax.set_ylim(-0.1, 1)
        ax.set_xlim(0, matrix.shape[1])
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        interval = self.fullindlist[0].interval
        hrtick = np.arange(matrix.shape[1] * interval / 60 / 60).astype(int)
        ax.set_xticks(hrtick * 60 * 60 / interval)
        ax.set_xticklabels(hrtick - 1)
        if not alignhead:
            shift = matrix.shape[1] - (hrtick * 60 * 60 / interval)[-1]
            ax.set_xticks(hrtick * 60 * 60 / interval + shift)
            ax.set_xticklabels(hrtick - max(hrtick) + 1)
        for adata in matrix:
            ax.plot(adata, linestyle=":", linewidth=0.5, color="gray")

        if representtype == "mean":
            idx = np.sum(~np.isnan(matrix), axis=0) != 0
            mean = np.full(matrix.shape[1], fill_value=np.nan)
            mean[idx] = np.nanmean(matrix[:, idx], axis=0)
            sd = np.full(matrix.shape[1], fill_value=np.nan)
            sd[idx] = np.nanstd(matrix[:, idx], axis=0)
            ax.plot(mean + sd, linestyle="--", linewidth=1, color="gray")
            ax.plot(mean - sd, linestyle="--", linewidth=1, color="gray")
            ax.plot(mean, linestyle="-", linewidth=1, color="black")
        elif representtype == "median":
            idx = np.sum(~np.isnan(matrix), axis=0) != 0
            median = np.full(matrix.shape[1], fill_value=np.nan)
            median = np.nanmedian(matrix, axis=0)
            sd = np.full(matrix.shape[1], fill_value=np.nan)
            sd[idx] = np.nanstd(matrix[:, idx], axis=0)
            ax.plot(median, linestyle="-", linewidth=1, color="black")

        fig.tight_layout()
        label = f"{self.groupname} ({str(matrix.shape[0])})"
        ax.annotate(
            label,
            xy=(0.01, 0.9),
            xycoords="axes fraction",
            fontsize=8,
            horizontalalignment="left",
            verticalalignment="bottom",
        )
        return fig

    def makeallfigure(self):
        fig, axlist = self.prepallfigureformat(
            [0, len(self.indlist[0].rawdata)],
            [0, 1],
            self.indlist,
        )

        for ind, ax in zip(self.indlist, axlist):
            ax = ind.plot_foq(ax)
            if ind.has_valid_lethargus:
                final_lt = max(ind.letharguslist, key=lambda lt: lt.totalq)
                ax = ind.plot_lethargus(ax, final_lt)
        return fig

    def makeallfigurefoq(
        self,
        indlist: List[Individual],
        xlim=None,
        ylim=(0, 1.1),
        lethargus="off",
        **kwargs,
    ):
        xlim = (0, len(indlist[0].rawdata)) if xlim is None else xlim

        fig, axlist = self.prepallfigureformat(xlim, ylim, indlist, **kwargs)

        for ax, ind in zip(axlist, indlist):
            ax = ind.plot_foq(ax)
            if lethargus.lower() == "line" and ind.has_valid_lethargus:
                for lt in ind.letharguslist:
                    ax.axvline(x=lt.start, color="black")
                    ax.axvline(x=lt.end, color="black", linestyle="-")

        return fig

    def makeallfigurewithinitialq(self, indlist: List[Individual]):
        window = 60
        fig = self.prepallfigureformat(
            [0, len(indlist[0].rawdata)],
            [0, 1.1],
            indlist,
        )
        # ax = ind.plotfoq(ax)
        # ax = ind.plotlethargus(ax, ind.fq_finallethargusperiod)
        axlist = list(fig.get_axes())
        i = 0
        for ind, ax in zip(indlist, axlist):
            ax = ind.plot_foq(ax)
            # ax = ind.plotnormalizedarea(ax,window)
            # ax.plot(ind.normalizedarea.rolling(window=int(window/ind.interval), center= True).median(), color = "black", linewidth =0.5)
            ##ax = ind.plotlethargus(ax, ind.fq_finallethargusperiod)

            # initical q
            ax.axvline(x=ind.fq_onsetcandidates[0], color="black")
            initialafterlethargus = ind.fq_onsetcandidates[
                ind.fq_onsetcandidates > ind.fq_finallethargusperiod[1]
            ][0]
            ax.axvline(x=ind.fq_finallethargusperiod[1], color="black", linestyle="--")
            ax.axvline(x=initialafterlethargus, color="black")

        return fig

    def calcinitialqlatency(self):
        logger.info("initial reactivate")
        for ind in self.indlist:
            initiallatency = ind.fq_onsetcandidates[0]

            initialafterlethargus = ind.fq_onsetcandidates[
                ind.fq_onsetcandidates > ind.fq_finallethargusperiod[1]
            ][0]
            init_ = initiallatency * ind.interval
            post_ = (
                initialafterlethargus - ind.fq_finallethargusperiod[1]
            ) * ind.interval
            logger.info(f"{init_} {post_}")

    # _xlim and _ylim two vlues list, _figsize is taple
    # def prepallfigureformat(self,_xlim, _ylim, _figsize):
    # _xlim and _ylim two vlues list, _indlist could be self.indlist or fullindlist
    # kwargs: col row
    def prepallfigureformat(
        self,
        _xlim: Tuple[int, int],
        _ylim: Tuple[int, int],
        _indlist: List[Individual],
        width: float = 8,
        col: int = 1,
        row: Optional[int] = None,
        labeloff: bool = False,
        lethargus_candidates: Optional[List[Lethargus]] = None,
        figsize=None,
    ) -> Tuple[Figure, List[Axes]]:
        # fig = plt.figure(figsize = _figsize)
        row = len(_indlist) if row is None else row
        figsize = (width, row) if figsize is None else figsize

        fig = plt.figure(figsize=figsize)

        for i, ind in enumerate(_indlist):
            # ind = self.indlist[i]
            ax = fig.add_subplot(row, col, i + 1)
            # ax.set_ylim(0,1)
            # ax.set_xlim(0, len(ind.rawdata))
            ax.set_xlim(_xlim[0], _xlim[1])
            ax.set_ylim(_ylim[0], _ylim[1])
            ax.spines["right"].set_visible(False)
            ax.spines["top"].set_visible(False)
            # hrtick = np.arange(len(ind.rawdata)*ind.interval/60/60).astype(int)
            hrtick = np.arange(_xlim[1] * ind.interval / 60 / 60, dtype=int)
            ax.set_xticks(hrtick * 60 * 60 / ind.interval)
            ax.set_xticklabels([])
            if labeloff:
                ax.set_yticklabels([])
            if lethargus_candidates is not None:
                for aol in lethargus_candidates:
                    ostart = aol.start * 60 * 60 / ind.interval
                    oend = aol.end * 60 * 60 / ind.interval
                    rect = plt.Rectangle(
                        (ostart, 0),
                        oend - ostart,
                        1,
                        alpha=0.2,
                        color="blue",
                    )
                    ax.add_patch(rect)
            # ax = ind.plotfoq(ax)
            # ax = ind.plotlethargus(ax, ind.fq_finallethargusperiod)

            ax.annotate(
                ind.label_str,
                xy=(0.01, 0.9),
                xycoords="axes fraction",
                fontsize=8,
                horizontalalignment="left",
                verticalalignment="bottom",
            )

        # ax.set_xticklabels([int(x) for x in ax.get_xticks()])
        fig.tight_layout()
        ax.set_xticklabels(hrtick)
        return fig, list(fig.get_axes())

    def makeallfigureofarea(
        self,
        window: int,
        _indlist: List[Individual],
        *,
        xlim=None,
        ylim=None,
        **kwargs,
    ):
        xlim = (0, len(_indlist[0].rawdata)) if xlim is None else xlim
        ylim = (0, 1.1) if ylim is None else ylim

        fig, axlist = self.prepallfigureformat(
            xlim,
            ylim,
            _indlist,
            **kwargs,
        )
        # ax = ind.plotfoq(ax)
        # ax = ind.plotlethargus(ax, ind.fq_finallethargusperiod)
        for ax, ind in zip(axlist, _indlist):
            # ax = ind.plotnormalizedarea(ax,window)
            ax.plot(
                ind.calc_normalized_area_with_rolling_median(window),
                color="black",
                linewidth=0.5,
            )
            # ax = ind.plotlethargus(ax, ind.fq_finallethargusperiod)
        return fig

    def make_all_figure_of_area(
        self,
        window: int,
        ind_list: List[Individual],
        figure_type: str = "norm",
        *,
        xlim=None,
        ylim=None,
        **kwargs,
    ) -> Optional[Figure]:
        if not ind_list:
            logger.warn(f"No individuals for plotting: {ind_list}")
            return

        xlim = (0, len(ind_list[0].rawdata)) if xlim is None else xlim
        ylim = (0, 1.1) if ylim is None else ylim

        fig, axlist = self.prepallfigureformat(xlim, ylim, ind_list, **kwargs)

        for ax, ind in zip(axlist, ind_list):
            if figure_type.lower() == "norm":
                data = ind.calc_normalized_area_with_rolling_median(window)
            elif figure_type.lower() == "median":
                data = (
                    ind.calc_raw_area(window)
                    .rolling(window=int(window / ind.interval), center=True)
                    .median()
                )
            elif figure_type.lower() == "mean":
                data = (
                    ind.calc_raw_area(window)
                    .rolling(window=int(window / ind.interval), center=True)
                    .mean()
                )
            else:
                logger.warn(
                    f"Figure type is not supported: {figure_type}. (Only support norm, median, mean)"
                )
                return fig

            ax.plot(
                data,
                color="black",
                linewidth=0.5,
            )

        return fig

    def makeallfigureofareacropped(
        self,
        indlist: List[Individual],
        maxduration=None,
        ordinalnum=0,
    ):
        if self.fq_maxduration is None:
            self.calcmaxduration()

        maxduration = self.fq_maxduration if maxduration is None else maxduration
        # plus minus 1 h
        margin = int(60 * 60 / self.interval)
        # foq lethargus correction. it is box car 10 min. so, eliminate 10min
        corcoef = int(60 * 10 / self.interval)

        xlim = [0, maxduration + margin * 2 - corcoef * 2]
        ylim = [0, 1.2]
        fig, axlist = self.prepallfigureformat(
            xlim, ylim, indlist, maxduration / 60 / 30 * 6
        )

        for ax, ind in zip(axlist, indlist):
            # if (ind.letharguslist is not None) and (len(ind.letharguslist) > 0):
            ltstart = ind.letharguslist[ordinalnum].start
            ltend = ind.letharguslist[ordinalnum].end
            # plus minus 1 h
            # margin = int(60*60/ind.interval)
            # foq lethargus correction. it is box car 10 min. so, eliminate 10min
            # corcoef = int(60*10/ind.interval)
            plotstart = ltstart - margin + corcoef
            plotend = ltend + margin - corcoef
            # ax = ind.plotnormalizedarea(ax,window)
            tempdata = ind.normalized_area[plotstart:plotend]
            ax.plot(tempdata.values, linestyle="-", linewidth=0.5, color="gray")

            ax.plot(
                tempdata.rolling(window=int(60 / ind.interval), center=True)
                .median()
                .values,
                linestyle="-",
                linewidth=1,
                color="black",
            )

        # x tick is now 0,1,2... so start it from -1 hr.
        shiftedtick = [int(a.get_text()) - 1 for a in ax.get_xticklabels()]
        ax.set_xticklabels(shiftedtick)
        xlimit = maxduration + margin * 2 - corcoef * 2
        (maxduration / 60 / 30 * 6, len(indlist))
        fig.set_size_inches(xlimit / 60 / 60 * indlist[0].interval * 3, len(indlist))

        return fig

    def saveafig(
        self,
        figure: Figure,
        operation_name: str,
        **kwargs,
    ) -> "Samplegroup":
        suffix = f"_{operation_name}.png"
        figpath = self.homepath.joinpath(self.label_str + suffix)
        figure.savefig(figpath, **kwargs)
        return self

    def saveltmatrix(
        self,
        matrix: np.ndarray,
        operation_typename: str,
        ordinalnum: int = 0,
    ) -> "Samplegroup":
        labels = (
            self.summarydf.astype(
                {
                    "ltnum": int,
                    "date": str,
                    "groupname": str,
                    "expnum": str,
                    "samplenum": str,
                }
            )[self.summarydf.ltnum == ordinalnum][
                ["date", "groupname", "expnum", "samplenum"]
            ].T.agg(
                lambda s: "_".join(s)
            )
        ).values

        df = pd.DataFrame(matrix.T, columns=labels)

        suffix = f"_{self.interval}_{operation_typename}_df.csv"
        csvpath = self.homepath.joinpath(self.label_str + suffix)
        df.to_csv(csvpath, index=False)
        return self

    def showcandidates(self, samplenum: int):
        if samplenum > len(self.fullindlist) or samplenum < 1:
            raise IndexError()
        ind = self.fullindlist[samplenum - 1]
        fig = ind.preparefig()
        ax = ind.plotfoq(fig.get_axes()[0])
        # ax = ind.plotcandidates(ax, ind.fq_oescreendmatrix)
        for x in ind.fq_onsetcandidates:
            ax.axvline(x=x, color="red", linewidth=0.5, linestyle="--")
        for x in ind.fq_exitcandidates:
            ax.axvline(x=x, color="blue", linewidth=0.5, linestyle="--")
        self.printcandidates(ind)
        return fig

    def printcandidates(self, ind: Individual):
        start_str = ",".join(map(str, ind.fq_onsetcandidates))
        exit_str = ",".join(map(str, ind.fq_exitcandidates))
        logger.info(f"Start candidates: {start_str}")
        logger.info(f"End candidates: {exit_str}")

    def plot_foq_heatmap(
        self,
        xlim=None,
        align=False,
        cmap="binary",
        facecolor="darkkhaki",
        figsize=(5, 2),
    ) -> Figure:
        """
        PARAMS
        ::sg Samplegroup:
        ::xlim int: maximum frame number
        ::align bool: alignment the map at onset of lethargus
        ::cmap str: colormap for heatmap
        """

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 1, 1)
        if align:
            foq_matrix = self.ltfoqmatrix
        else:
            foq_matrix = self.makeamatrix("foq")
        sample_num = foq_matrix.shape[0]
        # plot heatmap
        im = ax.imshow(foq_matrix, aspect="auto", cmap=cmap, interpolation="none")
        # x_axis
        xlim = foq_matrix.shape[1] if xlim is None else int(xlim)

        # figure configuration
        xlim_hr = xlim // 1800
        xticklabels = np.linspace(0, xlim_hr, xlim_hr + 1, dtype=int)

        # x ticks
        ax.set_xlim(0, xlim)
        ax.set_xticks(np.linspace(0, xlim_hr * 1800, xlim_hr + 1))
        ax.set_xticklabels(xticklabels)
        # y ticks
        ax.set_ylim(-0.5, sample_num - 0.5)
        ax.set_yticks(np.linspace(0, sample_num - 1, sample_num // 5 + 1) - 0.5)
        yticslabels = np.linspace(0, sample_num - 1, sample_num // 5 + 1, dtype=int) + 1
        ax.set_yticklabels(yticslabels)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_ylabel("Ind.")
        # plot figure
        ax.set_title(self.groupname)
        ax.set_facecolor(facecolor)
        return fig
