import logging
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from .individual import Individual

logger = logging.getLogger(__name__)


class LethargusDetector:
    # foq threshod, 0.05 old way, 0.2 for 48x holes may be 0.1 is for 24x?
    # minimum duration. 1 or 2hrs?
    # minimul interval between lethargus. 6hrs? 10 hrs
    def __init__(self, foqthreshold, minduration=1.5, mininterval=10):
        self.foqthreshod = foqthreshold
        self.minduration = minduration
        self.mininterval = mininterval

    def setadata(self, ind: Individual) -> "LethargusDetector":
        self.ind = ind
        self.interval = self.ind.interval
        return self

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

    def plotlethargus(self, _ax, _finallethargusperiod):
        # ax = _fig.get_axes()[0]
        # rect = plt.Rectangle((_finallethargusperiod[0],0), _finallethargusperiod[1]-_finallethargusperiod[0],1,alpha = 0.2,hatch="/", color="blue")
        # rect = plt.Rectangle((_finallethargusperiod[0],0), _finallethargusperiod[1]-_finallethargusperiod[0],0.05, alpha = 0.2, color="blue")
        rect = plt.Rectangle(
            (_finallethargusperiod[0], 0),
            _finallethargusperiod[1] - _finallethargusperiod[0],
            1,
            alpha=0.2,
            color="gray",
        )
        _ax.add_patch(rect)
        lthresh = 0.05
        _ax.axhline(y=lthresh, color="gray", linewidth=0.5)

        return _ax

    # foq based or area rate based are defined by input arg _oescreendmatrix.
    # self.ind.fq_oescreendmatrix
    def detectlethargus(self):
        # if there is one, [[aaa, bbb]], two [[aaa,bbb],[ccc,ddd]]
        finallethargusperiods = []
        if not len(self.ind.fq_oescreendmatrix):
            logger.info("This sample doesn't match lethargus criteria")
            # dataname = thedate + "_" + genotype+ "_" + str(samplenum +1) +  "_" +str(interval)
            # fig.savefig(dataname+"foqandlethargus.png",figsize=(8,2),dpi=100)
            return finallethargusperiods
        # if ther are multiple lethargus state, choose by human eye?
        # if there are more than 2, need to check interval between them

        if len(self.ind.fq_oescreendmatrix) > 1:
            # .ginput return [()]
            """
            logger.info("please click lethargus period")
            ax.annotate("please click lethargus period",\
                        xy=(0.05, 0.8),\
                        xycoords='axes fraction', fontsize=8,\
                        horizontalalignment='left', verticalalignment='bottom')
            clickpos = plt.ginput(n=1)[0]
            logger.info(clickpos)
            for oei in self.ind.fq_oescreendmatrix:
                if oei[0] < clickpos[0] < oei[1]:
                    finallethargusperiod=oei.copy()
                    rect = plt.Rectangle((finallethargusperiod[0],0), 
                                finallethargusperiod[1]-finallethargusperiod[0],
                                1,alpha = 0.2,hatch="/", color="blue")
                    ax.add_patch(rect)
            if finallethargusperiod == []:
                logger.info("you didnt choose any candidate")
                ax.annotate("you didnt choose any candidate",\
                        xy=(0.05, 0.7),\
                        xycoords='axes fraction', fontsize=8,\
                        horizontalalignment='left', 
                        verticalalignment='bottom', color = "red")
            """

            # here is some filter to check if they are actual lethargus.
            # interval, mean foq? etc...
            # 1st, filter by rawdata.rolling > threshold 0.01 rate
            # -> seems not work well l3 lethargus tend have high rate?
            """
            rawdata = self.ind.calcrawarea(60)
            rmrawdata = rawdata.rolling(window = 60).median()
            threshold = 0.01
            for lt in self.ind.fq_oescreendmatrix:
                logger.info(lt)
                #1min running median over the thredhold
                sumoverth = sum((rmrawdata[lt[0]: lt[1]]) > threshold)
                #if overthe threshold time is more than 80 %
                if (sumoverth/(lt[1]-lt[0])) > 0.8:
                    logger.info(sumoverth/(lt[1]-lt[0]))
                    logger.info("This has high motion rate")
                else:
                    logger.info("This one seems lethargus")
                    finallethargusperiods.append(lt)
            """

        else:
            finallethargusperiods = self.ind.fq_oescreendmatrix.copy()

        logger.info(f"finallethargusperiods {finallethargusperiods}")
        # rect = plt.Rectangle((finallethargusperiod[0],0), finallethargusperiod[1]-finallethargusperiod[0],1,alpha = 0.2,hatch="/", color="blue")
        # ax.add_patch(rect)
        return finallethargusperiods

        # dataname = thedate + "_" + genotype+ "_" + str(samplenum+1) +  "_" +str(interval)
        # fig.savefig(dataname+"foqandlethargus.png",figsize=(8,2),dpi=100)

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
