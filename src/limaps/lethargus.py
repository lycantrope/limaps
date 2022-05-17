import logging
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class Lethargus:
    interval: int = field()
    start: int = field()
    end: int = field()
    meanquiescent: float = field(init=False, default=np.nan)
    totalq: int = field(init=False, default=-1)
    numberofbout: float = field(init=False, default=np.nan)
    qmean: float = field(init=False, default=np.nan)
    amean: float = field(init=False, default=np.nan)
    qmedian: float = field(init=False, default=np.nan)
    amedian: float = field(init=False, default=np.nan)
    fq_duration: float = field(init=False, default=np.nan)
    qfreq: float = field(init=False, default=np.nan)

    qaboutdataframe: pd.DataFrame = field(
        init=False,
        repr=False,
        default_factory=pd.DataFrame,
    )
    lethargusperiodqadf: pd.DataFrame = field(
        init=False,
        repr=False,
        default_factory=pd.DataFrame,
    )
    distanceproblem: bool = field(init=False, repr=False, default=False)

    def calc_measurements(self, qaboolean: pd.Series) -> "Lethargus":
        # lethargus start and end
        ltstart = self.start
        ltend = self.end

        totalq = qaboolean[ltstart:ltend].sum() * self.interval / 60
        # totalq all imaging duration
        totalqall = qaboolean.sum() * self.interval / 60
        # totalq out of lethargus
        totalqout = totalqall - totalq
        # here cause some trouble, is arawdata was seriise totalq is not list.s
        logger.info(f"totalq {int(totalq)} min")

        # calc q and a bout duration
        qabooleandiff = qaboolean.astype(int).diff()
        qstart = np.where(qabooleandiff == 1)[0]
        qend = np.where(qabooleandiff == -1)[0]
        if not len(qstart) or not len(qend):
            return self
        # fix always qstart < qend
        if qstart[0] > qend[0]:
            qend = qend[1:]
        if qstart[-1] > qend[-1]:
            qstart = qstart[:-1]
        qaboutdataframe = pd.DataFrame().assign(
            qstart=qstart,
            qend=qend,
            qduration=qend - qstart,
            aduration=np.append(qstart[1:] - qend[:-1], np.nan),
        )
        lethargusperiodqadf = qaboutdataframe.query(
            f"qstart > {ltstart} & qend < {ltend}"
        )

        lastrow_idx = lethargusperiodqadf.tail(1).index
        lastrow = lethargusperiodqadf.loc[lastrow_idx, :]
        qend_aduration = int(lastrow["qend"]) + int(lastrow["aduration"])

        logger.info(f"qend + aduration {qend_aduration}")
        logger.info(f"ltend {ltend}")

        if qend_aduration > ltend:
            logger.info("trim the lastrow data")
            lastrow["aduration"] = np.nan

        numberofbout = len(lethargusperiodqadf)
        logger.info(f"numberofbout {numberofbout}")
        info = lethargusperiodqadf[["qduration", "aduration"]].agg(["mean", "median"])
        # sec
        qmean = info["qduration"]["mean"] * self.interval
        amean = info["aduration"]["mean"] * self.interval
        # frame count
        qmedian = info["qduration"]["median"]
        amedian = info["aduration"]["median"]

        logger.info(f"qmean {np.round(qmean, decimals=2)} sec")
        logger.info(f"amean {np.round(amean, decimals=2)} sec")

        logger.info(f"qmedian {qmedian}frame")
        logger.info(f"amedian {amedian}frame")

        # frequency /hour
        fq_duration = (ltend - ltstart) * self.interval / 60 / 60
        qfreq = numberofbout / fq_duration
        logger.info(f"qfreq {np.round(qfreq, 2)} num/hour")
        # quiescent / time neary= mean foq
        meanquiescent = totalq / ((ltend - ltstart) * self.interval / 60)
        logger.info(f"meanquiescent {meanquiescent}")
        # 180226 meanquiescentout is not correct way if multiple lethargus are deteced
        meanquiescentout = totalqout / (
            (len(qaboolean) - (ltend - ltstart)) * self.interval / 60
        )
        logger.info(f"meanquiescentout {meanquiescentout}")

        self.meanquiescent = meanquiescent
        self.meanquiescentout = meanquiescentout
        self.totalq = totalq
        self.numberofbout = numberofbout
        self.qmean = qmean
        self.amean = amean
        self.qmedian = qmedian
        self.amedian = amedian
        self.fq_duration = fq_duration
        self.qfreq = qfreq
        self.qaboutdataframe = qaboutdataframe
        self.lethargusperiodqadf = lethargusperiodqadf

        return self

    def get_measures(self) -> dict:
        return {
            "fqlt_start": self.start,
            "fqlt_end": self.end,
            "fqlt_duration": self.fq_duration,
            "totalq": self.totalq,
            "meanquiescent": self.meanquiescent,
            "meanquiescentout": self.meanquiescentout,
            "numberofbout": self.numberofbout,
            "qmean": self.qmean,
            "amean": self.amean,
            "qmedian": self.qmedian,
            "amedian": self.amedian,
            "qfreq": self.qfreq,
        }

    def crop_and_pad_foq(
        self,
        foq: pd.Series,
        padding: float = 1,
        mode: str = "both",
        flip: bool = False,
    ) -> pd.Series:
        modes = ("both", "head", "tail")
        if mode.lower() not in modes:
            raise TypeError(f"mode only support {modes}")

        padding_size = max(int(60 * 60 * padding / self.interval), 0)
        arr_size = self.end - self.start + 1
        if mode.lower() == "both":
            arr_size += padding_size * 2
        elif mode.lower() in ("head", "tail"):
            arr_size += padding_size

        padding_foq = np.full(arr_size, fill_value=np.nan)
        start = self.start
        end = self.end
        if mode.lower() in ("both", "head"):
            start = max(start - padding_size, 0)
        if mode.lower() in ("both", "tail"):
            end = min(self.end + padding_size, len(foq))

        crop_foq = foq.values[start:end]
        shift = 0
        if mode.lower() in ("both", "head"):
            shift = max(0, padding_size - self.start)
        padding_foq[shift : len(crop_foq) + shift] = crop_foq
        if flip:
            return pd.Series(np.fliplr(padding_foq[:, None]).flatten())
        else:
            return pd.Series(padding_foq)
