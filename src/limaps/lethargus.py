import logging
from dataclasses import dataclass, field
from typing import Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


LETHARGUS_DTYPES = np.dtype(
    [
        ("interval", "f8"),
        ("start", "i8"),
        ("end", "i8"),
        ("meanquiescent", "f8"),
        ("meanquiescentout", "f8"),
        ("totalq", "i8"),
        ("numberofbout", "f8"),
        ("qmean", "f8"),
        ("amean", "f8"),
        ("qmedian", "f8"),
        ("amedian", "f8"),
        ("fq_duration", "f8"),
        ("qfreq", "f8"),
    ]
)


@dataclass
class Lethargus:
    interval: int = field()
    start: int = field()
    end: int = field()
    meanquiescent: float = field(init=False, default=np.nan)
    meanquiescentout: float = field(init=False, default=np.nan)
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

    def calc_measurements(
        self,
        qaboolean: pd.Series,
        heatshock: bool = False,
    ) -> "Lethargus":
        # lethargus start and end
        ltstart = self.start
        ltend = self.end
        # duration
        fq_duration = (ltend - ltstart) * self.interval / 60 / 60

        # here cause some trouble, is arawdata was seriise totalq is not list.s
        totalq = qaboolean[ltstart:ltend].sum() * self.interval / 60
        logger.info(f"totalq {int(totalq)} min")

        # totalq all imaging duration
        totalqall = qaboolean.sum() * self.interval / 60
        # totalq out of lethargus
        totalqout = totalqall - totalq

        # quiescent / time neary= mean foq
        meanquiescent = totalq / ((ltend - ltstart) * self.interval / 60)
        logger.info(f"meanquiescent {meanquiescent}")

        # 180226 meanquiescentout is not correct way if multiple lethargus are deteced
        meanquiescentout = totalqout / (
            (len(qaboolean) - (ltend - ltstart)) * self.interval / 60
        )
        logger.info(f"meanquiescentout {meanquiescentout}")

        self.fq_duration = fq_duration
        self.totalq = totalq
        self.meanquiescent = meanquiescent
        self.meanquiescentout = meanquiescentout

        if heatshock:
            return self

        # calc q and a bout duration
        qabooleandiff = qaboolean.astype("i1").diff()
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
        if len(lastrow):
            qend_aduration = float(lastrow["qend"]) + float(lastrow["aduration"])

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
        qfreq = numberofbout / fq_duration
        logger.info(f"qfreq {np.round(qfreq, 2)} num/hour")

        self.numberofbout = numberofbout
        self.qmean = qmean
        self.amean = amean
        self.qmedian = qmedian
        self.amedian = amedian
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
        padding: Tuple[float, float] = (1.0, 1.0),
        flip: bool = False,
    ) -> pd.Series:
        if not (isinstance(padding, tuple) and len(padding) == 2):
            raise ValueError(f"padding is not a size-2 tuple")

        pre, post = padding

        pre = max(int(60 * 60 * pre / self.interval), 0)
        post = max(int(60 * 60 * post / self.interval), 0)
        arr_size = self.end - self.start + 1 + pre + post

        padding_foq = np.full(arr_size, fill_value=np.nan)
        start = self.start
        end = self.end

        start = max(start - pre, 0)
        end = min(self.end + post, len(foq))

        crop_foq = foq.values[start:end]

        shift = max(0, pre - self.start)
        padding_foq[shift : len(crop_foq) + shift] = crop_foq
        if flip:
            return pd.Series(np.fliplr(padding_foq[:, None]).flatten())
        else:
            return pd.Series(padding_foq)

    def to_numpy_arr(self) -> np.ndarray:
        data = [
            self.interval,
            self.start,
            self.end,
            self.meanquiescent,
            self.meanquiescentout,
            self.totalq,
            self.numberofbout,
            self.qmean,
            self.amean,
            self.qmedian,
            self.amedian,
            self.fq_duration,
            self.qfreq,
        ]
        data = tuple(val if not np.isnan(val) else -1 for val in data)
        return np.array(
            [data],
            dtype=LETHARGUS_DTYPES,
        )

    @classmethod
    def from_numpy_arr(
        cls,
        interval: int,
        start: int,
        end: int,
        meanquiescent: float,
        meanquiescentout: float,
        totalq: int,
        numberofbout: float,
        qmean: float,
        amean: float,
        qmedian: float,
        amedian: float,
        fq_duration: float,
        qfreq: float,
    ):
        leth = cls(interval, start, end)
        leth.meanquiescent = meanquiescent
        leth.meanquiescentout = meanquiescentout
        leth.totalq = totalq
        leth.numberofbout = numberofbout
        leth.qmean = qmean
        leth.amean = amean
        leth.qmedian = qmedian
        leth.amedian = amedian
        leth.fq_duration = fq_duration
        leth.qfreq = qfreq
        return leth
