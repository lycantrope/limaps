# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 13:14:31 2017

"""
"""
###########################################################################
lethargus and inter mold duration analysis python script

Copyright (c) 2018, Taizo kawano

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
###########################################################################
20180223
start develop based on laps.


This script process one experiment at a time and output the several measurements;
lethargus duration, total quiescent, number of bout etc.
For futher analysis and figure prepareation, 
you should use your own program, excel or whatever you like.

How to use
Before run this script, Prepare .xls file using imageJ's Multi mesure function.
Name the .xls file as "date_groupname_experimentnumber_interval.xls" format.
e.g. 170719_oy59_1_2.xls

Open Spyder and laps.py (this script).
Run this script (click triangle or push F5). 
It will show file choose dialog.
Choose the .xls file you made.
The script read the file and shows foq graphs of each sample.

In some case, it will detect multiple possible lethargus periods.
You have to choose the one that most likely to be correct by clicking the graph.
Also sometimes it mistakenly detect quiescent perid as lethargus.
You should always inspect the result by your eyes,
 and eliminate such data from futher analysis.

The script outputs following files.

1. date_groupname_experimentnumber_samplenumber_foq.png
    Each sample's graph of fractoin of q

2. date_groupname_experimentnumber_samplenumber_foq.csv
    Each sample's fractoin of q.
    The second and third rows indicate beginning and end of lethargus

3. date_groupname_experimentnumber_summary.csv
    Table of measurements.
    fqlt_start: start point of lethargus
    fqlt_end: end point of lethargus
    fqlt_duration: duration (hrs) of lethargus
    totalq: total quiescence during lethargus(min)
    meanquiescent: quiescent time/lethargus duration (totalq/fqlt_duration/60)
    meanquiescentout: quiescent time/non-lethargus time
    numberofbout: quiescent number during lethargus
    qmean: mean duration of quiescent bout (sec)
    amean: mean duration of active bout(sec)
    qmedian: median duration of quiescent bout (frame)
    amedian: median duration of active bout (frame)
    qfreq: numberofbout/fqlt_duration (number/hr)

4. date_groupname_experimentnumber_fqalll_lethargus.png
    A graph contains foq of each sample.
    
5. date_groupname_experimentnumber_fqaligned_lethargus.png
    foq aligned at beginning of lethargus
    dotted lines: each sample
    black line: mean
    dashed line: mean +- sd
    
6. date_groupname_experimentnumber_fqalignedtail_lethargus.png
    foq aligned at end of lethargus

7. date_groupname_experimentnumber__fqlt_df.csv
    foq table during lethargus +- 1hr
    
As mentioned above, the script sometimes fails lethargus detection,
so these files contains incorrect data.
You have to eliminate the incorrect data for publication quality analysis.    

20180123 ver.1.5.1
    lethargus detection criteria; longer than 2hrs.

20171211 ver.1.5
    1. foqthreshold is changeable at the top of the script.
    2. fqlt_df.csv and summary.csv contains interval data at the filename

20171122 ver.1.4
    1. foq plot of all sample 
    2. area rate of all sample 
    are saved

20170823 ver.1.3 
    1. save foq csv all data
    2. plot and save area figure 
    3. samplegroup functions for manual analysis;

    showcandidates(self, _samplenum): This method create a plot of 
    foq and onset exit candidate. also print the number of frame.
    
    manualanalysis(self, _samplenum, _start, _end): input the frame number.
    
    eg.
    : sg.showcandidates(1)
    start 764,7497
    end 235,6026,8018 
    :sg.manualanalysis(1,764,6026)
    totalq 95 min
    qend + aduration 7372.0
     ltend+inoutshift 6326
    trim the lastrow data
    numberofbout 306
    qmean 18.61sec
    amean 15.27sec
    qmedian 3.0frame
    amedian 2.0frame
    qfreq 104.68/hour
    meanquiescent 0.542189281642
    meanquiescentout 0.0067892790191
    ------------------
    man_start 764
    man_end 6026
    man_duration 2.9233333333333333
    totalq 95.1
    meanquiescent 0.542189281642
    meanquiescentout 0.0067892790191
    numberofbout 306
    qmean 18.607843137254903
    amean 15.265573770491804
    qmedian 3.0
    amedian 2.0
    qfreq 104.67502850627137
    
20170810 ver.1.2 clicking outside of lethargus period ignore the sample
20170807 ver.1.1 handle .csv as well
20170720 ver.1

"""
# %%
import os
import sys
import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from joblib import Parallel, delayed
from modules import Lethargus, Individual, LethargusDetector, Samplegroup, dotplots
from utils import (
    getTargetfile,
    getSampleIndex,
    deloutlier,
    parserFilename,
    readData,
    saveFigure,
)
from pathlib import Path

plt.rcParams["figure.max_open_warning"] = 100


# 180309 following usr defined parameters has to be set
foqthreshold = 0.2
# l3l4 24x hole use 0.05? rpi also 0.05 might better
# foqthreshold = 0.05
interval = 2
# 180309 to handle an experiments cont8ains multiple genotype..
# 6x8 config
colnum = 8
rownum = 6

# TODO
grouporder = "h"
# grouporder = "v"
uniquegroupnames = [r"aptf-1;pkc-1;dgk-1", r"k183x3F1"]
# uniquegroupnames = ["kk"]
# uniquegroupnames = [r"N2",r"rem34",r"rem34;remEx331","rem34;remEx335"]

# in this case from 1st to 4th colums are initial group etc
# gindex = [(1,4),(5,8)]
gindex = [(1, 4), (5, 6)]
# gindex = [(1),(2),(3),(4),(5),(6),(7),(8)]
# gindex = [(1,2),(3),(4),(5,6),(7,8)]
lethargus_aligned = 1
_ems = 1
if _ems:
    uniquegroupnames = ["aptf-1;pkc-1;dgk-1"]
    gindex = [(1, 8)]

assert len(gindex) == len(
    uniquegroupnames
), "[ERROR] The sizes of uniquegroupnames and gindex are not equal"

# 180309 above usr defined parameters has to be set

#################################################################
indexgrid = np.array(np.arange(colnum * rownum) + 1).reshape(rownum, colnum)
if grouporder == "h":
    indexgrid = np.fliplr(indexgrid.T)
"""
gsamplenumlist = []
for i in range(len(gindex)):
    ilist = list(range(gindex[i][0]-1,gindex[i][1]))C
    templist = []
    for j in ilist:
        templist.extend(indexgrid[j,])
    print(uniquegroupnames[i]+" "+str(templist))
    gsamplenumlist.append(templist)
"""
# 181226 mod for single column group
gsamplenumlist = getSampleIndex(gindex, uniquegroupnames, colnum, rownum, grouporder)


scriptname = "limaps"
print(scriptname)
# print(str(time.time()))
print(str(datetime.datetime.today()))
print(os.getcwd())

targetfile, targetdir, separater = getTargetfile()
thedate, groupname, expnum, interval = parserFilename(targetfile)
runmedwindow = int(60 / interval) + 1
areadf = readData(targetfile, sep=separater)

###############################################################################
def process(gn, indice, rawdata, interval):
    print(gn, indice)
    # this may not good? just use groupname?
    uniquegroupname = "_".join([thedate, gn, str(expnum)])
    sg = Samplegroup(uniquegroupname)
    sg.setDirectory(targetdir)

    sg.processagroup(thedate, gn, expnum, interval, rawdata, threshold=foqthreshold)

    # summary output
    df = sg.makedf(ltlist=sg.fullindlist)
    sg.savesummarydf()

    ltfoqmatrix = sg.makeltmatrix(ordinalnum=0)
    if ltfoqmatrix is not None:
        tempfig = sg.makealignfigure(sg.ltfoqmatrix, True, "mean")
        sg.saveafig(tempfig, "fqaligned")
        sg.saveltmatrix(sg.ltfoqmatrix, "fqlt", ordinalnum=0)

    ltfoqmatrixtail = sg.makeltmatrix(align="tail", ordinalnum=0)
    if ltfoqmatrixtail is not None:
        tempfig = sg.makealignfigure(sg.ltfoqmatrixtail, False, "mean")
        sg.saveafig(tempfig, "fqalignedtail")

    if not df.empty:
        areacropfig = sg.makeallfigureofareacropped(sg.indlist, ordinalnum=0)
        sg.saveafig(areacropfig, "areacropped")
    # sg.saveafig(areacropfig, "areacropped_del21")

    fullfoqfig = sg.makeallfigurefoq(sg.fullindlist)  #
    fullfoqfig.subplots_adjust(hspace=0.3)
    fullfoqfig.subplots_adjust(wspace=0.01)
    sg.saveafig(fullfoqfig, "foqfull")

    areafig = sg.makeallfigureofrawarea(60, sg.fullindlist, "median")
    areafig.subplots_adjust(hspace=0.1)
    areafig.subplots_adjust(wspace=0.01)
    sg.saveafig(areafig, "rawareafull")

    plt.close("all")
    return sg


def outputsummary(sg: Samplegroup):
    params = ["qaboolean", "foq"]
    print("sg.groupname: " + sg.groupname)
    for param in params:
        try:
            getattr(sg.fullindlist[0], param)
        except AttributeError as e:
            print(f"AttributeError: {e}")
            continue
        tempqadf = pd.DataFrame(
            data={
                "_".join([ind.date, ind.groupname, str(ind.samplenum)]): getattr(
                    ind, param
                )
                for ind in sg.fullindlist
            }
        )
        suffix = f"_summary_{param}.csv"
        csvpath = os.path.join(sg.targetdir, sg.groupname + suffix)
        tempqadf.to_csv(csvpath, index=False)


def foq_heatmap(
    sg: Samplegroup,
    xlim=None,
    align=False,
    cmap="binary",
    facecolor="darkkhaki",
    figsize=(5, 2),
):
    """
    PARAMS
    ::sg Samplegroup:
    ::xlim int: maximum frame number
    ::align bool: alignment the map at onset of lethargus
    ::cmap str: colormap for heatmap
    """

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(1, 1, 1)
    data = sg.makeamatrix("foq")
    data = data[:, np.where(np.sum(~np.isnan(data), axis=0))[0]]
    sample_num = data.shape[0]
    if align:
        max_start = 0
        max_end = 0
        for i, ind in enumerate(sg.fullindlist):
            if isinstance(ind.letharguslist, list):
                if len(ind.letharguslist) > 0:
                    lt = ind.letharguslist[0]
                    max_start = max(0, min(max_start, lt.start - 1800))
                    max_end = min(data.shape[1], max(max_end, lt.end + 3600))
            else:
                continue
        new_data = np.full((sample_num, max_end - max_start), fill_value=np.nan)
        for i, ind in enumerate(sg.fullindlist):
            if isinstance(ind.letharguslist, list):
                if len(ind.letharguslist) > 0:
                    lt = ind.letharguslist[0]
                    start = max(0, lt.start - 1800)
                    end = max_end
                    if (lt.start - 1800) >= 0:
                        new_data[i, : end - start] = data[i, start:end]
                    else:
                        new_data[i, (1800 - lt.start) : max_end] = data[
                            i, : (max_end - 1800 + lt.start)
                        ]

        new_data = new_data[~np.isnan(new_data).all(axis=1), :]
        sample_num = new_data.shape[0]
        foq_matrix = new_data
        outputname = f"{sg.groupname}_foq_heatmap_align.png"
    else:
        foq_matrix = data
        outputname = f"{sg.groupname}_foq_heatmap.png"

    # plot heatmap
    im = ax.imshow(foq_matrix, aspect="auto", cmap=cmap, interpolation="none")
    # x_axis
    xlim = foq_matrix.shape[1] if xlim is None else int(xlim)

    # figure configuration
    xlim_hr = xlim // 1800
    xticklabels = np.linspace(0, xlim_hr, xlim_hr + 1, dtype=np.int8)
    if align:
        xticklabels -= 1
    # x ticks
    ax.set_xlim(0, xlim)
    ax.set_xticks(np.linspace(0, xlim_hr * 1800, xlim_hr + 1))
    ax.set_xticklabels(xticklabels)
    # y ticks
    ax.set_ylim(-0.5, sample_num - 0.5)
    ax.set_yticks(np.linspace(0, sample_num - 1, sample_num // 5 + 1, dtype=np.int8))
    yticslabels = np.linspace(0, sample_num - 1, sample_num // 5 + 1, dtype=np.int8) + 1
    ax.set_yticklabels(yticslabels)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylabel("Ind.")
    # plot figure
    ax.set_title(sg.groupname.split("_")[1])
    ax.set_facecolor(facecolor)
    fig.tight_layout()
    figpath = os.path.join(sg.targetdir, outputname)
    saveFigure(figpath, fig)
    plt.show()
    return fig


# gridindex; start at top left and goes right (row)
# same order of pyplot subplot number
def gridtosgindex(gridindex):
    gnum = None
    theindex = None
    for i, alist in enumerate(gsamplenumlist):
        if gridindex in alist:
            gnum = i
            theindex = alist.index(gridindex)
            # returns group index and sample index within the group
            return gnum, theindex


def plotgridfig(
    samplegroups,
    plotdata="foq",
    window=60,
    overlayparam=None,
    meanduration=2.64,
    sdduration=0.11,
    xlim=None,
):
    gridfig = plt.figure(figsize=(16, 6))
    assert plotdata in {"foq", "area"}, f"{plotdata}: Not implemented"
    for i in range(colnum * rownum):
        # figmultifoq.clear()
        ax = gridfig.add_subplot(rownum, colnum, i + 1)
        sgn, ti = gridtosgindex(i + 1)
        _ind = samplegroups[sgn].fullindlist[ti]
        hasLethargus = (_ind.letharguslist is not None) and (
            len(_ind.letharguslist) > 0
        )
        if plotdata == "foq":
            ax.plot(_ind.foq, linestyle="-", linewidth=1, color="black")
            if hasLethargus:
                for lt in _ind.letharguslist:
                    ax.axvline(x=lt.start, color="black")
                    ax.axvline(x=lt.end, color="black", linestyle="-")
        elif plotdata == "area":
            tempraw = _ind.calcrawarea(60)
            ax.plot(
                tempraw.rolling(window=int(window / _ind.interval), center=True)
                .median()
                .values,
                linestyle="-",
                linewidth=0.5,
                color="black",
            )
        if hasLethargus:
            if overlayparam is not None:
                if overlayparam == "fq_duration":
                    lt = _ind.letharguslist[0]
                    fontcolor = "black"
                    if lt.fq_duration > meanduration + 3 * sdduration:
                        fontcolor = "red"
                    elif lt.fq_duration < meanduration - 3 * sdduration:
                        fontcolor = "blue"
                    ax.annotate(
                        "{0:4.3}".format(lt.fq_duration),
                        xy=(lt.start, 0.1),
                        fontsize=10,
                        color=fontcolor,
                    )

        ax.set_ylim(0, 1.1)
        # ax.set_xlim(0, 5*60*60/interval)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        hrtick = np.arange(len(_ind.foq) * _ind.interval / 60 / 60).astype(int)
        ax.set_xticks(hrtick * 60 * 60 / _ind.interval)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        if xlim is not None:
            ax.set_xlim(0, xlim)
            # print("xlim set")
        # ax.set_xlabel("time (h)")
        # ax.legend(groupnames)
        label = "_".join(
            [_ind.date, _ind.groupname, str(_ind.expnum), str(_ind.samplenum)]
        )

        ax.annotate(
            label,
            xy=(0.01, 0.9),
            xycoords="axes fraction",
            fontsize=8,
            horizontalalignment="left",
            verticalalignment="bottom",
        )

    ax.set_xticklabels(hrtick)
    gridfig.tight_layout()
    gridfig.subplots_adjust(hspace=0.1)
    gridfig.subplots_adjust(wspace=0.01)
    return gridfig


#############################################################################
if __name__ == "__main__":
    sys.path.append(str(Path(__file__).absolute.parent))
    jobs_num = min(os.cpu_count(), len(uniquegroupnames))
    samplegroups = []

    try:
        samplegroups = Parallel(n_jobs=jobs_num, verbose=50)(
            delayed(process)(gn, indice, areadf.iloc[:, np.array(indice) - 1], interval)
            for (gn, indice) in zip(uniquegroupnames, gsamplenumlist)
        )
    except Exception as e:
        print("[ERROR] %s" % e)

    if len(samplegroups) != 0:
        jobs_num = min(os.cpu_count(), len(samplegroups))
        Parallel(n_jobs=jobs_num)(delayed(outputsummary)(sg) for sg in samplegroups)
        Parallel(n_jobs=jobs_num)(
            delayed(foq_heatmap)(sg, align=lethargus_aligned) for sg in samplegroups
        )
        objfilename = targetfile.stem + "_samplegroups.pkl"
        print(objfilename)
        objfilepath = os.path.join(targetdir, objfilename)
        pd.to_pickle(samplegroups, objfilepath)

    # plot grid fig
    for plotdata in ["foq", "area"]:
        window = 120 if plotdata == "area" else 60
        gridfig = plotgridfig(samplegroups, plotdata="plotdata", window=window)
        samplegroups[0].targetdir
        figpath = os.path.join(samplegroups[0].targetdir, f"grid{plotdata}.png")
        saveFigure(figpath, gridfig, dpi=150)

    # for heat induced sleep and hsp over expression
    stimflag = False
    # stimflag = True
    crophrs = 6
    # crophrs = 3
    cropframe = int(crophrs * 60 * 60 / 2)
    if stimflag:
        for plotdata in ["foq", "area"]:
            window = 120 if plotdata == "area" else 60
            gridfigxh = plotgridfig(
                samplegroups, plotdata="plotdata", window=window, xlim=cropframe
            )
            samplegroups[0].targetdir
            filename = f"grid{plotdata}{crophrs}h.png"
            figpath = os.path.join(samplegroups[0].targetdir, filename)
            saveFigure(figpath, gridfigxh, dpi=100)

        foqmatlist = []
        areamatlist = []
        for sg in samplegroups:
            # sg.indlist[0].
            foqfigcrop = sg.makeallfigurefoq(sg.fullindlist, xlim=(0, cropframe))
            foqfigcrop.subplots_adjust(hspace=0.3)
            foqfigcrop.subplots_adjust(wspace=0.01)
            sg.saveafig(foqfigcrop, "foq{0}hrs".format(crophrs))

            # temmat = sg.makeamatrix(datatype = "foq", xlim =(0, cropframe))
            # 180904
            # to eliminate some failed saples.... need to impliment later
            temmat = sg.makeamatrix(datatype="foq", xlim=(0, cropframe), autoelim=True)

            tempmatfig = sg.makealignfigure(temmat, True, "mean")
            sg.saveafig(tempmatfig, "fqalignedat0_{0}hr".format(crophrs))
            foqmatlist.append(temmat)

            temmat = sg.makeamatrix(datatype="area", xlim=(0, cropframe))
            tempmatfig = sg.makealignfigure(temmat, True, "mean")
            sg.saveafig(tempmatfig, "areaalignedat0_{0}hr".format(crophrs))
            areamatlist.append(temmat)

        # at time point 5min, 10 min...
        timepoints = np.linspace(5, 60, 12, dtype=int)
        for t in timepoints:
            data0 = foqmatlist[0][:, 30 * t]
            data1 = foqmatlist[1][:, 30 * t]
            stat, pval = stats.ttest_ind(data0, data1, equal_var=False)
            print("foq {0}min pval ".format(t), str(pval))
        for t in timepoints:
            data0 = areamatlist[0][:, 30 * t]
            data0 = data0[~np.isnan(data0)]
            data1 = areamatlist[1][:, 30 * t]
            data1 = data1[~np.isnan(data1)]
            stat, pval = stats.ttest_ind(data0, data1, equal_var=False)
            print("area {0}min pval ".format(t), str(pval))

    param = "fqlt_duration"

    drdf = pd.DataFrame(columns=["groupname", param])
    for sg in samplegroups:
        if sg.summarydf is not None:
            drdf = drdf.append(sg.summarydf.loc[:, ["groupname", param]])
    if len(drdf) != 0:
        # grouplist = anadf.groupname.unique()

        dfwool, upperlimit, lowerlimit = deloutlier(drdf, "fqlt_duration")
        meantheex = dfwool.loc[:, ["groupname", param]].mean().values[0]
        stdtheex = dfwool.loc[:, ["groupname", param]].std().values[0]
        gridfig = plotgridfig(
            samplegroups, overlayparam="fq_duration", meanduration=meantheex
        )
        gridfigpath = os.path.join(targetdir, "gridfoqwithflag_meancalc.png")
        saveFigure(gridfigpath, gridfig, dpi=100)

        # import importlib
        # importlib.reload(dotplot)

        dotplotfig = dotplots(
            drdf,
            outlier="oc",
            ylim=(0, max(drdf[param])),
            binnum=25,
            size=5,
            thickness=1.5,
            figsize=(4, 8),
        )
        dotplotfig.axes[0].set_ylabel(param)

        dotplotfig.tight_layout()
        dotplotname = f"dotplot_{param}_all.png"
        dotplotfigpath = os.path.join(targetdir, dotplotname)
        saveFigure(dotplotfigpath, dotplotfig, dpi=100)

sys.exit()
# %%
# TODO
dgk_like = True

if dgk_like:
    from letharguschooser import Letharguschooser

    lc = Letharguschooser(foqthreshold)
    ld = LethargusDetector(foqthreshold)
    new_samplegroups = []
    for sg in samplegroups:
        indlist = []
        lc = Letharguschooser(foqthreshold)
        for ind in sg.fullindlist:
            if not isinstance(ind.letharguslist, list) or len(ind.letharguslist) == 0:
                datalabel = "_".join(
                    [ind.date, ind.groupname, str(ind.expnum), str(ind.samplenum)]
                )
                lc.setadata(ind)
                ld.setadata(ind)
                ret = lc.selectlethargus()
                if ret != -1:
                    ind.savefoqincsvfile(targetdir)
                    fig, ax = ld.processdata()
                    ax.annotate(
                        datalabel,
                        xy=(0.01, 0.9),
                        xycoords="axes fraction",
                        fontsize=8,
                        horizontalalignment="left",
                        verticalalignment="bottom",
                    )
                    ind.saveafig(targetdir, fig)
                    plt.close(fig)
            indlist.append(ind)
        groupname = sg.groupname
        sg.__init__(groupname)
        sg.setDirectory(targetdir)
        sg.fullindlist = indlist
        sg.indlist = indlist
        df = sg.makedf(ltlist=indlist)
        if not df.empty:
            sg.savesummarydf()

        temp = sg.makeltmatrix(ordinalnum=0)
        tempfig = sg.makealignfigure(sg.ltfoqmatrix, True, "mean")
        sg.saveafig(tempfig, "fqaligned")

        temp = sg.makeltmatrix(align="tail", ordinalnum=0)
        temp = sg.ltfoqmatrixtail
        tempfig = sg.makealignfigure(sg.ltfoqmatrixtail, False, "mean")
        sg.saveafig(tempfig, "fqalignedtail")

        if not df.empty:
            sg.saveltmatrix(sg.ltfoqmatrix, "fqlt", ordinalnum=0)

            areacropfig = sg.makeallfigureofareacropped(sg.indlist, ordinalnum=0)
            sg.saveafig(areacropfig, "areacropped")
        # sg.saveafig(areacropfig, "areacropped_del21")

        fullfoqfig = sg.makeallfigurefoq(sg.fullindlist)  #
        fullfoqfig.subplots_adjust(hspace=0.3)
        fullfoqfig.subplots_adjust(wspace=0.01)
        sg.saveafig(fullfoqfig, "foqfull")

        areafig = sg.makeallfigureofrawarea(60, sg.fullindlist, "median")
        areafig.subplots_adjust(hspace=0.1)
        areafig.subplots_adjust(wspace=0.01)
        sg.saveafig(areafig, "rawareafull")
        plt.close("all")

    # plot grid figure
    gridfig = plotgridfig(plotdata="foq")
    gridfig.savefig(targetdir + "/" + "gridfoq.png", dpi=100, transparent=True)

    gridfig = plotgridfig(plotdata="area", window=120)
    gridfig.savefig(targetdir + "/" + "gridarea.png", dpi=100, transparent=True)

    """
    ['date', 'groupname', 'expnum', 'samplenum', 'ltnum', 'fqlt_start',
     'fqlt_end', 'fqlt_duration', 'totalq', 'meanquiescent',
     'meanquiescentout', 'numberofbout', 'qmean', 'amean', 'qmedian',
   'amedian', 'qfreq']
    """
    param = "fqlt_duration"

    drdf = pd.DataFrame(columns=["groupname", param])
    for sg in samplegroups:
        if sg.summarydf is not None:
            drdf = drdf.append(sg.summarydf.loc[:, ["groupname", param]])

    # grouplist = anadf.groupname.unique()

    dfwool, upperlimit, lowerlimit = deloutlier(drdf, param)
    meantheex = dfwool.loc[:, ["groupname", param]].mean().values[0]
    stdtheex = dfwool.loc[:, ["groupname", param]].std().values[0]
    gridfig = plotgridfig(overlayparam="fq_duration", meanduration=meantheex)

    gridfig.savefig(
        targetdir + "/" + "gridfoqwithflag_meancalc.png", dpi=100, transparent=True
    )

    sys.path.append("C:\\Users\\Hayashi_Lab\\Documents\\programs\\tkmodules")
    # import importlib
    # importlib.reload(dotplot)
    import dotplot

    dotplotfig = dotplot.dotplots(
        drdf,
        outlier="oc",
        ylim=(0, max(drdf[param])),
        binnum=25,
        size=5,
        thickness=1.5,
        figsize=(4, 8),
    )
    dotplotfig.axes[0].set_ylabel(param)
    # dotplotfig.axes[0].annotate("p={0:<8.3g}".format(pval),
    # uplvalue = upperlimit.values[0]
    # dotplotfig.axes[0].annotate("upper={0:<3.2g}".format(uplvalue),
    #                 xy=(0.03, 0.92),\
    #                 xycoords='axes fraction', fontsize=10,\
    #                 horizontalalignment='left', verticalalignment='bottom')
    dotplotfig.tight_layout()
    dotplotfig.savefig(targetdir + "/" + "dotplot" + "_" + param + "_all.png", dpi=100)

    Parallel(n_jobs=jobs_num)(delayed(outputsummary)(sg) for sg in samplegroups)
    figlist = Parallel(n_jobs=jobs_num)(
        delayed(foq_heatmap)(sg, align=True, xlim=6 * 1800) for sg in samplegroups
    )
    objfilename = os.path.splitext(targetfile)[0] + "_samplegroups.pkl"
    print(os.path.split(objfilename)[1])
    objfilepath = os.path.join(targetdir, objfilename)
    pd.to_pickle(samplegroups, objfilepath)
else:
    pass


sys.exit()


# %%

# 180501 for mt screen.

# "C:/Users/Hayashi_Lab/Documents/imageprocessed/180426rem1x2random/180426_rem1hetz_1_2.csv"

basedir = "C:/Users/Hayashi_Lab/Documents/imageprocessed/180426rem1x2random"

summaryfilename = "180426_rem1hetz_1_2_summary.csv"

sys.path.append("C:\\Users\\Hayashi_Lab\\Documents\\programs\\tkmodules")
import dotplot
import fqltheatmap

filepath = basedir + "/" + summaryfilename

tempdf = pd.read_csv(filepath, sep=",")


durationdf = tempdf.loc[:, ["samplenum", "fqlt_duration"]]

medianduration = durationdf.loc[:, ["fqlt_duration"]].median()
meanduration = durationdf.loc[:, ["fqlt_duration"]].mean()
sdduration = durationdf.loc[:, ["fqlt_duration"]].std()
# quantileduration = durationdf.loc[:,["fqlt_duration"]].quantile([.25, .5, .75])
# quantileduration.loc[0.25]-1.5*(quantileduration.loc[0.75]-quantileduration.loc[0.25])
# quantileduration.loc[0.75]+1.5*(quantileduration.loc[0.75]-quantileduration.loc[0.25])


# 180313 deleate #1,0, #4,0
# sg.deletealethargus(1,0)
# sg.deletealethargus(4,0)
# sg.deletealethargus(11,0)
# sg.deletealethargus(24,0)
# sg.calcmaxduration()
# 180309 aligned foq fig


align = "tail"
# aligned at head
align = "head"

figmultifoq = plt.figure(figsize=(8, 4))
# figmultifoq.clear()
ax = figmultifoq.add_subplot(1, 1, 1)
colorlist = ["black", "magenta", "cyan", "green", "red", "blue", "orange", "olive"]

# maxlen = max([a.shape[1] for a in fqmatrixlist])
maxlen = 0
for sg, col in zip(samplegroups, colorlist):
    if sg.ltfoqmatrixtail.shape[1] > maxlen:
        maxlen = sg.ltfoqmatrixtail.shape[1]
    # np.fliplr(sg.ltfoqmatrixtail)

j = 0
for sg, col in zip(samplegroups, colorlist):
    # reverse order of group. when put control at right side of chip
    # for sg, col in zip(samplegroups[::-1], colorlist):
    if align == "tail":
        tempmatrix = sg.ltfoqmatrixtail
        emptymatrix = np.zeros(shape=(tempmatrix.shape[0], maxlen))
        emptymatrix.fill(np.nan)
        for i, ar in enumerate(tempmatrix):
            x = ar[~np.isnan(ar)]
            if len(x) != 0:
                # emptymatrix[i][-len(x):] = x
                emptymatrix[i][-len(ar) :] = ar
        tempmatrix = emptymatrix
    else:
        tempmatrix = sg.ltfoqmatrix

    tempmedian = np.nanmedian(tempmatrix, axis=0)
    ax.plot(tempmedian, linestyle="-", linewidth=1, color=col)
    print(sg.groupname)
    ax.annotate(
        sg.groupname,
        xy=(0.8, 1 - j * 0.1),
        color=col,
        xycoords="axes fraction",
        fontsize=8,
        horizontalalignment="left",
        verticalalignment="bottom",
    )
    j = j + 1


ax.set_xticklabels([])
hrtick = np.arange(len(tempmedian) * interval / 60 / 60).astype(int)

if align != "tail":
    ax.axvline(x=1800, color="black")
    ax.set_xticks(hrtick * 60 * 60 / interval)
    ax.set_xticklabels(hrtick - 1)
    ax.set_xlim(0, 5 * 60 * 60 / interval)
else:
    ax.axvline(x=maxlen - 1800, color="black")
    shift = maxlen - (hrtick * 60 * 60 / interval)[-1]
    ax.set_xticks(hrtick * 60 * 60 / interval + shift)
    ax.set_xticklabels(hrtick - max(hrtick) + 1)
    ax.set_xlim(maxlen - 5 * 60 * 60 / interval, maxlen)


ax.set_ylim(0, 1)
# ax.set_xlim(0, 11*60*60/interval)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.set_xlabel("time (h)")
# ax.legend(groupnames)

figmultifoq.tight_layout()

# figmultifoq.savefig(targetdir +"/"+"foqalignsummarygraph.png",dpi=100)
# figmultifoq.savefig(targetdir +"/"+"foqaligntailsummarygraph.png",dpi=100)
# filename = "foqalignsummarygraph_{0}_2.png".format(align)
filename = "foqalignsummarygraph_{0}.png".format(align)
figfilepath = os.path.join(targetdir, filename)
figmultifoq.savefig(figfilepath, dpi=150, transparent=True)


sys.exit()


# summarydf analysis
# dot polots and t.test stat
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd

sys.path.append("C:\\Users\\Hayashi_Lab\\Documents\\programs\\tkmodules")
import dotplot

"""
cpdf.columns
['name', 'fqlt_start', 'fqlt_end', 'fqlt_duration', 'totalq',
'meanquiescent', 'meanquiescentout', 'numberofbout', 'qmean', 'amean',
'qmedian', 'amedian', 'qfreq', 'group']
"""
# 180313 make new summarydf after the sg.deletealethargus(1,0) etc...
# df = sg.makedf(ltlist = sg.fullindlist)

param = "fqlt_duration"
uniquegroupnames

anadf = pd.DataFrame(columns=["groupname", param])
for sg, col in zip(samplegroups, colorlist):
    if sg.summarydf is not None:
        anadf = anadf.append(sg.summarydf.loc[:, ["groupname", param]])

grouplist = anadf.groupname.unique()
# reverse order
# grouplist = grouplist[::-1]


# also put effect size in the future
# stats
data0 = anadf.loc[anadf.groupname == grouplist[0], param]
data1 = anadf.loc[anadf.groupname == grouplist[1], param]

stat, pval = stats.ttest_ind(data0, data1, equal_var=False)

stats.mannwhitneyu(data0, data1, alternative="two-sided")

# tukey
# statresult = pairwise_tukeyhsd(anadf[param],
# 190822 add .astype("float") from some time ago, this cause errer only about numofbout
statresult = pairwise_tukeyhsd(
    anadf[param].astype("float"), anadf["groupname"], alpha=0.05
)
print(statresult)


# dotplot
# dotplotfig = dotplot.dotplots(anadf, groupnames = ["cont", "stim"],
dotplotfig = dotplot.dotplots(
    anadf,
    groupnames=grouplist,
    # dotplotfig = dotplot.dotplots(anadf,
    ylim=(0, max(anadf[param])),
    binnum=25,
    size=5,
    thickness=1.5,
    figsize=(4, 8),
)
dotplotfig.axes[0].set_ylabel(param)
dotplotfig.axes[0].annotate(
    "p={0:<8.3g}".format(pval),
    xy=(0.03, 0.02),
    xycoords="axes fraction",
    fontsize=10,
    horizontalalignment="left",
    verticalalignment="bottom",
)
dotplotfig.tight_layout()
dotplotfig.savefig(targetdir + "/" + "dotplot" + "_" + param + ".png", dpi=100)


####################################################################
# multi lethargus L3 L4 analysis


namelist = []
ftoslist = []
# not start of lethargus. end of lethargus?
startof0list = []
startof1list = []
for gn in uniquegroupnames:
    # multiple lethargus sampole group
    sgm = Samplegroup("_".join([thedate, gn, str(expnum)]))
    for ind in sg.fullindlist:
        numoflet = len(ind.letharguslist)
        # print("numoflet " + str(numoflet))
        if numoflet > 1:
            first_to_sec = ind.letharguslist[1].end - ind.letharguslist[0].end
            print(first_to_sec)
            ind.interlethargus.append(first_to_sec * ind.interval / 60 / 60)
            print(ind.interlethargus)
            sgm.appendanindividual(ind)
            thename = "_".join(
                [ind.date, ind.groupname, str(ind.expnum), str(ind.samplenum)]
            )
            namelist.append(thename)
            ftoslist.append(first_to_sec)
            startof0list.append(ind.letharguslist[0].end)
            startof1list.append(ind.letharguslist[1].end)


len(sgm.fullindlist)

# df[[False,multiltboolean[1:]]]
# mfig = sgm.prepallfigureformat([0,len(sgm.fullindlist[0].rawdata)],
#                                        [0,1.1],
#                                        sgm.fullindlist)
mfig = sgm.makeallfigurefoq(sgm.fullindlist)
mfig.subplots_adjust(hspace=0.3)
mfig.subplots_adjust(wspace=0.01)
sgm.saveafig(mfig, "mfoqfull")


multiltboolean = df.ltnum.values.astype(int) > 0
# multiltboolean[1:]
multidf = df[multiltboolean]
# df.loc[:,["samplenum"]] == multidf.samplenum[2]

interlethargusdf = pd.DataFrame(
    {
        "name": namelist,
        "first_to_sec": ftoslist,
        "start_0": startof0list,
        "start_1": startof1list,
    }
)


for ind in sgm.fullindlist:
    print(ind.interlethargus)


subdf = interlethargusdf[["name", "first_to_sec"]]
subdf.columns = ["group", "duration"]
dotplotfig = dotplot.dotplots(
    subdf,
    ylim=(0, max(subdf["duration"])),
    binnum=25,
    size=5,
    thickness=1.5,
    figsize=(4, 8),
)

sdfig = plt.figure(figsize=(8, 8))
ax = sdfig.add_subplot(1, 1, 1)
# ax.set_ylim(-0.1,1)
# ax.set_xlim(0, 10)
# ax.set_xlim(0, 10)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
# interval = self.indlist[0].interval
# hrtick = np.arange(_amatrix.shape[1]*interval/60/60).astype(int)
# ax.set_xticks(hrtick*60*60/interval)
# ax.set_xticklabels(hrtick-1)
ax.scatter(
    interlethargusdf.start_1 * 4 / 60 / 60,
    interlethargusdf.first_to_sec * 4 / 60 / 60,
    linestyle=":",
    linewidth=0.5,
    color="gray",
)
ax.scatter(
    interlethargusdf.start_0 * 4 / 60 / 60,
    interlethargusdf.first_to_sec * 4 / 60 / 60,
    linestyle=":",
    linewidth=0.5,
    color="gray",
)
sdfig.tight_layout()
sdfig.savefig(targetdir + "/" + "start0_l4periodscatter.png", dpi=100)
sdfig.savefig(targetdir + "/" + "start1_l4periodscatter.png", dpi=100)

sys.exit()


basedir = (
    "C:/Users/Hayashi_Lab/Documents/imageprocessed/180226n2l3l4light100/multi/compare"
)

summaryfilename = "full.csv"

sys.path.append("C:\\Users\\Hayashi_Lab\\Documents\\programs\\tkmodules")
import dotplot
import fqltheatmap


filepath = basedir + "/" + summaryfilename
tempdf = pd.read_csv(filepath, sep=",")

intermoltduration = []
groupname2 = []
for i in range(len(tempdf)):
    arow = tempdf.iloc[i]
    if arow.ltnum == 1:
        # subdf = subdf.append(tempdf.iloc[i-1])
        # subdf = subdf.append(tempdf.iloc[i])
        durationframe = tempdf.iloc[i].fqlt_end - tempdf.iloc[i - 1].fqlt_end
        durationhrs = durationframe * interval / 60 / 60
        intermoltduration.append(durationhrs)
    else:
        intermoltduration.append(0)
    if arow.date == 180222:
        groupname2.append("cont")
    elif arow.date == 180226:
        groupname2.append("stim")

imdf = pd.DataFrame({"group": groupname2, "duration": intermoltduration})

subdf = imdf[imdf.duration > 0]
subdf = subdf[["group", "duration"]]

import dotplot

dotplotfig = dotplot.dotplots(
    subdf,
    ylim=(0, max(subdf["duration"])),
    binnum=25,
    size=5,
    thickness=1.5,
    figsize=(4, 8),
)
dotplotfig.axes[0].set_ylabel("duration")
dotplotfig.tight_layout()
from scipy import stats

stat, pval = stats.ttest_ind(
    subdf[subdf.group == "cont"].duration,
    subdf[subdf.group == "stim"].duration,
    equal_var=False,
)
dotplotfig.axes[0].annotate(
    pval,
    xy=(0.02, 0.02),
    xycoords="axes fraction",
    fontsize=8,
    horizontalalignment="left",
    verticalalignment="bottom",
)

dotplotfig.savefig(targetdir + "/" + "intermoltdurationdotplot.png", dpi=100)


tempdf["intermoltduration"] = intermoltduration

tempdf.iloc[tempdf.groupname == "n2"]

tempdf.iloc[:, tempdf.intermoltduration > 0]
subdf = tempdf[tempdf.intermoltduration > 0]
contdf = tempdf[tempdf.date == 180222]
contdf[contdf.intermoltduration > 0]

tempdf.iloc[:, [tempdf.ltnum > 0]]


tempdf.index([tempdf.ltnum == 1])
list(tempdf.ltnum == 1).index(True)

# .set_index("name")
#    groupnamelist = [groupnames[i] for n in range(len(tempdf))]
#    tempdfwgname = tempdf.assign(group = groupnamelist)
#
#    summarydflist.append(tempdfwgname)


"""
#6x8 config 12 
rotetedindex = np.fliplr(np.array(np.arange(6*8)+1).reshape(6,8).T)
#4x6 config
rotetedindex = np.fliplr(np.array(np.arange(4*6)+1).reshape(6,4).T)

rotetedindex.flatten()
rotetedindex[0:2,].flatten()
rotetedindex[2:4,].flatten()
rotetedindex[4:6,].flatten()
rotetedindex[6:8,].flatten()
rotetedindex[0:4,].flatten()
rotetedindex[4:8,].flatten()

reorderlist = [sg.fullindlist[i-1] for i in rotetedindex.flatten()]
rorderfullfoqfig = sg.makeallfigurefoq(reorderlist)
sg.saveafig(rorderfullfoqfig, "foqfullreorder")


samplenumofindlist = [ind.samplenum for ind in sg.indlist]
#indexlist = []
reorderdltpluslist = []
#for a in rotetedindex[0:2,].flatten():
for a in rotetedindex.flatten():
    if a in samplenumofindlist:
        #indexlist.append(samplenumofindlist.index(a))
        indofsample = samplenumofindlist.index(a)
        reorderdltpluslist.append(sg.indlist[indofsample])

areacropfigsub = sg.makeallfigureofareacropped(reorderdltpluslist)
#areacropfigsub = sg.makeallfigureofareacropped(reorderdltpluslist[9:])

sg.saveafig(areacropfigsub, "areacroppedsub")


#maxduration
areacropsub3h = sg.makeallfigureofareacropped(reorderdltpluslist[9:],
                                            maxduration = 3*60*60/interval)
sg.saveafig(areacropsub3h, "areacroppedsub3h")



"""

"""
"""

# 180110 heat stim compare two group on a same tip
contindex = rotetedindex[
    0:4,
].flatten()
stimlindex = rotetedindex[
    4:8,
].flatten()
# 180222 fasting 24x well, so upper 12 and lower 12
contindex = list(range(1, 13, 1))
stimlindex = list(range(13, 25, 1))
stimlindex.remove(21)

contdurationlist = []
for index in contindex:
    tempind = sg.fullindlist[index - 1]
    if tempind.fq_duration is not None:
        contdurationlist.append(tempind.fq_duration)
groupnameseries = pd.Series(["cont" for a in range(len(contdurationlist))])
contdf = pd.DataFrame({"group": groupnameseries, "duration": contdurationlist})
stimdurationlist = []
for index in stimlindex:
    tempind = sg.fullindlist[index - 1]
    if tempind.fq_duration is not None:
        stimdurationlist.append(tempind.fq_duration)
groupnameseries = pd.Series(["stim" for a in range(len(stimdurationlist))])
stimdf = pd.DataFrame({"group": groupnameseries, "duration": stimdurationlist})

durationdf = pd.concat([contdf, stimdf])
durationdf = durationdf[["group", "duration"]]
sys.path.append("C:\\Users\\Hayashi_Lab\\Documents\\programs\\tkmodules")

import dotplot

dotplotfig = dotplot.dotplots(
    durationdf,
    ylim=(0, max(durationdf["duration"])),
    binnum=25,
    size=5,
    thickness=1.5,
    figsize=(4, 8),
)
dotplotfig.axes[0].set_ylabel("duration")
dotplotfig.tight_layout()
from scipy import stats

stat, pval = stats.ttest_ind(contdurationlist, stimdurationlist, equal_var=False)
dotplotfig.axes[0].annotate(
    pval,
    xy=(0.02, 0.02),
    xycoords="axes fraction",
    fontsize=8,
    horizontalalignment="left",
    verticalalignment="bottom",
)

dotplotfig.savefig(targetdir + "/" + "durationdotplot.png", dpi=100)

from scipy import stats

stat, pval = stats.ttest_ind(contdurationlist, stimdurationlist, equal_var=False)
# pvalue=0.028
stats.mannwhitneyu(contdurationlist, stimdurationlist, alternative="two-sided")
# pvalue=0.00013


alist = [0, 1, 2, 3, 4, 5, 6]
alist[2:4]
temp = range(2, 5)
alist[temp]

# 180112 area of heat shocked
rawlistcont = [ind.calcrawarea(60) for ind in sg.fullindlist[0:24]]
rawcontdf = pd.DataFrame(rawlistcont)
rawlistmt = [ind.calcrawarea(60) for ind in sg.fullindlist[24:48]]
rawmtdf = pd.DataFrame(rawlistmt)

# foq
foqlistcont = [ind.foq for ind in sg.fullindlist[0:24]]
foqcontdf = pd.DataFrame(foqlistcont)
foqlistmt = [ind.foq for ind in sg.fullindlist[24:48]]
foqmtdf = pd.DataFrame(foqlistmt)
contfoqmed = foqcontdf.median()
mtfoqmed = foqmtdf.median()
contfoqmean = foqcontdf.mean()
mtfoqmean = foqmtdf.mean()

contmed = rawcontdf.median()
mtmed = rawmtdf.median()

contmean = rawcontdf.mean()
mtmean = rawmtdf.mean()

fig = plt.figure(figsize=(4, 2))

ax = fig.add_subplot(1, 1, 1)
ax.set_ylim(0, 1.1)
ax.set_xlim(0, len(areacont))
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
# ax.set_xticklabels([])

ax.plot(contmed.rolling(window=30).median(), color="black")
ax.plot(mtmed.rolling(window=30).median(), color="red")
ax.plot(contmean.rolling(window=30).median(), color="black")
ax.plot(mtmean.rolling(window=30).median(), color="red")
ax.plot(contfoqmed.rolling(window=30).median(), color="black")
ax.plot(mtfoqmed.rolling(window=30).median(), color="red")
ax.plot(contfoqmean.rolling(window=30).median(), color="black")
ax.plot(mtfoqmean.rolling(window=30).median(), color="red")

# mintick = np.arange(25000*interval/60).astype(int)
mintick = np.arange(0, 61, 10).astype(int)
ax.set_xticks(mintick * 60 / interval)
ax.set_xticklabels([])
ax.set_xticklabels(mintick)
ax.set_xlim(0, 1 * 60 * 60 / interval + 1)

fig.tight_layout()
fig.savefig(targetdir + "/" + "areamean.png", dpi=100)
fig.savefig(targetdir + "/" + "areamed.png", dpi=100)
fig.savefig(targetdir + "/" + "foqmed.png", dpi=100)
fig.savefig(targetdir + "/" + "foqmean.png", dpi=100)

"""
"""


"""
i=0
for i in range (len(areadf.columns)):
    #i=6
    #samplenum start from 1
    ind = Individual(thedate, groupname, expnum, interval, i+1, 
                     areadf.loc[:,areadf.columns[i]])
    datalabel = "_".join([ind.date,ind.groupname, str(ind.expnum), 
                          str(ind.samplenum)])
    print("---------------------")
    print(datalabel)
    #low level process
    ind.qaboolean = ind.calcqa(ind.rawdata)
    ind.foq = ind.calcfoq(ind.rawdata)
    ind.normalizedarea = ind.normalizearea()
    ind.arearate = ind.calcarearate()
    #lethargus detectoin by foq
    #find slice where foq over/go under the threshold
    ind.detectlethargus1stscreen_fq()
    fig = ind.preparefig()
    ax = ind.plotlowlevel(fig.get_axes()[0])
    if len(ind.fq_onsetcandidates) >0:
        #screen periods that have 1h pre/post lethargus and over 1h duration
        ind.fq_oescreendmatrix = ind.detectlethargus2ndscreen(ind.fq_qbooleanvec, ind.fq_onsetcandidates, ind.fq_exitcandidates)
        if len(ind.fq_oescreendmatrix) > 0:
            ax = ind.plotcandidates(ax, ind.fq_oescreendmatrix)
            #even after above, still may exist multiple. use human eyes.
            ind.fq_finallethargusperiod = ind.detectlethargus(ind.fq_oescreendmatrix)
            if len(ind.fq_finallethargusperiod) > 0:
                ax = ind.plotlethargus(ax, ind.fq_finallethargusperiod)
                ltstart = ind.fq_finallethargusperiod[0]
                ltend = ind.fq_finallethargusperiod[1]        
                ind.calcfoqlethargusmeasures(ltstart, ltend)
                #ind.savefoqincsvfile()
                #sg.appendanindividual(ind)
    #fig.tight_layout()
    ind.savefoqincsvfile()
    sg.appendanindividual(ind)
        
    ax.annotate(datalabel,\
                xy=(0.01, 0.9),\
                 xycoords='axes fraction', fontsize=8,\
                 horizontalalignment='left', verticalalignment='bottom')

    
    ind.saveafig(fig)




df = sg.makedf()
sg.savesummarydf()

fig = sg.makeallfigure()
sg.saveafig(fig, "fqalll")

temp = sg.makeltmatrix(True)
#tempfig = sg.makealignfigure(sg.ltfoqmatrix, True, "mean")
tempfig = sg.makealignfigure(sg.ltfoqmatrix, True, "median")
#sg.saveafig(tempfig, "fqaligned")
sg.saveafig(tempfig, "fqalignedmed")

temp = sg.makeltmatrix(False)
temp = sg.ltfoqmatrixtail
#tempfig = sg.makealignfigure(sg.ltfoqmatrixtail, False, "mean")
tempfig = sg.makealignfigure(sg.ltfoqmatrixtail, False, "median")
#sg.saveafig(tempfig, "fqalignedtail")
sg.saveafig(tempfig, "fqalignedtailmed")

sg.saveltmatrix(sg.ltfoqmatrix, "fqlt")


areacropfig = sg.makeallfigureofareacropped()
sg.saveafig(areacropfig, "areacropped")

"""

"""
to eliminate some sample from the samplegroup
sg.deleteanindividual(6)
sg.deleteanindividual(21)
sg.calcmaxduration()

"""


"""

figlatency = sg.makeallfigurewithinitialq(sg.indlist)
sg.saveafig(figlatency, "withlatency")


fullfoqfig = sg.makeallfigurefoq(sg.fullindlist)#
sg.saveafig(fullfoqfig, "foqfull")

areafig = sg.makeallfigureofarea(60, sg.indlist)#60 sec
sg.saveafig(areafig, "areafull")
areafig = sg.makeallfigureofarea(60, sg.fullindlist)#60 sec
sg.saveafig(areafig, "areafullallsample")

areafig = sg.makeallfigureofrawarea(60, sg.fullindlist, "median")#60 sec
areafig = sg.makeallfigureofrawarea(240, sg.fullindlist, "median")#180 sec
sg.saveafig(areafig, "rawareafullallsample")

areafig = sg.makeallfigureofrawarea(600, sg.fullindlist, "mean")#60 sec


sg.calcinitialqlatency()
sg.makeallfigurewithinitialq(sg.indlist)

sg.saveafig(fig, "testfig")
sg.saveafig(fig, "testfig2", dpi=300)
sg.saveafig(fig, "testfig3", dpi=300)
sg.saveafig(fig, "testfig4", dpi=300)
fig.savefig(targetdir+"/test3.png")
fig.savefig(targetdir+"/test5.png", figsize = (2,8))
fig.savefig(targetdir+"/test20.png", figsize = (8,8))
fig.get_size_inches()
fig.set_size_inches(0.5,0.1)
fig.set_size_inches(5,1)
fig.set_size_inches(5,30)
fig.get_size_inches()

fig.tight_layout()
fig.savefig(targetdir+"/test6.png", dpi = 300)
fig.savefig(targetdir+"/test8.png", dpi = 100)
 
fig.savefig(targetdir+"/test.png")

fig.savefig(targetdir+"/test2_300dpi.png", dpi = 300)
#plot for green blue on 
gon = np.array([[300,750,990,1230,1470,1710,1950,2190],[600,840,1080,1320,1560,1800,2040,2280]])
bon = np.array([[600,840,1080,1320,1560,1800,2040],[750,990,1230,1470,1710,1950,2190]])

for i in range(len(gon[0])):
    rect = plt.Rectangle((gon[0,i],0), gon[1,i]- gon[0,i],1, alpha = 0.2, color="green")
    ax.add_patch(rect)
for i in range(len(bon[0])):
    rect = plt.Rectangle((bon[0,i],0), bon[1,i]- bon[0,i],1, alpha = 0.2, color="blue")
    ax.add_patch(rect)

fig = sg.prepallfigureformat([0,len(sg.fullindlist[0].rawdata)],
                                [0,1.1],sg.fullindlist)
fig = sg.prepallfigureformat([0,3600],
                                [0,1.1],sg.fullindlist)
axlist = fig.get_axes()
i=0
for ax in axlist:
    ind = sg.fullindlist[i]
    ax.plot(ind.normalizedarea.rolling(window=int(60/ind.interval), center= True).median(), color = "black", linewidth =0.5)
    for j in range(len(gon[0])):
        rect = plt.Rectangle((gon[0,j],0), gon[1,j]- gon[0,j],1, alpha = 0.2, color="green")
        ax.add_patch(rect)
    for j in range(len(bon[0])):
        rect = plt.Rectangle((bon[0,j],0), bon[1,j]- bon[0,j],1, alpha = 0.2, color="blue")
        ax.add_patch(rect)
    i = i + 1
sg.saveafig(fig, "areafullallsamplewithbluegreen")
sg.saveafig(fig, "areafullallsamplewithbluegreen2h")

ind = sg.fullindlist[0]
fig = ind.preparefig()
ax = fig.get_axes()[0]
ax = ind.plotlowlevel(fig.get_axes()[0])
ax = ind.plotnormalizedarea(ax, 60)
ax = ind.plotlethargus(ax, ind.fq_finallethargusperiod)
ax.set_ylim(0, max(ind.rawdata))
ax.plot(ind.rawdata.rolling(window = 30).median(),  color = "gray", linewidth =0.5)
ax.clear()
ax.set_ylim(0, 1.2)
ax.plot(ind.qaboolean,  color = "gray", linewidth =0.5)
ax.plot(ind.qaboolean.rolling(window = 30).mean(),  color = "red", linewidth =0.5)
#ax.plot(ind.qaboolean.rolling(window = 30).median(),  color = "blue", linewidth =0.5)
ax.plot(ind.normalizedarea.rolling(window = 3).mean(),  color = "black", linewidth =0.5)
ax.plot(ind.normalizedarea.rolling(window = 30).mean(),  color = "black", linewidth =0.5)
ax.plot(ind.normalizedarea.rolling(window = 60).mean(),  color = "blue", linewidth =0.5)
ax.plot(ind.normalizedarea.rolling(window = 120).mean(),  color = "green", linewidth =0.5)
ax.plot(ind.normalizedarea.rolling(window = 300).mean(),  color = "green", linewidth =0.5)
ax.plot(ind.normalizedarea.rolling(window = 300).median(),  color = "green", linewidth =0.5)
ax.plot(ind.normalizedarea.rolling(window = 1800).mean(),  color = "green", linewidth =0.5)
ax.plot(ind.normalizedarea.rolling(window = 1800).median(),  color = "gray", linewidth =0.5)
ax.plot(ind.normalizedarea.rolling(window = 30).min(),  color = "gray", linewidth =0.5)
ax.plot(ind.normalizedarea.rolling(window = 30).max(),  color = "gray", linewidth =0.5)
rolmax = ind.normalizedarea.rolling(window = 30).max()
rolmin = ind.normalizedarea.rolling(window = 30).min()
ax.plot(rolmax+rolmin,  color = "gray", linewidth =0.5)
ax.plot(ind.normalizedarea.rolling(window = 300).std(),  color = "gray", linewidth =0.5)
ax.plot(ind.normalizedarea.rolling(window = 30).var(),  color = "gray", linewidth =0.5)
ax.plot(ind.normalizedarea.rolling(window = 30).corr(),  color = "gray", linewidth =0.5)
ax.plot(ind.normalizedarea.rolling(window = 300).cov(),  color = "gray", linewidth =0.5)
ax.plot(ind.normalizedarea.rolling(window = 30).skew(),  color = "black", linewidth =0.5)
ax.plot(ind.normalizedarea.rolling(window = 300).skew(),  color = "gray", linewidth =0.5)
ax.plot(ind.normalizedarea.rolling(window = 30).kurt(),  color = "gray", linewidth =0.5)
ax.plot(ind.normalizedarea.rolling(window = 300).quantile(0.25),  color = "gray", linewidth =0.5)
ax.axhline(y = 0.0, color = "gray", linewidth = 0.5, linestyle="--")
ax.axhline(y = 0.5, color = "gray", linewidth = 0.5, linestyle="--")
ax.axhline(y = 0.75, color = "gray", linewidth = 0.5, linestyle="--")


ind.calcfoqlethargusmeasures()

ax.clear()

#see the worm size. compare wt mt.
ax.plot(ind.curmaxarea)
for anind in sg.indlist:
    ax.plot(anind.curmaxarea)

ax.set_ylim(0,1200)
fig.savefig(targetdir +"/"+sg.groupname+"_curmaxarea.png",figsize=(8,2),dpi=100)
 



ind.qaboolean.astype(int).diff()


temp = ind.foq>0.05
temp.ix[:,0]#too many indexers error
temp.iloc[0,:]#too many indexers error
temp.iloc[:,0]
temp.loc[0,:]
temp.astype(int).diff()


ind.detectlethargus_fq()


fig = plt.figure(figsize = (8,2))

ax = fig.add_subplot(1,1,1)
ax.set_ylim(0,1.1)
ax.set_xlim(0, len(ind.rawdata))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xticklabels([]) 
ax.plot(ind.qaboolean)
ax.plot(ind.qaboolean.diff())
ax.plot(ind.qaboolean.astype(int).diff())
ax.plot(ind.foq)
ax.plot((ind.foq>0.05).astype(int).diff())
ax.set_ylim(-1.1,1.1)


ax.plot(temp[0])
ax.plot(temp[1])
ax.plot(temp[2])
ax.plot(temp[3])
ax.plot(temp[4])
ax.plot(temp[5])

tempmean = np.nanmean(temp, axis=0)
ax.plot(tempmean, color = "black")
tempsd = np.nanstd(temp, axis=0)
ax.plot(tempmean + tempsd,linestyle=":", linewidth =0.5, color = "black")
ax.plot(tempmean - tempsd,linestyle=":", linewidth =0.5, color = "black")

"""


"""
#dotplot for initial q latency

targetfile = "qlatency.csv"
df = pd.read_csv(targetfile, sep=",")

#df.rename(columns = {"Unnamed: 0":"group"}, inplace = True)


dotplots(df, (0,max(df["initial"])))
dotplots(df.loc[:,["group","reactivate"]], (0,max(df["reactivate"])))

dfmin = df.copy()
dfmin.loc[:,["initial","reactivate"]] = dfmin.loc[:,["initial","reactivate"]]/60

dotplots(dfmin.loc[:,["group","initial"]], (0,max(dfmin["initial"])), 10)
tempfig = dotplots(dfmin.loc[:,["group","initial"]], (0,max(dfmin["initial"])),2)

dotplots(dfmin.loc[:,["group","reactivate"]], (0,max(dfmin["reactivate"])),10)


#tukey-kramer hsd
from statsmodels.stats.multicomp import pairwise_tukeyhsd
statresult = pairwise_tukeyhsd(dfmin["initial"], dfmin["group"], alpha = 0.05)
print(statresult)

statresult = pairwise_tukeyhsd(dfmin["reactivate"], dfmin["group"], alpha = 0.05)
print(statresult)


"""

sys.exit()

20171013
inleth = [1, 3, 4, 5]
outleth = [2, 6, 9, 10, 12, 14, 15, 16, 17, 18]
20171016
inleth = [1, 10, 12, 13, 16, 17]
outleth = [2, 3, 4, 15]

inlist = []
for i in inleth:
    ind = sg.fullindlist[i - 1]
    inlist.append(ind)

inrawdata = np.zeros((len(inlist), len(inlist[0].rawdata)))
i = 0
for anind in inlist:
    inrawdata[i] = anind.normalizedarea
    i = i + 1

# ind = inlist[11]

plt.plot(pd.Series(inrawdata[0]).rolling(window=30).median())
plt.plot(pd.Series(np.median(inrawdata, axis=0)).rolling(window=300).median())


outlist = []
for i in outleth:
    ind = sg.fullindlist[i - 1]
    outlist.append(ind)

outrawdata = np.zeros((len(outlist), len(outlist[0].rawdata)))
i = 0
for anind in outlist:
    outrawdata[i] = anind.normalizedarea
    i = i + 1

# ind = sg.fullindlist[11]
fig = plt.figure(figsize=(8, 4))
# fig.add_subplot(1,1,1)
fig.add_subplot(2, 1, 1)
fig.add_subplot(2, 1, 2)
axes = fig.get_axes()
window = 30

ax = axes[0]
# 4-5 hrs + 1hr
ax.set_xlim(7200 - 1800, 9000 + 1800)
ax.set_ylim(0, 1)
# ax.set_xlim(_xlim[0], _xlim[1])
# ax.set_ylim(_ylim[0],_ylim[1])
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.plot(
    pd.Series(np.mean(outrawdata, axis=0)).rolling(window=window, center=True).median(),
    color="black",
)
rect = plt.Rectangle((7200, 0), 9000 - 7200, 1, alpha=0.2, color="blue")
ax.add_patch(rect)
ax.set_xticklabels([])

ax = axes[1]
ax.set_xlim(7200 - 1800, 9000 + 1800)
ax.set_ylim(0, 1)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.plot(
    pd.Series(np.mean(inrawdata, axis=0)).rolling(window=window, center=True).median(),
    color="black",
)
rect = plt.Rectangle((7200, 0), 9000 - 7200, 1, alpha=0.2, color="blue")
ax.add_patch(rect)
# ax.set_xticklabels([])

ax.set_xlim(7200 - 1800, 9000 + 1800)

# ax.set_xlim(1800, 1800*4)
# ax.plot(ind.normalizedarea.rolling(window=window, center = True).median(),color = "black")
ax.plot(ind.foq, color="black")
rect = plt.Rectangle((7200, 0), 9000 - 7200, 1, alpha=0.2, color="gray")
ax.add_patch(rect)


ax.set_xlim(1800, 1800 * 4)
ax.plot(ind.normalizedarea.rolling(window=window, center=True).median(), color="black")
# ax.plot(ind.foq,color = "black")
rect = plt.Rectangle((1800 * 2, 0), 1800, 1, alpha=0.2, color="gray")
ax.add_patch(rect)


sg.deleteanindividual(9)
sg.deleteanindividual(20)
sg.deleteanindividual(25)
sg.deleteanindividual(31)
sg.calcmaxduration()

temp = sg.makeltmatrix(True)
sg.calcmaxduration()


# heat map
heatmapfig = plt.figure(figsize=(8, 2))
# fig.add_subplot(1,1,1)
heatmapfig.add_subplot(1, 1, 1)
axes = heatmapfig.get_axes()
ax = axes[0]
heatmapfig.subplots_adjust(bottom=0.25, left=0.01, top=1, right=0.99)
ax.cla()
ax.set_yticklabels([])
plt.gca().get_yaxis().set_ticks_position("none")

# full sg.fullindlist
emptymatrix = np.zeros(shape=(2, 3))
emptymatrix = np.zeros(shape=(len(sg.fullindlist), len(ind.rawdata)))
emptymatrix.fill(np.nan)
i = 0
for ind in sg.fullindlist:
    ind = sg.fullindlist[i]
    emptymatrix[i, 0 : len(ind.foq)] = np.array(ind.foq).flatten()
    i = i + 1

tempmatrix = emptymatrix

# lethargusonly sg.ltfoqmatrix
tempmatrix = sg.ltfoqmatrix
ltfoqmatrixoy59 = sg.ltfoqmatrix.copy()
# lethargusonly align tail sg.ltfoqmatrixtail
tempmatrix = sg.ltfoqmatrixtail


hrtick = np.arange(tempmatrix.shape[1] * interval / 60 / 60).astype(int)
ax.set_xticks(hrtick * 60 * 60 / interval)
# ax.set_xticklabels(hrtick-1)
ax.set_xticklabels(hrtick)
ax.set_xlabel("time (h)")


# plt.imshow(tempmatrix,interpolation='none', aspect="auto",cmap=plt.cm.cubehelix, vmin=0, vmax=1)
# plt.imshow(tempmatrix,interpolation='none', aspect="auto",cmap=plt.cm.gray, vmin=0, vmax=1)
plt.imshow(
    tempmatrix, interpolation="none", aspect="auto", cmap=plt.cm.Greys, vmin=0, vmax=1
)

# drop row all nan
# mask = np.all(np.isnan(tempmatrix), axis=1)
heatmapfig.tight_layout()

heatmapfig.savefig(targetdir + "/" "foqheatmapgray.png", dpi=100)
heatmapfig.savefig(targetdir + "/" "foqheatmapgreyvirtical.png", dpi=100)
heatmapfig.savefig(targetdir + "/" "foqheatmapgrey.png", dpi=100)
heatmapfig.savefig(targetdir + "/" "foqheatmapcrop.png", dpi=100)
heatmapfig.savefig(targetdir + "/" "foqheatmapcrop2.png", dpi=100)
heatmapfig.savefig(targetdir + "/" "foqheatmapcrop2tail.png", dpi=100)

# scale holizontal
scalefig = plt.figure(figsize=(8, 1))
scaleax = scalefig.add_subplot(1, 1, 1)
scaleax.set_yticklabels([])
plt.gca().get_yaxis().set_ticks_position("none")
scalearray = np.arange(0, 1.01, 0.01)
scaleax.imshow(
    np.expand_dims(scalearray, axis=0),
    interpolation="none",
    aspect="auto",
    cmap=plt.cm.Greys,
    vmin=0,
    vmax=1,
)
# scalefig.tight_layout()
scaleax.set_xticks(np.arange(0, 101, 20))
# scaleax.set_xticklabels(np.arange(0,1.2,0.2))
# scaleax.set_xticklabels(np.arange(0,1.2,0.2),rotation=0, size = 20)
scaleax.set_xticklabels(np.arange(0, 1.2, 0.2), rotation=-90, size=20)
# scaleax.set_xticklabels([0,1,2,3,4,5,6])
# scaleax.set_xticklabels(scaleax.get_xticks()/14400)
scalefig.tight_layout()


# full sg.fullindlist
maxwidth = ltfoqmatrixoy59.shape[1]

paddingmatrix = np.zeros(shape=(len(tempmatrix), maxwidth - tempmatrix.shape[1]))
paddingmatrix.fill(np.nan)
ltfoqmatrixn2 = np.concatenate((sg.ltfoqmatrix, paddingmatrix), axis=1)

heatmapfig = plt.figure(figsize=(8, 2))
# fig.add_subplot(1,1,1)
ax = heatmapfig.add_subplot(2, 1, 1)
# axes = heatmapfig.get_axes()
# ax = axes[0]
# heatmapfig.subplots_adjust( bottom = 0.25,left = 0.01,top=1, right=0.99)
ax.set_yticklabels([])
plt.gca().get_yaxis().set_ticks_position("none")
hrtick = np.arange(ltfoqmatrixn2.shape[1] * interval / 60 / 60).astype(int)
ax.set_xticks(hrtick * 60 * 60 / interval)
ax.set_xticklabels([])
# ax.set_xticklabels(hrtick-1)
plt.imshow(
    ltfoqmatrixn2,
    interpolation="none",
    aspect="auto",
    cmap=plt.cm.Greys,
    vmin=0,
    vmax=1,
)

ax = heatmapfig.add_subplot(2, 1, 2)
ax.set_yticklabels([])
plt.gca().get_yaxis().set_ticks_position("none")
hrtick = np.arange(ltfoqmatrixn2.shape[1] * interval / 60 / 60).astype(int)
ax.set_xticks(hrtick * 60 * 60 / interval)
# ax.set_xticklabels([])
ax.set_xticklabels(hrtick - 1)
plt.imshow(
    ltfoqmatrixoy59,
    interpolation="none",
    aspect="auto",
    cmap=plt.cm.Greys,
    vmin=0,
    vmax=1,
)

heatmapfig.tight_layout()
heatmapfig.savefig(targetdir + "/" "foqheatmapn2oy59crop2.png", dpi=100)


# 20171024 3 min deprive 50% duty

stimstart = 60 / interval * 60
stimend = stimstart + 60 / interval * 60 * 4
# even index on, odd off
stimindexs = np.arange(stimstart, stimend, 60 / interval * 3)

ind = sg.fullindlist[15]
# ind = sg.fullindlist[19]
fig = plt.figure(figsize=(8, 4))
fig.clear()
fig.add_subplot(1, 1, 1)
# fig.add_subplot(2,1,1)
# fig.add_subplot(2,1,2)
axes = fig.get_axes()
window = 90

ax = axes[0]
# 4-5 hrs + 1hr
# ax.set_xlim(7200-1800, 9000+1800)
ax.set_xlim(0, 60 / interval * 60 * 3)
ax.set_xlim(0, 60 / interval * 60 * 4)
ax.set_ylim(0, 1)
# ax.set_xlim(_xlim[0], _xlim[1])
# ax.set_ylim(_ylim[0],_ylim[1])
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.plot(ind.normalizedarea.rolling(window=window, center=True).median(), color="black")

for i in range(int(len(stimindexs) / 2)):
    rect = plt.Rectangle(
        (stimindexs[i * 2], 0), 60 / interval * 3, 1, alpha=0.2, color="gray"
    )
    ax.add_patch(rect)

ax.set_xticklabels([])

hrtick = np.arange(len(ind.normalizedarea) * interval / 60 / 60).astype(int)
ax.set_xticks(hrtick * 60 * 60 / interval)
# ax.set_xticklabels(hrtick-1)
ax.set_xticklabels(hrtick)
ax.set_xlabel("time (h)")

fig.tight_layout()


# 20171024 for full control. aligned, mean graph
indindex = [1, 2, 4, 6, 8]
# 20171010 for full control. aligned, mean graph
indindex = [10, 12, 13, 14, 15, 17, 18]

# 20171025 for full stimulated samples. aligned, mean graph
indindex = [3, 6, 10, 13, 16, 21, 27]
# 20171024 for full stimulated samples. aligned, mean graph
indindex = [13, 15, 20, 21, 22]
indlist = []
for i in indindex:
    indlist.append(sg.fullindlist[i - 1])

window = 90
fig = plt.figure(figsize=(8, 4))
fig.clear()
# fig.add_subplot(1,1,1)
# fig.add_subplot(2,1,1)
# fig.add_subplot(len(indlist),1,1)
axes = fig.get_axes()
ax = axes[0]
# ind = indlist[0]
# lethargus detected ind list
indlistlt = []
i = 0
for ind in indlist:
    # ax = fig.add_subplot(len(indlist),1,i+1)
    # ax.spines['right'].set_visible(False)
    # ax.spines['top'].set_visible(False)
    ind.detectlethargus1stscreen_ar()
    fig = ind.preparefig()
    ax = ind.plotlowlevel(fig.get_axes()[0])
    ind.ar_oescreendmatrix = ind.detectlethargus2ndscreen(
        ind.ar_qbooleanvec, ind.ar_onsetcandidates, ind.ar_exitcandidates
    )
    if len(ind.ar_oescreendmatrix) > 0:
        ax = ind.plotcandidates(ax, ind.ar_oescreendmatrix)
        # even after above, still may exist multiple. use human eyes.
        ind.ar_finallethargusperiod = ind.detectlethargus(ind.ar_oescreendmatrix)
        if len(ind.ar_finallethargusperiod) > 0:
            ax = ind.plotlethargus(ax, ind.ar_finallethargusperiod)
            ltstart = ind.ar_finallethargusperiod[0]
            ltend = ind.ar_finallethargusperiod[1]
            ind.calcfoqlethargusmeasures(ltstart, ltend)
            indlistlt.append(ind)
    ax.plot(
        ind.normalizedarea.rolling(window=window, center=True).median(), color="black"
    )
    for a in ind.ar_onsetcandidates:
        ax.axvline(x=a, color="red")
    for a in ind.ar_exitcandidates:
        ax.axvline(x=a, color="blue")
    i = i + 1

fig.tight_layout()

len(indlistlt)

###
# def calcmaxduration(self):
ar_maxduration = 0
for ind in indlistlt:
    if ind.ar_finallethargusperiod is not None:
        ltstart = ind.ar_finallethargusperiod[0]
        ltend = ind.ar_finallethargusperiod[1]
        if ltend - ltstart + 1 > ar_maxduration:
            ar_maxduration = ltend - ltstart + 1

# adataseries could be foq, arearate etc. alignhead True or tail False
# def makeltmatrix(self, adataseries, start, end, alignhead):
# def makeltmatrix(self, alignhead):

emptymatrix = np.zeros(
    shape=(len(indlistlt), ar_maxduration + int(60 * 60 / indlistlt[0].interval) * 2)
)
emptymatrix.fill(np.nan)
ltareamatrix = emptymatrix.copy()

j = 0
for ind in indlistlt:
    if ind.ar_finallethargusperiod is not None:
        # print(ind.ar_finallethargusperiod)
        ltstart = ind.ar_finallethargusperiod[0]
        ltend = ind.ar_finallethargusperiod[1]
        sampleend = int(ltend + 60 * 60 / ind.interval)
        if len(ind.normalizedarea) < sampleend:
            sampleend = len(ind.normalizedarea)
        rawltnormalizedarea = ind.normalizedarea[
            int(ltstart - 60 * 60 / ind.interval) : sampleend
        ]
        print(rawltnormalizedarea)
        # ltfilterd = rawltnormalizedarea.rolling(window = 90).mean()
        ltfilterd = rawltnormalizedarea.rolling(window=90).median()
        ltareamatrix[j][0 : len(ltfilterd)] = np.array(ltfilterd).flatten()
    j = j + 1

# return ltareamatrix
ltareamatrix
ltareamatrixcontrol = ltareamatrix.copy()
ltareamatrixdep = ltareamatrix.copy

tempfig = sg.makealignfigure(ltareamatrix, True, "mean")
# tempfig = sg.makealignfigure(sg.ltfoqmatrix, True, "median")
sg.saveafig(tempfig, "depareaaligned")


fig = plt.figure(figsize=(8, 4))
fig.clear()
ax = fig.add_subplot(1, 1, 1)
mean = np.nanmean(ltareamatrix, axis=0)
mean2 = np.nanmean(ltareamatrixcontrol, axis=0)
# sd = np.nanstd(_amatrix, axis = 0)
# ax.plot(mean+sd, linestyle = "--", linewidth = 1, color = "gray")
# ax.plot(mean-sd, linestyle = "--", linewidth = 1, color = "gray")
ax.plot(mean[900:-1], linestyle="--", linewidth=1, color="red")
ax.plot(mean, linestyle="--", linewidth=1, color="black")
ax.plot(mean2, linestyle="-", linewidth=1, color="black")
ax.axvline(x=1800, color="black")

ax.set_ylim(0, 1)
ax.set_xlim(0, ar_maxduration + 3600)
ax.set_xlim(0, ar_maxduration + 1800)
ax.set_xlim(0, ar_maxduration)
# ax.set_xlim(_xlim[0], _xlim[1])
# ax.set_ylim(_ylim[0],_ylim[1])
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)


ax.set_xticklabels([])

hrtick = np.arange(len(ind.normalizedarea) * interval / 60 / 60).astype(int)
ax.set_xticks(hrtick * 60 * 60 / interval)
# ax.set_xticklabels(hrtick-1)
ax.set_xticklabels(hrtick - 1)
ax.set_xlabel("time (h)")

fig.tight_layout()


# 171124 to compare light stimulated elongated lethargus
# 20171130 to summarize light stimulation effect on duration
# foq threshold 0.2
basedir = "C:/Users/Hayashi_Lab/Documents/imageprocessed/171130summaryoflightstim"
# 20180411 IB16
basedir = "C:/Users/Hayashi_Lab/Documents/imageprocessed/" "180412ceh17light/compare"
# write up experiment use analysis here
expnames = [
    "20171121_n2_1",
    "20171128_n2_1",
    "20171129_n2_1",
    "20171205_n2_1",
    "171122_n2_1",
    "171127_n2_1",
    "20171130_n2_1",
    "20171204_n2_1",
    "20171206_n2_1",
]
# group names
groupnames = [
    "cont1",
    "cont2",
    "cont3",
    "cont4",
    "stim1",
    "stim2",
    "stim3",
    "stim4",
    "stim5",
]

expnames = ["180410_ib16_1", "180402_ib16_1", "180412_ib16_1"]
# group names
groupnames = ["cont1", "stim1", "stim2"]

# 180411 add _2 for interval containig limaps format output
summaryext = "_2_summary.csv"
fqext = "_fqlt_df.csv"

sys.path.append("C:\\Users\\Hayashi_Lab\\Documents\\programs\\tkmodules")
import dotplot
import fqltheatmap


# read summary csvs
summarydflist = []
for i in range(len(expnames)):
    filepath = basedir + "/" + expnames[i] + summaryext
    # tempdf = pd.read_csv(filepath, sep = ",").set_index("name")
    tempdf = pd.read_csv(filepath, sep=",")
    groupnamelist = [groupnames[i] for n in range(len(tempdf))]
    tempdfwgname = tempdf.assign(group=groupnamelist)

    summarydflist.append(tempdfwgname)

"""
cpdf.columns
['name', 'fqlt_start', 'fqlt_end', 'fqlt_duration', 'totalq',
'meanquiescent', 'meanquiescentout', 'numberofbout', 'qmean', 'amean',
'qmedian', 'amedian', 'qfreq', 'group']
"""
# duration
param = "meanquiescent"
subdflist = []
for i in range(len(summarydflist)):
    tempdf = summarydflist[i]
    subdf = tempdf.loc[:, [param, "group"]]
    subdf.columns = [param, "group"]
    subdflist.append(subdf.loc[:, ["group", param]])

melteddf = pd.concat(subdflist)
# dotplotfig = dotplots(melteddf, (0,max(melteddf["duration"])),
#                      size = 5, binnum = 25)
# dotplotfig = dotplots(melteddf, (0,max(melteddf["duration"])),
#                       binnum = 25, size = 5, thickness = 1.5, sort = True)

dotplotfig = dotplot.dotplots(
    melteddf,
    ylim=(0, max(melteddf[param])),
    # dotplotfig = dotplot.dotplots(melteddf, ylim = (0,4),
    binnum=25,
    size=5,
    thickness=1.5,
    figsize=(4, 8),
)
dotplotfig.axes[0].set_ylabel(param)
dotplotfig.tight_layout()
dotplotfig.savefig(basedir + "/" + param + "_dotplot.png", dpi=100)

from statsmodels.stats.multicomp import pairwise_tukeyhsd

statresult = pairwise_tukeyhsd(melteddf[param], melteddf["group"], alpha=0.05)
print(statresult)
"""
Multiple Comparison of Means - Tukey HSD,FWER=0.05
============================================
group1 group2 meandiff  lower  upper  reject
--------------------------------------------
cont1  cont2   0.0891  -0.1263 0.3044 False 
cont1  stim1   0.2874    0.07  0.5048  True 
cont1  stim2   0.2194   0.002  0.4368  True 
cont2  stim1   0.1983  -0.0191 0.4157 False 
cont2  stim2   0.1303  -0.0871 0.3477 False 
stim1  stim2  -0.0681  -0.2875 0.1514 False 
--------------------------------------------
"""

# combine cont as cont
"""
melteddfcombined = melteddf.replace(["cont1","cont2","cont3","cont4",
                                     "stim1","stim2","stim3","stim4"],
                                     ["cont","cont","cont","cont",
                                      "stim","stim","stim","stim"])
"""
melteddfcombined = melteddf.replace(
    ["cont1", "cont2", "cont3", "cont4", "stim1", "stim2", "stim3", "stim4", "stim5"],
    ["cont", "cont", "cont", "cont", "stim", "stim", "stim", "stim", "stim"],
)

data0 = melteddfcombined.loc[melteddfcombined.group == "cont", param]
data1 = melteddfcombined.loc[melteddfcombined.group == "stim", param]

from scipy import stats

stat, pval = stats.ttest_ind(data0, data1, equal_var=False)
stats.mannwhitneyu(data0, data1, alternative="two-sided")

dotplotcombfig = dotplot.dotplots(
    melteddfcombined,
    ylim=(0, max(melteddf[param])),
    size=5,
    binnum=25,
    thickness=1.5,
    figsize=(4, 8),
)
dotplotcombfig = dotplot.dotplots(
    melteddfcombined, ylim=(1, 4), size=5, binnum=25, thickness=1.5, figsize=(4, 8)
)

# sub = list(_thedata[_thedata.columns[1]][_thedata[_thedata.columns[0]]==x])
# datalist.append(sub)

# dotplotcombfig = dotplots(melteddfcombined,figsize = (4,8))
# dotplotcombfig = dotplots(melteddfcombined, (0,max(melteddf["duration"])),
#                      size = 5, binnum = 50, thickness = 1.5, sort = True)
dotplotcombfig.axes[0].set_ylabel(param)
dotplotcombfig.tight_layout()

# contvals = melteddfcombined[melteddfcombined.group == "cont"].duration.values
# stimvals = melteddfcombined[melteddfcombined.group == "stim"].duration.values


dotplotcombfig.axes[0].annotate(
    pval,
    xy=(0.02, 0.02),
    xycoords="axes fraction",
    fontsize=8,
    horizontalalignment="left",
    verticalalignment="bottom",
)

dotplotcombfig.savefig(basedir + "/" + param + "_combinedotplot.png", dpi=100)

# stats.ttest_ind(contvals,stimvals,equal_var = False)
# pvalue=0.000795
# stats.mannwhitneyu(contvals,stimvals,alternative = "two-sided")
# pvalue=6.4e-08


# anova
stats.f_oneway(*[a.duration.values for a in subdflist])
# pvalue=0.003

c1 = melteddf[melteddf.group == "cont1"].duration.values
c2 = melteddf[melteddf.group == "cont2"].duration.values
s1 = melteddf[melteddf.group == "stim1"].duration.values
s2 = melteddf[melteddf.group == "stim2"].duration.values

# welch two-side
stats.ttest_ind(c1, s1, equal_var=False)
# pvalue=0.01
stats.ttest_ind(c1, s2, equal_var=False)
# pvalue=1.4e-05
stats.ttest_ind(c2, s1, equal_var=False)
# pvalue=0.083
stats.ttest_ind(c2, s2, equal_var=False)
# pvalue=0.032

stats.mannwhitneyu(c1, s1, alternative="two-sided")
# pvalue=0.00036
stats.mannwhitneyu(c1, s2, alternative="two-sided")
# pvalue=1.2e-07
stats.mannwhitneyu(c2, s1, alternative="two-sided")
# pvalue=0.02
stats.mannwhitneyu(c2, s2, alternative="two-sided")
# pvalue=5.5e-05


# see other parameters
testparmeters = ["totalq", "meanquiescent", "numberofbout", "qmean", "amean", "qfreq"]

for j in range(len(testparmeters)):
    subdflist = []
    for i in range(len(summarydflist)):
        tempdf = summarydflist[i]
        subdf = tempdf.loc[:, [testparmeters[j], "group"]]
        subdf.columns = [testparmeters[j], "group"]
        subdflist.append(subdf.loc[:, ["group", testparmeters[j]]])

    melteddf = pd.concat(subdflist)
    dotplotfig = dotplots(melteddf, (0, max(melteddf[testparmeters[j]])), size=5)
    dotplotfig.axes[0].set_ylabel(testparmeters[j])
    dotplotfig.tight_layout()

    statresult = pairwise_tukeyhsd(
        melteddf[testparmeters[j]], melteddf["group"], alpha=0.05
    )
    print(testparmeters[j])
    print(statresult)

# all false
# eg. meanquiescent
"""
Multiple Comparison of Means - Tukey HSD,FWER=0.05
=============================================
group1 group2 meandiff  lower   upper  reject
---------------------------------------------
cont1  cont2  -0.0342  -0.0669 -0.0015  True 
cont1  stim1  -0.0413  -0.0744 -0.0083  True 
cont1  stim2  -0.0044  -0.0374  0.0287 False 
cont2  stim1  -0.0071  -0.0402  0.0259 False 
cont2  stim2   0.0298  -0.0032  0.0629 False 
stim1  stim2   0.0369   0.0036  0.0703  True 
---------------------------------------------
with in cont or stim have differ. not good parameter
all above are like this
"""


# median fracton of q graph
# fqext ="_fqlt_df.csv"
# read fraction of q csvs
fqmatrixlist = []
for i in range(len(expnames)):
    filepath = basedir + "/" + expnames[i] + fqext
    tempdf = pd.read_csv(filepath, sep=",")
    tempmat = np.array(tempdf).T
    fqmatrixlist.append(tempmat)
# group names
"""
groupnames = ["cont1",
              "cont2",
              "cont3",
              "stim1",
              "stim2"]
"""
"""
colorlist = ["black",
             "gray",
             "red",
             "orange"
             ]
"""
contnum = 4
# controls are gray
colorlist = cm.Greys(1 - 1 / (contnum + 1) * np.arange(contnum))
stimnum = 5
# colorlist = np.append(colorlist, cm.Reds(1-1/stimnum*np.arange(stimnum)),axis = 0)
# colorlist = np.append(colorlist, cm.winter(1-1/(stimnum+1)*np.arange(stimnum)),axis = 0)
colorlist = np.append(
    colorlist, cm.cool(1 - 1 / (stimnum + 1) * np.arange(stimnum)), axis=0
)

# aligned at head
figmultifoq = plt.figure(figsize=(8, 4))
# figmultifoq.clear()
ax = figmultifoq.add_subplot(1, 1, 1)

for i in range(len(fqmatrixlist)):
    tempmedian = np.nanmedian(fqmatrixlist[i], axis=0)
    ax.plot(tempmedian, linestyle="-", linewidth=1, color=colorlist[i])

ax.axvline(x=1800, color="black")

ax.set_ylim(0, 1)
ax.set_xlim(0, 5 * 60 * 60 / interval)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)

ax.set_xticklabels([])

hrtick = np.arange(len(tempmedian) * interval / 60 / 60).astype(int)
ax.set_xticks(hrtick * 60 * 60 / interval)
# ax.set_xticklabels(hrtick-1)
ax.set_xticklabels(hrtick - 1)
ax.set_xlabel("time (h)")
ax.set_xlim(0, 5 * 60 * 60 / interval)
ax.legend(groupnames)
figmultifoq.tight_layout()

figmultifoq.savefig(basedir + "/" + "foqcontandlightstim.png", dpi=100)


# aligned at tale
figmultifoqtail = plt.figure(figsize=(8, 4))
# figmultifoq.clear()
ax = figmultifoqtail.add_subplot(1, 1, 1)
maxlen = max([a.shape[1] for a in fqmatrixlist])
for i in range(len(fqmatrixlist)):
    temparray = fqmatrixlist[i]
    emptymatrix = np.zeros(shape=(temparray.shape[0], maxlen))
    emptymatrix.fill(np.nan)
    ltfoqmatrix = emptymatrix.copy()
    flipedarray = np.fliplr(temparray)
    j = 0
    for anarray in flipedarray:
        x = anarray[~np.isnan(anarray)]
        # ltstart = ind.fq_finallethargusperiod[0]
        # ltend = ind.fq_finallethargusperiod[1]
        # ltfoq = ind.foq[int(ltstart-60*60/ind.interval):int(ltend+60*60/ind.interval)]
        # if not alignhead:
        # ltfoq = ltfoq.reverse()
        # ltfoq = ltfoq.iloc[::-1]
        ltfoqmatrix[j][0 : len(x)] = x
        # ltfoqmatrix[j][0:len(ltfoq)] = np.array(ltfoq).flatten()
        j = j + 1

    tempmedian = np.nanmedian(np.fliplr(ltfoqmatrix), axis=0)
    ax.plot(tempmedian, linestyle="-", linewidth=1, color=colorlist[i])


ax.axvline(x=maxlen - 1800, color="black")

ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)

hrtick = np.arange(len(tempmedian) * interval / 60 / 60).astype(int)
ax.set_xticklabels([])
shift = maxlen - (hrtick * 60 * 60 / interval)[-1]
ax.set_xticks(hrtick * 60 * 60 / interval + shift)
ax.set_xticklabels(hrtick - max(hrtick) + 1)

ax.set_xlabel("time (h)")

ax.set_ylim(0, 1)
ax.set_xlim(maxlen - 5 * 60 * 60 / interval, maxlen)
ax.legend(groupnames)

figmultifoqtail.tight_layout()
figmultifoqtail.savefig(basedir + "/" + "alignedtailfoqgraph.png", dpi=100)


# int(hrs)*60*60/interval
maxshow = int(5) * 60 * 60 / interval
_kwargs = {
    "figsize": (8, 8),
    "groupnames": groupnames,
    "interval": interval,
    "maxshow": maxshow,
}
# heatmap of foq aligned at head
heatmapfig = fqltheatmap.heatmap(fqmatrixlist, **_kwargs)
# heatmapfig.tight_layout()
# heatmaptailfig.subplots_adjust(hspace = .05)
heatmapfig.savefig(basedir + "/" + "heatmap.png", dpi=100)


_kwargs.update({"head": False})
# heatmap of foq aligned at tail
heatmaptailfig = fqltheatmap.heatmap(fqmatrixlist, **_kwargs)
heatmaptailfig.savefig(basedir + "/" + "heatmaptail.png", dpi=100)


scalefig = fqltheatmap.makescale(figsize=(8, 1.2))
scalefig.tight_layout()
scalefig.savefig(basedir + "/" + "scale.png", dpi=100)


# see if spent time in the chamber affect
# scatter plot and linear regression
startdurationfig = plt.figure(figsize=(8, 8))
ax = startdurationfig.add_subplot(1, 1, 1)
# duration
for i in range(len(summarydflist)):
    tempdf = summarydflist[i]
    x = tempdf.fqlt_start
    y = tempdf.fqlt_duration
    ax.plot(x, y, "o", color=colorlist[i], markersize=2)

    a, b = np.polyfit(x, y, 1)
    y2 = a * x + b
    ax.plot(x, y2, color=colorlist[i])

ax.set_xlim(0, 25000)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
hrtick = np.arange(25000 * interval / 60 / 60).astype(int)
ax.set_xticks(hrtick * 60 * 60 / interval)
ax.set_xticklabels([])
ax.set_xticklabels(hrtick)
startdurationfig.tight_layout()

startdurationfig.savefig(basedir + "/" + "startdurationscatter.png", dpi=100)

ax.clear()

stats.linregress(summarydf2.fqlt_start, summarydf2.fqlt_duration)
stats.linregress(subsumdf1.fqlt_start, subsumdf1.fqlt_duration)
stats.linregress(subsumdf2.fqlt_start, subsumdf2.fqlt_duration)
stats.linregress(subsumdfs2.fqlt_start, subsumdfs2.fqlt_duration)
