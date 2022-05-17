"""
lethargus and inter mold duration analysis python script

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



20220516 ver 3.0.0 handle the heatshock quiescence, update author
"""

__version__ = "3.0.0"


from .dotplot import dotplots
from .individual import Individual
from .lethargus import Lethargus
from .lethargusdetector import LethargusDetector
from .project import Project
from .samplegroup import Samplegroup

__all__ = [
    "dotplots",
    "Individual",
    "Lethargus",
    "LethargusDetector",
    "Samplegroup",
    "Project",
]
