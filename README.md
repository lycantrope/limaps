# limaps
A python package for analyzing nematode lethargus and inter-molting behavior

# Installation

After activate a python virtual environment, run:

```console
(.venv) % pip install -U git+https://github.com/lycantrope/limaps

```

# Usage
## Import module
```python
from limaps import Project
```
---


## Create a project and process the data
```python
proj = (
    Project(
        foqthreshold=0.2,
        colnum=8,
        rownum=6,
        grouporder="h",
    )
    .set_targetfile(targetfile="./data/220426_remi_1_2.csv")
    .set_groupnames(
        uniquegroupnames=["N2"],
        gindex=[(1, 6)],
    )
    .read_dataframe()
    .process_samplegroups()
)
```

    2022-04-29 14:57:41.003 project.py:116 INFO 
    Basic Parameters:
    Threshold: 0.2
    number of column: 8
    number of row: 6
    group orientation: h
    
    2022-04-29 14:57:41.004 project.py:130 INFO 
    File information:
    filename: 220426_remi_1_2.csv
    date: 220426
    experiment name: remi
    number of experiment: 1
    interval between frame (sec): 2.0
    
    2022-04-29 14:57:41.005 project.py:191 WARNING [6, 7] are assigned as `unknown` group
    2022-04-29 14:57:41.005 project.py:202 INFO 
     Remi grid
    2022-04-29 14:57:41.006 project.py:203 INFO 
    N2_6       N2_12      N2_18      N2_24      N2_30      N2_36      unknown_6  unknown_12 
    N2_5       N2_11      N2_17      N2_23      N2_29      N2_35      unknown_5  unknown_11 
    N2_4       N2_10      N2_16      N2_22      N2_28      N2_34      unknown_4  unknown_10 
    N2_3       N2_9       N2_15      N2_21      N2_27      N2_33      unknown_3  unknown_9  
    N2_2       N2_8       N2_14      N2_20      N2_26      N2_32      unknown_2  unknown_8  
    N2_1       N2_7       N2_13      N2_19      N2_25      N2_31      unknown_1  unknown_7  
    2022-04-29 14:57:41.106 samplegroup.py:116 INFO ---------------------
    2022-04-29 14:57:41.107 samplegroup.py:117 INFO 220426_N2_1_1
    ...


---


## Ploting grid figure and dotplot

```python
proj = (
    proj.saveafig(
        figure=proj.plot_samplegroups_grid("foq", 60, "fq_duration"),
        filename="gridfoq.png",
        dpi=150,
    )
    .saveafig(
        figure=proj.plot_samplegroups_grid("area", 120, "fq_duration"),
        filename="gridarea.png",
        dpi=150,
    )
    .saveafig(
        figure=proj.dotplots("fqlt_duration"),
        filename="fq_duration_dotplot.png",
        dpi=150,
    )
)
```
    
...

---


## Save result as slide and pickle file

```python
proj  = (proj
        .create_summary_slide() # generate pptx
        .to_pickle() # save Project as pkl.gz
    )

```

    2022-04-29 14:58:16.528 project.py:398 INFO Save project at: data/220426_remi_1_2.pkl.gz

---

## Load a project from pickle gzip file

```python
proj = Project.from_pickle("./data/220426_remi_1_2.pkl.gz")
proj
```




    Project(foqthreshold=0.2, colnum=8, rownum=6, grouporder='h', uniquegroupnames=['N2'], date='220426', expname='remi', expnum='1', interval=2.0, datapath=PosixPath('data/220426_remi_1_2.csv'), homepath=PosixPath('data'))



---
## Display the instance within a Project
```python
proj.samplegroups
```




    [Samplegroup(date='220426', groupname='N2', expnum='1', threshold=0.2, homepath=PosixPath('data')),
     Samplegroup(date='220426', groupname='unknown', expnum='1', threshold=0.2, homepath=PosixPath('data'))]


---

```python
proj.samplegroups[0].fullindlist
```




    [Individual(date='220426', groupname='N2', expnum='1', interval=2.0, samplenum=1),
     Individual(date='220426', groupname='N2', expnum='1', interval=2.0, samplenum=2),
     ...
     ,
     Individual(date='220426', groupname='N2', expnum='1', interval=2.0, samplenum=36)]



---
```python
proj.samplegroups[0].fullindlist[0].letharguslist
```




    [Lethargus(interval=2.0, start=4796, end=9513, meanquiescent=0.5829976680093281, totalq=91.66666666666667, numberofbout=530, qmean=10.362264150943396, amean=6.935849056603773, qmedian=3.0, amedian=2.0, fq_duration=2.6205555555555553, qfreq=202.24719101123597)]


