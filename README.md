fitclock
===============

The `fitclock` package implements related algorithms to fit
reaction times and rewards from the clock task developed by
Michael Frank and colleagues. 

Installation
------------

To install the package from github, run the following in R:
 
```r
if (!require(devtools)) { install.packages("devtools") }

install_github("MplusAutomation", "michaelhallquist")
```

Examples
--------

The package includes an example dataset of a single subject who completed
eight runs of the clock task (50 trials per run). This dataset is called
`clocksubject_fMRI_008jh`, and you can see the basic structure:

```r
str(clocksubject_fMRI_008jh)
head(clocksubject_fMRI_008jh)
```

The fitclock package can fit a vector of parameters for a group (optimize
parameters of multiple subjects with multiple runs), a subject (many runs),
or a single run. The data storage structure for the package defines three object:
`clockdata_group`, `clockdata_subject`, and `clockdata_run`. Not surprisingly,
`clockdata_group` objects are composed of multiple `clockdata_subject` objects,
which are in turn composed of multiple `clockdata_run` objects.

The basic time-clock algorithm is designed to fit a series of trials with a
consistent contingency (i.e., a run). Thus, observed reaction times and rewards
are stored at the `clockdata_run` level, whereas the `clockdata_subject` and
`clockdata_group` objects are more storage lists for multiple runs or subjects,
respectively.

At this time, the `clockdata_group` functionality is not fully implemented,
but scaffolding is in place. Fitting subjects or runs, however, is complete,
and examples for each are provided below.

Importing data into fitclock
----------------------------

To import behavioral data into the fitclock package, write the behavioral data
for all runs as a .csv file, perhaps using the ClockToCSV.m script. This csv file should
have the following fields