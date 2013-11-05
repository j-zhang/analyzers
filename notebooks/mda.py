# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# **Imports the necessary packages**

# <codecell>

import datetime
import time
import urllib
import numpy as np
from pandas import read_csv
from dateutil import parser
from pprint import pprint

# <markdowncell>

# **Gets some sample or fake data, perhaps from hw2's Alaska?**

# <codecell>

url = 'http://earthquake.usgs.gov/earthquakes/catalogs/eqs7day-M1.txt'
data = read_csv(urllib.urlopen(url))
clean_data = data.dropna(axis=0, how='any')
test_data = clean_data[0:10]
test_data

# <markdowncell>

# **Get the test magnitude and test datetime from the USGS code**

# <codecell>

test_mag = test_data['Magnitude']
test_mag = test_mag.tolist()
test_mag

# <codecell>

test_dt = test_data['Datetime']
test_dt = test_dt.tolist()
test_dt = [parser.parse(dt) for dt in test_dt]
test_dt

# <markdowncell>

# **Extracts the alarm length for each possible quake**

# <codecell>

def mda(mag, dt):
    """
    Uses basic MDA model (tau*u^mag) to predict earthquakes.
    Returns tuple of (start, end) representing date range when alarm should be on.
    MAG: list of earthquake magnitudes
    DT: list of earthquake datetimes in python datetime format
    (MAG and DT have the same length and come from earthquakes data frame)
    """
    assert len(mag) == len(dt), "Dude are you mad?"
    
    tau = 0.7 # we will figure out what tau is later
    u = 4 # we will add fancy tuning functionality later
    alarms = []
    
    for i in range(0, len(mag)):
        alarm_length = tau * (u ** mag[i])
        val = datetime.timedelta(seconds=alarm_length)
        rng = (dt[i], dt[i]+val)
        alarms.append(rng)
    return alarms

# <markdowncell>

# **A sample run over the test magnitude and test datetime**

# <codecell>

alarm_ranges = mda(test_mag, test_dt)
alarm_ranges

# <codecell>

for alarm in alarm_ranges:
    print alarm[0], "\t", alarm[1]

# <codecell>


