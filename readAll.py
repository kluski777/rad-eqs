import numpy as np
import matplotlib.pyplot as plt

from csv import DictReader
from datetime import datetime
from sys import path

path2EQ, path2CR, path2Sun = '/home/stud/praktyki/Data_analysis/EQs/query.csv', '/home/stud/praktyki/Data_analysis/pierre_auger/scalers.csv', '/home/stud/praktyki/Data_analysis/sunspots/SIDC.csv'
path += [path2EQ, path2CR, path2Sun]
EQDateFormat = "%Y-%m-%dT%H:%M:%S.%fZ"

def sec2date(arg, formatS=EQDateFormat):
    """ conversion from seconds to date (in string) """
    return datetime.fromtimestamp(int(arg)).strftime(formatS)

def date2sec(arg, formatS=EQDateFormat):
    """ convert date (in string) to seconds """
    return int(datetime.strptime(arg, formatS).timestamp())


class SunClass:
    timeRecorded = []
    spots = []
    
    def __init__(self, path, startDate = '1970-01-01T00:00:00.00Z', stopDate = '2023-08-30T15:17:00.00Z'):
        data = DictReader(open(path))
        startDate = date2sec(startDate)
        stopDate = date2sec(stopDate)
        
        for row in data:
            if int(row['Year']) >= 1970:
                timestamp = date2sec(row['Year']+'-'+row['Month']+'-'+row['Day']+'T12:00:00.00Z')
                if timestamp >= startDate and timestamp < stopDate and row['sunspots'] != '-1':
                    self.timeRecorded.append(timestamp)
                    self.spots.append(int(row['sunspots']))
        
        self.timeRecorded = np.array(self.timeRecorded)
        self.spots = np.array(self.spots)

class EQClass:
    timeRecorded = []
    mag = []
    magType = []
    
    def __init__(self, path, startDate = '1970-01-01T00:00:00.00Z', stopDate = '2023-08-30T15:17:00.00Z'):
        data = DictReader(open(path))
        startDate = date2sec(startDate)
        stopDate = date2sec(stopDate)
        
        for row in data:
            timestamp = date2sec(row['time'])
            if timestamp >= startDate and timestamp < stopDate:
                self.timeRecorded.append(timestamp)
                self.mag.append(float(row['mag']))
                self.magType.append(row['magType'])
        
        self.timeRecorded = np.array(self.timeRecorded[::-1])
        self.mag = np.array(self.mag[::-1])
        self.magType = np.array(self.magType[::-1])
        
class CRClass:
    timeRecorded = []
    rateCorr = []
    
    def __init__(self, path, startDate = '1970-01-01T00:00:00.00Z', stopDate = '2023-08-30T15:17:00.00Z'):
        data = DictReader(open(path))
        startDate = date2sec(startDate)
        stopDate = date2sec(stopDate)
        
        for row in data:
            timestamp = int(row['time'])
            if timestamp >= startDate and timestamp < stopDate:
                self.timeRecorded.append(timestamp)
                self.rateCorr.append(float(row['rateCorr']))
        
        self.timeRecorded = np.array(self.timeRecorded)
        self.rateCorr = np.array(self.rateCorr)
        
def makeWindow(var, time, timeWidth, operation = 'average'):
    toRet = []
    binTimes = []
    curTime = time[0]
    
    while curTime < time[-1]:
        binTimes.append(curTime)
        if operation == 'average':
            toRet.append(np.mean(var[((time >= curTime) & (time < curTime + timeWidth))]))
        elif operation == 'sum':
            toRet.append(np.sum(var[(( time >= curTime) & (time < curTime + timeWidth))]))
        else:
            raise "Wrong operation argument sent to a function 'makeWindow'"
        curTime += timeWidth

    return np.array(toRet), np.array(binTimes)
    
SunData = SunClass(path2Sun)
EQData = EQClass(path2EQ)
CRData = CRClass(path2CR)

sunBinWidth = 30*24*60*60
magBinWidth = 5*24*60*60
crBinWidth = 5*24*60*60 # ew. 0.9915 dnia

spotsAvg, sunBinTimes = makeWindow(SunData.spots, SunData.timeRecorded, sunBinWidth)
magSum, magBinTimes = makeWindow(EQData.mag, EQData.timeRecorded, magBinWidth, 'sum')
rateCorrAvg,  crBinTimes = makeWindow(CRData.rateCorr, CRData.timeRecorded, crBinWidth, 'average')


plt.plot(sunBinTimes[425:], spotsAvg[425:], color='black', label='Sun spots')
plt.plot(crBinTimes, rateCorrAvg, label='Pierre Auger')
plt.legend(loc='upper left')

ax2 = plt.gca().twinx()
ax2.plot(magBinTimes,magSum, linewidth=.5, color='tab:orange', label='Earthquakes')
ax2.set_ylabel('Magnitude sum', color='tab:orange')
ax2.legend()
