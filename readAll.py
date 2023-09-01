import numpy as np
import matplotlib.pyplot as plt

from csv import DictReader
from datetime import datetime
from sys import path
from math import factorial

path2EQ, path2CR, path2Sun = '/home/stud/praktyki/Data_analysis/EQs/wwEQ.csv', '/home/stud/praktyki/Data_analysis/pierre_auger/scalers.csv', '/home/stud/praktyki/Data_analysis/sunspots/SIDC.csv'
path += [path2EQ, path2CR, path2Sun]
EQDateFormat = "%Y-%m-%dT%H:%M:%S.%fZ"

def sec2fraction(arg, formatS=EQDateFormat):
    """ conversion from timestamp to year fraction """
    yrFrac = 0.0
    secInYear = 365*24*60*60
    secIn4Years = 4*secInYear
    while arg > secIn4Years:
        arg -= secIn4Years
        yrFrac += 4
    while arg > secInYear:
        arg -= secInYear
        yrFrac += 1
    yrFrac += arg/secInYear
    return yrFrac+1970

""" conversion from seconds to date (in string) """        
sec2date = lambda arg, formatS=EQDateFormat: datetime.fromtimestamp(int(arg)).strftime(formatS)

""" convert date (in string) to seconds """
date2sec = lambda arg, formatS=EQDateFormat: int(datetime.strptime(arg, formatS).timestamp())

def Pcdf(n, k):
    if k <= n and k >= 0:
        return factorial(n)/(factorial(n-k)*factorial(k)*2**n) + Pcdf(n, k-1)
    else:
        return 0
    
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
        
    def median(self, step, period, startTime):
        toRet = []
        temp = startTime
        
        while temp < startTime + period:
            toRet.append(np.mean(self.mag[((self.timeRecorded < temp + step) & (self.timeRecorded > temp))]))
            temp += step
        
        return np.median(toRet)
        
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
        
    def median(self, step, period, startTime):
        toRet = []
        temp = startTime
        
        while temp < startTime + period:
            toRet.append(sum(self.rateCorr[((self.timeRecorded < temp + step) & (self.timeRecorded > temp))]))
            temp += step
        
        return np.median(toRet)
        
def retWindow(var, time, timeWidth, timeStart, timeStop, operation = 'average'):
    """ 
    var - variable over which we'll calculate means / sums,
    time - vector of time which specifies during what time certain var value was obtained,
    timeWidth - width of 1 window
    """
    toRet = []
    binTimes = []
    curTime = timeStart
    
    while curTime < timeStop:
        binTimes.append(curTime)
        if operation == 'average':
            toRet.append(np.mean(var[((time >= curTime) & (time < curTime + timeWidth))])) # jakim cude tu 0 może być?
        elif operation == 'sum':
            toRet.append(np.sum(var[(( time >= curTime) & (time < curTime + timeWidth))]))
        else:
            raise "Wrong operation argument sent to a function 'makeWindow'"
        curTime += timeWidth

    return np.array(toRet), np.array(binTimes)

def makeWin42vars(EQ, CR, binWidth, startTime, stopTime): # CR and EQ are the classes
    curTime = startTime
    toRetEQ, toRetCR = [], []
    
    while curTime < stopTime:
        toRetCR.append(np.mean(CR.rateCorr[((CR.timeRecorded > curTime) & (CR.timeRecorded < curTime + binWidth))]))
        toRetEQ.append(np.sum(EQ.mag[((EQ.timeRecorded > curTime) & (EQ.timeRecorded < curTime + binWidth))]))
        curTime += d
        
    toRetEQ = np.array(toRetEQ)
    toRetCR = np.array(toRetCR)
    indexes = ~np.isnan(toRetCR + toRetEQ)
    return toRetEQ[indexes], toRetCR[indexes]
    
SunData = SunClass(path2Sun)
EQData = EQClass(path2EQ)
CRData = CRClass(path2CR)


P = 1675*24*60*60
d = 5*24*60*60 # length of a window (I suppose)
res = np.zeros(61)

for i in range(-30,30):
    print(f"Starting with day delay equal to {i}")
    EQData.timeRecorded = EQData.timeRecorded + i*24*60*60 # shifting EQ data
    
    curTime     = date2sec('2014-04-02T22:07:12.00Z')  # CRData.timeRecorded[0]
    periodTime  = date2sec('2014-04-02T22:07:12.00Z')
    c, probability, times = [], [], []
    stopTime = CRData.timeRecorded[-1] - P
    
    
    while curTime < stopTime:
        sumMag, meanCR = makeWin42vars(EQData, CRData, d, curTime, curTime + P)
        # if np.median(sumMag) != 0 and np.median(meanCR) != 0: # zaburzać czas zaczyna lol
        c = np.where((sumMag/np.median(sumMag)-1)*(meanCR/np.median(meanCR)-1) > 0, 1, 0)
        probability.append(Pcdf(len(c), sum(c)))
        times.append(sec2fraction(curTime))
        curTime += d

    res[i+30] = np.min(probability)

plt.title("Minimum probabilities with dependency to the number of days EQ data was shifted")
plt.plot(np.linspace(-30, 30, num=len(res)), np.log10(res))
plt.xlabel("Number of days that EQ data was shifted")
plt.ylabek("Minimum probability for given shift")


# """ #Poniżej wykres 1 """
# sunBinWidth = 30*24*60*60
# magBinWidth = 5*24*60*60
# crBinWidth = 5*24*60*60 # ew. 0.9915 dnia

# spotsAvg, sunBinTimes = retWindow(SunData.spots, SunData.timeRecorded, sunBinWidth, 'average')
# magSum, magBinTimes = retWindow(EQData.mag, EQData.timeRecorded, magBinWidth, 'sum')
# rateCorrAvg,  crBinTimes = retWindow(CRData.rateCorr, CRData.timeRecorded, crBinWidth, 'average')
# sunBinTimes = [sec2fraction(i) for i in sunBinTimes]
# crBinTimes = [sec2fraction(i) for i in crBinTimes]
# magBinTimes = [sec2fraction(i) for i in magBinTimes]

# plt.plot(sunBinTimes[425:], spotsAvg[425:], color='black', label='Sun spots')
# plt.xticks(np.linspace(sunBinTimes[425],sunBinTimes[-1], num=4, dtype=int))
# plt.xlabel('year')
# plt.ylabel('Sun spots | rateCorr ')
# plt.plot(crBinTimes, rateCorrAvg, label='Pierre Auger')
# plt.legend(loc='upper left')

# ax2 = plt.gca().twinx()
# ax2.plot(magBinTimes, magSum, linewidth=.5, color='tab:orange', label='Earthquakes')
# ax2.set_ylabel('Magnitude sum', color='tab:orange')
# ax2.legend(loc='upper right')