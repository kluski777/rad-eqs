import numpy as np
import matplotlib.pyplot as plt

from csv import DictReader
from datetime import datetime
from sys import path
from math import factorial

path2EQ, path2CR, path2Sun = '/home/stud/praktyki/Data_analysis/EQs/query.csv', '/home/stud/praktyki/Data_analysis/pierre_auger/scalers.csv', '/home/stud/praktyki/Data_analysis/sunspots/SIDC.csv'
path += [path2EQ, path2CR, path2Sun]
EQDateFormat = "%Y-%m-%dT%H:%M:%S.%fZ"

def sec2fraction(timestamp):
    """ conversion from timestamp to year fraction """
    secondsIn4Years = 4*365.25*24*60*60
    secondsInYear = 365*24*60*60
    rest = timestamp % secondsIn4Years
    yearFraction = (timestamp - rest)/secondsInYear+ 1970
    return yearFraction + rest/(365*24*60*60)

""" conversion from seconds to date (in string) """        
sec2date = lambda arg, formatS=EQDateFormat: datetime.fromtimestamp(int(arg)).strftime(formatS)

""" convert date (in string) to seconds """
date2sec = lambda arg, formatS=EQDateFormat: int(datetime.strptime(arg, formatS).timestamp())

""" Calculate probabilty N >= k """
Pcdf = lambda n, k: factorial(n)/(factorial(n-k)*factorial(k)*2**n) + Pcdf(n, k-1) if k <= n and k >= 0 else 0
    
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
        
def makeWin42vars(EQ, CR, binWidth, startTime, stopTime): 
    binStarts = np.arange(startTime, stopTime, binWidth)
    binStops = binStarts + binWidth
    toRetEQ, toRetCR = [], []
    
    for i,j in zip(binStarts, binStops):
        toRetCR.append(np.mean( CR.rateCorr[((CR.timeRecorded > i) & (CR.timeRecorded < j)) ] ))
        toRetEQ.append(np.sum( EQ.mag[((EQ.timeRecorded > i) & (EQ.timeRecorded < j)) ] ))
        
    toRetEQ, toRetCR = np.array(toRetEQ), np.array(toRetCR)
    indexes = ~np.isnan(toRetCR + toRetEQ)
    return toRetEQ[indexes], toRetCR[indexes]

def dayShiftProbability(startDay, stopDay, shiftStop):
    P = 1675*24*60*60
    d = 5*24*60*60 # length of a window over which the sum / mean value will be calculated
    shift = np.arange(startDay,stopDay,shiftStop)
    res = np.zeros(len(shift))
    czas = EQData.timeRecorded # co ciekawe mają inne adresy
    
    for i in range(len(res)):
        EQData.timeRecorded = czas + shift[i]*24*60*60 # shifting EQ data
        
        binStarts   = np.arange(date2sec('2014-04-02T22:07:12.00Z'), EQData.timeRecorded[-1] - P, d)
        binEnds     = binStarts + P
        print(f"{i/len(res)*100}% completed\n", end = '\r', flush=True)
        probability, times = [], []
        
        for j,l in zip(binStarts,binEnds):
            sumMag, meanCR = makeWin42vars(EQData, CRData, d, j, l)
            if np.median(sumMag) != 0 and np.median(meanCR) != 0:
                c = np.where((sumMag/np.median(sumMag)-1)*(meanCR/np.median(meanCR)-1) > 0, 1, 0)
                probability.append(Pcdf(len(c), sum(c)))
                times.append(sec2fraction(l)) # tylko po to żeby wyświetlać lata na wykresie
        
        res[i] = np.min(probability[:-75])
    
    
    plt.title("Minimum probabilities with dependency\nto the number of days EQ data was shifted")
    plt.plot(shift, np.log10(res))
    plt.xlabel("Number of days that EQ data was shifted")
    plt.ylabel(r"$log_{10}(P_{cdf})$")

SunData = SunClass(path2Sun)
EQData = EQClass(path2EQ)
CRData = CRClass(path2CR)

# dobra tu się zaczynamy korelacją pearsona zajmować, to bardzo proste będzie

def PearsonCorr(daysBack):
    meanRateCorrs = medianRateCorrs = magnitudes = []

    for i in range(len(EQData.timeRecorded)):
        rateCorrs = CRData.rateCorr[((CRData.timeRecorded < EQData.timeRecorded[i]) & (CRData.timeRecorded > EQData.timeRecorded[i] - daysBack))]
        if len(rateCorrs) != 0:
            magnitudes.append(EQData.mag[i])
            meanRateCorrs.append(np.mean(rateCorrs))
            medianRateCorrs.append(np.median(rateCorrs))
        
    print(f"Długość danych: {len(magnitudes)}")
    return np.corrcoef(meanRateCorrs, magnitudes)[0,1], np.corrcoef(medianRateCorrs, EQData.mag)[0,1]


Range = np.arange(.1,3,.1)
Range *= 24*60*60
arrLength = len(Range)
meanCorr, medianCorr = np.zeros(arrLength), np.zeros(arrLength)


# for i in range(arrLength):
#     meanCorr[i], medianCorr[i] = PearsonCorr(Range[i])

meanRateCorrs, medianRateCorrs, magnitudes = [], [], []


for i in range(len(EQData.timeRecorded)):
    rateCorrs = CRData.rateCorr[((CRData.timeRecorded < EQData.timeRecorded[i]) & (CRData.timeRecorded > EQData.timeRecorded[i] - 5*24*60*60))]
    if len(rateCorrs) != 0:
        magnitudes.append(EQData.mag[i])
        meanRateCorrs.append(np.mean(rateCorrs))
        medianRateCorrs.append(np.median(rateCorrs))




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