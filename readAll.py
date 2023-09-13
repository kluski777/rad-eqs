import numpy as np
import matplotlib.pyplot as plt


from bisect import bisect_left, bisect_right # binary search, żeby szybciej znajdować indeksy
from csv import DictReader
from datetime import datetime
from sys import path
from math import factorial
from scipy.stats import spearmanr

# potrzebne ścieżki poniżej
path2EQ, path2CR, path2Sun = '/home/stud/praktyki/Data_analysis/EQs/', '/home/stud/praktyki/Data_analysis/pierre_auger/scalers.csv', '/home/stud/praktyki/Data_analysis/sunspots/SIDC.csv'
path2EQ = [ path2EQ+'1960_1969.csv', 
            path2EQ+'1970_1970.csv', 
            path2EQ+'1971_1971.csv', 
            path2EQ+'1972_1972.csv', 
            path2EQ+'1973_1973.csv', 
            path2EQ+'1974_1974.csv', 
            path2EQ+'1975_1975.csv', 
            path2EQ+'1976_1976.csv', 
            path2EQ+'1977_1977.csv', 
            path2EQ+'1978_1978.csv', 
            path2EQ+'1979_1979.csv', 
            path2EQ+'1980_1980.csv', 
            path2EQ+'1981_1981.csv', 
            path2EQ+'1982_1982.csv', 
            path2EQ+'1983_1983.csv', 
            path2EQ+'1984_1984.csv', 
            path2EQ+'1985_1985.csv',
            
            path2EQ+'1986_1986_0.csv',
            path2EQ+'1986_1986_1.csv',
            path2EQ+'1987_1987_0.csv',
            path2EQ+'1987_1987_1.csv',
            path2EQ+'1988_1988_0.csv',
            path2EQ+'1988_1988_1.csv',
            path2EQ+'1989_1989_0.csv',
            path2EQ+'1989_1989_1.csv',
            path2EQ+'1990_1990_0.csv',
            path2EQ+'1990_1990_1.csv',
            path2EQ+'1991_1991_0.csv',
            path2EQ+'1991_1991_1.csv',
            
            path2EQ+'1992_1992_0.csv',
            path2EQ+'1992_1992_1.csv',
            path2EQ+'1992_1992_2.csv',
            path2EQ+'1993_1993_0.csv',
            path2EQ+'1993_1993_1.csv',
            path2EQ+'1993_1993_2.csv',
            path2EQ+'1994_1994_0.csv',
            path2EQ+'1994_1994_1.csv',
            path2EQ+'1994_1994_2.csv',
            path2EQ+'1995_1995_0.csv',
            path2EQ+'1995_1995_1.csv',
            path2EQ+'1995_1995_2.csv',
            path2EQ+'1996_1996_0.csv',
            path2EQ+'1996_1996_1.csv',
            path2EQ+'1996_1996_2.csv',
            path2EQ+'1997_1997_0.csv',
            path2EQ+'1997_1997_1.csv',
            path2EQ+'1997_1997_2.csv',
            path2EQ+'1998_1998_0.csv',
            path2EQ+'1998_1998_1.csv',
            path2EQ+'1998_1998_2.csv',
            path2EQ+'1999_1999_0.csv',
            path2EQ+'1999_1999_1.csv',
            path2EQ+'1999_1999_2.csv',
            path2EQ+'2000_2000_0.csv',
            path2EQ+'2000_2000_1.csv',
            path2EQ+'2000_2000_2.csv',
            path2EQ+'2001_2001_0.csv',
            path2EQ+'2001_2001_1.csv',
            path2EQ+'2001_2001_2.csv',
            path2EQ+'2002_2002_0.csv',
            path2EQ+'2002_2002_1.csv',
            path2EQ+'2002_2002_2.csv',
            path2EQ+'2003_2003_0.csv',
            path2EQ+'2003_2003_1.csv',
            path2EQ+'2003_2003_2.csv',
            path2EQ+'2004_2004_0.csv',
            path2EQ+'2004_2004_1.csv',
            path2EQ+'2004_2004_2.csv',
            path2EQ+'2005_2005_0.csv',
            path2EQ+'2005_2005_1.csv',
            path2EQ+'2005_2005_2.csv',
            path2EQ+'2006_2006_0.csv',
            path2EQ+'2006_2006_1.csv',
            path2EQ+'2006_2006_2.csv',
            path2EQ+'2007_2007_0.csv',
            path2EQ+'2007_2007_1.csv',
            path2EQ+'2007_2007_2.csv',
            path2EQ+'2008_2008_0.csv',
            path2EQ+'2008_2008_1.csv',
            path2EQ+'2008_2008_2.csv',
            path2EQ+'2009_2009_0.csv',
            path2EQ+'2009_2009_1.csv',
            path2EQ+'2009_2009_2.csv',
            path2EQ+'2010_2010_0.csv',
            path2EQ+'2010_2010_1.csv',
            path2EQ+'2010_2010_2.csv',
            path2EQ+'2011_2011_0.csv',
            path2EQ+'2011_2011_1.csv',
            path2EQ+'2011_2011_2.csv',
            path2EQ+'2012_2012_0.csv',
            path2EQ+'2012_2012_1.csv',
            path2EQ+'2012_2012_2.csv',
            path2EQ+'2013_2013_0.csv',
            path2EQ+'2013_2013_1.csv',
            path2EQ+'2013_2013_2.csv',
            path2EQ+'2014_2014_0.csv',
            path2EQ+'2014_2014_1.csv',
            path2EQ+'2014_2014_2.csv',
            path2EQ+'2015_2015_0.csv',
            path2EQ+'2015_2015_1.csv',
            path2EQ+'2015_2015_2.csv',
            path2EQ+'2016_2016_0.csv',
            path2EQ+'2016_2016_1.csv',
            path2EQ+'2016_2016_2.csv',
            path2EQ+'2017_2017_0.csv',
            path2EQ+'2017_2017_1.csv',
            path2EQ+'2017_2017_2.csv',
            
            path2EQ+'2018_2018_0.csv',
            path2EQ+'2018_2018_1.csv',
            path2EQ+'2018_2018_2.csv',
            path2EQ+'2018_2018_3.csv',
            path2EQ+'2018_2018_4.csv',
            path2EQ+'2018_2018_5.csv',
            ]

path += [path2EQ, path2CR, path2Sun]
EQDateFormat = "%Y-%m-%dT%H:%M:%S.%fZ"

def sec2fraction(timestamp):
    """ conversion from timestamp to year fraction """
    dt = datetime.fromtimestamp(timestamp) # convert timestamp to datetime object
    year = dt.year # get the year
    fraction = dt.timetuple().tm_yday/365 # get the dat of the year
    return year + fraction

""" conversion from seconds to date (in string) """        
sec2date = lambda arg, formatS=EQDateFormat: datetime.fromtimestamp(int(arg)).strftime(formatS)

""" convert date (in string) to seconds """
date2sec = lambda arg, formatS=EQDateFormat: int(datetime.strptime(arg, formatS).timestamp())

""" Calculate probabilty N >= k """
Pcdf = lambda n, k: factorial(n)/(factorial(n-k)*factorial(k)*2**n) + Pcdf(n, k+1) if k <= n else 0
    


class SunClass:
    timeRecorded = []
    spots = []
    
    def __init__(self, path, startDate = '1970-01-01T00:00:00.00Z', stopDate = '2023-08-30T15:17:00.00Z'):
        data = DictReader(open(path))
        startDate = date2sec(startDate)
        stopDate = date2sec(stopDate)
        
        for row in data:
            if int(row['Year']) >= 1970: # warunek żeby ominąć błąd
                timestamp = date2sec(row['Year']+'-'+row['Month']+'-'+row['Day']+'T12:00:00.00Z')
                if timestamp >= startDate and timestamp < stopDate and row['sunspots'] != '-1':
                    self.timeRecorded.append(timestamp)
                    self.spots.append(int(row['sunspots']))
        
        self.timeRecorded = np.array(self.timeRecorded)
        self.spots = np.array(self.spots)

class EQClass:
    timeRecorded = []
    mag = []
    # magType = []
    
    def __init__(self, path, startDate = '1970-01-01T00:00:00.00Z', stopDate = '2023-08-30T15:17:00.00Z'):
        startDate = date2sec(startDate)
        stopDate = date2sec(stopDate)
        
        for i in path:
            data = DictReader(open(i))
            for row in data:
                timestamp = date2sec(row['time'])
                if timestamp >= startDate and timestamp < stopDate:
                    self.timeRecorded.append(timestamp)
                    self.mag.append(float(row['mag']))
                    # self.magType.append(row['magType']) # chyba jednak się nie przyda
                
        # self.timeRecorded = np.array(self.timeRecorded[::-1]) # W pliku csv z USGS dane są od najnowszych dlatego odwracam, ale ze skryptu od sławka wszystko jest git.
        # self.mag = np.array(self.mag[::-1])
        self.mag = np.array(self.mag)
        self.timeRecorded = np.array(self.timeRecorded)

        
class CRClass:
    timeRecorded = []
    rateCorr = []
    
    def __init__(self, path, startDate='1970-01-01T00:00:00.00Z', stopDate='2023-08-30T15:17:00.00Z'):
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



        
def makeWin42vars(EQtimeShifted, CR, binWidth, binStarts, binStops): 
    toRetEQ, toRetCR = [], []
    
    CR_start_indx = bisect_left(CR.timeRecorded, binStarts[0] - binWidth)
    CR_stop_indx = bisect_right(CR.timeRecorded, binStops[0] - binWidth)
    meanCR = np.mean(CR.rateCorr[CR_start_indx: CR_stop_indx])

    for begin,end in zip(binStarts, binStops):
        
        CR_start_indx = bisect_right(CR.timeRecorded, begin)
        CR_stop_indx = bisect_right(CR.timeRecorded, end)
        
        EQ_start_indx = bisect_right(EQtimeShifted, begin)
        EQ_stop_indx = bisect_right(EQtimeShifted, end)
        
        oldMeanCR = meanCR
        meanCR = np.mean( CR.rateCorr[ CR_start_indx: CR_stop_indx] )
        sumaEQ = np.sum( EQData.mag[ EQ_start_indx: EQ_stop_indx] )
        
        if not np.isnan(meanCR) and not np.isnan(sumaEQ) and not np.isnan(oldMeanCR):
            toRetCR.append(meanCR - oldMeanCR)
            toRetEQ.append(sumaEQ)
        
    toRetEQ, toRetCR = np.array(toRetEQ), np.array(toRetCR)
    medianEQ = np.median(toRetEQ)
    toRetCR = np.abs(toRetCR)
    medianCR = np.median(toRetCR)
    
    
    indices2remove = (toRetEQ == medianEQ) & (toRetCR == medianCR)
    toRetEQ = toRetEQ[~indices2remove]
    toRetCR = toRetCR[~indices2remove]
    
    return np.where( ( toRetEQ-medianEQ )*( toRetCR-medianCR ) > 0, 1, 0)


def dayShiftProbability(shift):        
    P           =   1675*24*60*60 # legth of a whole period over which the mean values will be gathered.
    d           =   5*24*60*60 # length of a window over which the sum / mean value will be calculated
    
    shift       *=  24*60*60
    shiftLen    =   len(shift)
    
    time        =   np.linspace(-20,20,32001)
    timeLen     =   len(time)
    
    probability =   np.zeros(timeLen)
    
    res         =   np.zeros(shiftLen)
    
    for j in range(timeLen):
        date = date2sec('2014-04-02T22:07:12.00Z') + time[j]*d
        binStarts   =   np.arange( date, date + P, d)
        binEnds     =   binStarts + d
        print(f"{j*100/timeLen:.2f}% completed.")
        
        for i in range(shiftLen):
            czas = EQData.timeRecorded + shift[i]
            
            c = makeWin42vars(czas, CRData, d, binStarts, binEnds)
            
            # print(f"{shift[i]/(24*60*60):.1f}. medianEQ = {medianEQ:.2f} \t medianCR = {medianCR:.2f}")
            
            res[i] = Pcdf(len(c), sum(c))
            
        probability[j] = np.min(res)
        
    shift = shift.astype('float') / (24*60*60)
    
    minimasSpotted = 41 # BARDZO WAŻNA LINIJKA!!! tj. ilość minimów.
    
    timeWins = np.linspace(0, timeLen-1, num=minimasSpotted+1, dtype=int) # indeksy między którymi będą minima (muszą być)
    minimas = np.zeros(minimasSpotted) # bo 5 minimów mamy, tutej
    
    for i in range(minimasSpotted):
        tempProb = probability[timeWins[i]:timeWins[i+1]] # mam nadzieję że id(temp) == id(probability)
        tempTime = time[timeWins[i]:timeWins[i+1]]
        minVal = np.min(tempProb)
        
        indices = np.where(tempProb == minVal)[0]
        
        minimas[i] = np.mean(tempTime[indices]) # hope it indicates the time when something happened.
        
    plt.title("Minimum probabilities with dependency\nto the number of days EQ data was shifted")
    plt.plot(time, np.log10(probability))
    plt.xlabel("Number of days that starting point was shifted")
    plt.ylabel(r"$log_{10}(P_{cdf})$")
    
    return minimas

def correlationIndicator(daysBack, shift, mode='Pearson'):
    czasCR = CRData.timeRecorded + shift # przesuwanie czasu w którym rateCorr zostało zaobserwowane
    meanRateCorrs, medianRateCorrs, magnitudes = [], [], []

    
    for i in range(len(EQData.timeRecorded)):
        rateCorrs = CRData.rateCorr[((czasCR < EQData.timeRecorded[i]) & (czasCR > EQData.timeRecorded[i] - daysBack))] # dane CR sprzed daysBack dni
        if len(rateCorrs) != 0:
            magnitudes.append(EQData.mag[i])
            meanRateCorrs.append(np.mean(rateCorrs))
            medianRateCorrs.append(np.median(rateCorrs))

    # tutaj leci funkcja jakiej korelacje chcemy policzyć.
    if mode=='Pearson':
        return np.corrcoef(meanRateCorrs, magnitudes)[0,1], np.corrcoef(medianRateCorrs, magnitudes)[0,1]

    elif mode=='Spearman':
        return spearmanr(meanRateCorrs, magnitudes)[1], spearmanr(medianRateCorrs, magnitudes)[1]

def corrDependency(shift, daysBack):
    shift       *= 24*60*60
    daysBack    *= 24*60*60
    shiftLength = len(shift)
    # daysBackLength = len(daysBack)
    meanCorr, medianCorr = np.zeros(shiftLength), np.zeros(shiftLength)

    for j in range(shiftLength):
        print(f"{(j/shiftLength)*100:.1f}% completed")
        meanCorr[j], medianCorr[j] = correlationIndicator(daysBack, shift[j], 'Spearman')

    shift /= 24*60*60
    daysBack /= 24*60*60

    return meanCorr, medianCorr

SunData = SunClass(path2Sun)
EQData = EQClass(path2EQ)
CRData = CRClass(path2CR)


# minimas = dayShiftProbability(np.arange(-16,-14, 0.1))

####### DO LICZENiA ŚREDNICH W CIĄGU DNiA
secsInDay       = .99915*24*60*60 # potem się zobaczy dla sidereal day'a
parts           = 24*4 # co 5 minut
magSumInDay     = [[] for _ in range(parts)] # deklaracja i alokacja array'a 2D
rateCorrs       = [[] for _ in range(parts)]
partsOfTheDay   = np.linspace(0, secsInDay, num=parts)

binStart = date2sec('2006-01-01T00:00:00.00Z')

CRData.rateCorr = np.diff(CRData.rateCorr)
CRData.timeRecorded = CRData.timeRecorded[:-1]

while binStart < date2sec('2008-01-01T00:00:00.00Z'):
    print(sec2fraction(binStart))
    for j in range(parts-1): # ale to się długo będzie musiało robić, lol 
        magSumInDay[j].extend(EQData.mag[((EQData.timeRecorded < binStart + partsOfTheDay[j+1]) & (EQData.timeRecorded > binStart))]) # na średnią magnitud
        condition = ((CRData.timeRecorded < binStart + partsOfTheDay[j+1]) & (CRData.timeRecorded > binStart))
        rateCorrs[j].extend(CRData.rateCorr[condition]) # na średnią rateCorr
    # magSumInDay[j] = np.sum(np.where((( EQData.timeRecorded > partsOfTheDay[j] + binStart ) & (EQData.timeRecorded < partsOfTheDay[j+1] + binStart )), 1, 0)) # na liczbę magnitud
    # magSumInDay[j] = np.sum(EQData.mag[ (( EQData.timeRecorded > partsOfTheDay[j] + binStart ) & (EQData.timeRecorded < partsOfTheDay[j+1] + binStart )) ]) # na sumę magnitud
    binStart += secsInDay

partsOfTheDay /= 60*60

print(rateCorrs)

rateCorrs   = [ np.mean(i) for i in rateCorrs]
magSumInDay = [ np.mean(i) for i in magSumInDay]

fig, ax1 = plt.subplots()

plt.title(f'Chart for the mean of the magnitudes in certain part of the day. Dividing the day on {parts} equal parts. Step is equal to 15 minutes. Start measuring on 01 Jan. 2000')

color = 'tab:green'
ax1.set_xlabel('Hour of the day')
ax1.set_ylabel('Sum of the magnitudes', color=color)
ax1.plot(partsOfTheDay[:-1], magSumInDay[:-1], linewidth=0, marker='o', label='Średnia magnitud', color=color)

ax2 = ax1.twinx()

color = 'tab:red'
ax2.plot(partsOfTheDay[:-1], rateCorrs[:-1], linewidth=0, marker='o', label='Średnia rateCorr', color=color)
ax2.set_ylabel('Mean value of dCR', color=color)
ax1.legend(loc='upper right')
ax2.legend(loc='lower left')



###### to jest do liczenia korelacji:   
# shift       =      np.arange(0, 15, 0.1)
# daysBack    =      1

# EQData.timeRecorded = EQData.timeRecorded[(EQData.mag > 4)]
# EQData.mag = EQData.mag[(EQData.mag > 4)]

# meanCorr, medianCorr = corrDependency(shift, daysBack)

# # wykres korelacji
# plt.plot(shift, meanCorr, label='Correlation for mean value of rateCorr')
# plt.plot(shift, medianCorr, label='Correlation for median value of rateCorr')
# plt.title(r'Pearson correlation between mean value of $rateCorr$ and $m^3$')
# # plt.xlabel('The number of days the magnitude\ndata was shifted forward')
# # plt.ylabel('Number of days the CR data was gathered\nbefore the days of recorded magnitude')
# plt.xlabel('The number of days the magnitude\ndata was shifted forward')
# plt.ylabel('Correlation')
# plt.legend(loc='upper right')

# print(f"Biggest correlation observed is equal to {np.amax(meanCorr)}")
# print(f"Smallest correlation observed is equal to {np.min(meanCorr)}")


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
