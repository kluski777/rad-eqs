import csv
from datetime import datetime
import numpy as np
path2eqs = "/home/stud/praktyki/Data_analysis/EQs/"
path2eqs += "query.csv"
path2pao = "/home/stud/praktyki/Data_analysis/pierre_auger/"
path2pao += "scalers.csv"
defFormt = "%Y-%m-%dT%H:%M:%S.%fZ"

# def sec2date(arg, formatS=defFormt):
#     """ conversion from seconds to date (in string) """
#     return datetime.fromtimestamp(int(arg)).strftime(formatS)
#

def date2sec(arg, formatS=defFormt):
    """ convert date (in string) to seconds """
    return int(datetime.strptime(arg, formatS).timestamp())

#### Tutaj czytam plik .csv z trzęsieniami ziemi
eqInfo= csv.DictReader(open(path2eqs, encoding='utf-8'))    

timeEQs = []
mag = []

for row in eqInfo: # tutaj się czyta od najnowszych do najstarszych odwrotnie jak z augera
    timeEQs.append(int(date2sec(row['time'])))
    mag.append(float(row['mag']))
    
mag = mag[::-1] # odwracam arraya bo się jakoś od tyłu czytają z tego pliku
timeEQs = timeEQs[::-1]
    
#### Koniec czytania trzęsień ziemi teraz czytanie pierre augiera    
auger = csv.DictReader(open(path2pao, encoding='utf-8'))

timeAug = []

for row in auger:
    timeAug.append(int(row['time']))
#### Koniec czytania pierre augera, liczenie korelacji teraz
timeEQs = np.array(timeEQs)
timeAug = np.array(timeAug)

def calcCorr(shift, minutes):
    global timeAug # jak tego uniknąć ten array jest ogromny?
    timeAug = timeAug + shift # ta linijka jest do bani
    
    binWidth = int(minutes*60) # szerokość bina w sekundach
    minTime = np.maximum(np.min(timeAug), np.min(timeEQs))
    maxTime = np.minimum(np.amax(timeAug), np.amax(timeEQs))
    # maxTime += maxTime%binWidth
    augWin = [] # liczba zarejstrowanych promieniowań w okienku czasowym
    eqsWin = [] # Największa magnituda w danym okienku.
    temp, oldTemp = 0, 0 # indeks dolny i górny
    maxBinVal = minTime + binWidth
    oldVal = 0
    
    while minTime < maxTime:
        temp += np.count_nonzero((minTime <= timeEQs) & (timeEQs < maxBinVal)) # chyba tą linijkę można poprawić żeby przyspieszyć skrypt
        if temp != oldTemp:
            eqsWin.append(np.amax(mag[oldTemp:temp]))
            augWin.append(np.count_nonzero((minTime <= timeAug) & (timeAug < maxBinVal)) - oldVal)
        oldVal = np.count_nonzero((minTime <= timeAug) & (timeAug < maxBinVal))
        oldTemp = temp
        maxBinVal += binWidth
        minTime += binWidth
        
    return (np.corrcoef(augWin, eqsWin))[0,1] # np.corrcoef returns a matrix [0][1] points out to corr(A,B) not corr(A,A) nor corr(B,B)      

days_start, days_stop = 0, 1
timeAug += days_start*24*60*60 # jak zaczynamy badań próbkę od dnia days_start to przesuwamy to
step = 0.01 # w dniach
length = int((days_stop-days_start)/step)
binWidth = np.linspace(1,150, num=150) # w minutach
coef = np.zeros([length,len(binWidth)])
step = int(step*24*60*60) # zamiana jednostek na sekundy

for i in range(length):
    for j in range(len(binWidth)):
        print(f"Starting loop for: day = {i} day and binWidth = {binWidth[j]}")
        coef[i] = calcCorr(step, binWidth[j])