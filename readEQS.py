import csv
from datetime import datetime
import numpy as np
path2eqs = "/home/stud/praktyki/Data_analysis/EQs/"
path2eqs += "query.csv"
path2pao = "/home/stud/praktyki/Data_analysis/pierre_auger/"
path2pao += "scalers.csv"
defFormt = "%Y-%m-%dT%H:%M:%S.%fZ"

def sec2date(arg, formatS=defFormt):
    """ conversion from seconds to date (in string) """
    return datetime.fromtimestamp(int(arg)).strftime(formatS)


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
    
mag = mag[::-1]
timeEQs = timeEQs[::-1]
    
#### Koniec czytania trzęsień ziemi teraz czytanie pierre augiera    
auger = csv.DictReader(open(path2pao, encoding='utf-8'))

timeAug = []

for row in auger:
    timeAug.append(int(row['time']))
#### Koniec czytania pierre augera, liczenie korelacji teraz
timeEQs = np.array(timeEQs)
timeAug = np.array(timeAug)

def calcCorr(days, minutes):
    shift = int(60*60*24*days) # 1 dzień shifta do przodu dla Augiera
    global timeAug
    timeAug = timeAug + shift
    
    binWidth = int(minutes*60) # szerokość bina w sekundach
    minTime = np.maximum(np.min(timeAug), np.min(timeEQs))
    maxTime = np.minimum(np.amax(timeAug), np.amax(timeEQs))
    maxTime += maxTime%binWidth # musi być podzielny przez minutes żeby liczba okien była całkowita
    augWin = [] # okienko czasowe dla augera, zawiera liczbe zarejestrowanych promieniowań
    eqsWin = [] # okienko dla trzęsień ziemi, największe trzęsienia ziemi w okienkach.
    temp = 0
    oldTemp = 0
    maxBinVal = minTime + binWidth
    
    while minTime < maxTime:
        augWin.append(np.count_nonzero((minTime <= timeAug) & (timeAug < maxBinVal)))
        temp += np.count_nonzero((minTime <= timeEQs) & (timeEQs < maxBinVal))
        if oldTemp != temp:
            eqsWin.append(np.amax(mag[oldTemp:temp]))
        else:
            eqsWin.append(0) # w sumie to to jest nieprawda najlepiej byłoby to 0 pominąć, ale rozmiary muszą się zgadzać dlatego muszę je zostawić
        oldTemp = temp
        maxBinVal += binWidth
        minTime += binWidth
        
    return (np.corrcoef(augWin, eqsWin))[0,1] # np.corrcoef returns a matrix [0][1] points out to corr(A,B) not corr(A,A) nor corr(B,B)      

days = np.linspace(4, 6, num=21) # shift of Auger data forward in time in days,
binWidth = np.linspace(105, 110, num=51) # binWidth in minutes
coef = np.zeros([len(days),len(binWidth)])

for i in range(len(days)):
    for j in range(len(binWidth)):
        print(f"Starting loop for: day = {i} day and binWidth = {binWidth[j]}")
        coef[i,j] = calcCorr( days[i], binWidth[j] )