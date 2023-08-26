import sys
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

for row in eqInfo:
    timeEQs.append(int(date2sec(row['time'])))
    mag.append(float(row['mag']))
    
#### Koniec czytania trzęsień ziemi teraz czytanie pierre augiera    
auger = csv.DictReader(open(path2pao, encoding='utf-8'))

timeAug = []

for row in auger:
    timeAug.append(int(row['time']))
#### Koniec czytania pierre augera, liczenie korelacji teraz
timeEQs = np.array(timeEQs)
timeAug = np.array(timeAug)

minutes = 60
minTime = np.min(timeAug) if np.min(timeAug) < np.min(timeEQs) else np.min(timeEQs)
maxTime = np.amax(timeAug) if np.amax(timeAug) > np.amax(timeEQs) else np.amax(timeEQs)
maxTime += maxTime%minutes # musi być podzielny przez minutes żeby liczba okien była całkowita
time = np.arange(minTime,maxTime, minutes*60) # krok co minutes*60
augWin = [] # okienko czasowe dla augera, zawiera liczbe zarejestrowanych promieniowań
eqsWin = [] # okienko dla trzęsień ziemi, największe trzęsienia ziemi w okienkach.
temp = 0
oldTemp = 0

for i in range(len(time)-1):
    augWin.append(np.count_nonzero((time[i] <= timeAug) & (timeAug < time[i+1])))
    temp += np.count_nonzero((time[i] <= timeEQs) & (timeEQs < time[i+1]))
    if oldTemp != temp:
        eqsWin.append(np.amax(mag[oldTemp:temp]))
    else:
        eqsWin.append(0)
    oldTemp = temp
    
np.corrcoef(augWin, eqsWin)
