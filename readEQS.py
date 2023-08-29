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
    # return datetime.fromtimestamp(int(arg)).strftime(formatS)

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

def calcCorr(shift, nDays):
    
    nDays *= 24*60*60
    shift *= 24*60*60
    augWin = np.zeros(len(timeEQs))
    
    for i in range(len(timeEQs)): # dodawać pochodną liczby rayów.
        augWin[i] = np.count_nonzero((timeAug > timeEQs[i] - nDays + shift) & (timeAug < timeEQs[i] + shift)) # zmiana liczby rayów
    
    return (np.corrcoef(augWin, mag))[0,1] # np.corrcoef returns a matrix [0][1] points out to corr(A,B) not corr(A,A) nor corr(B,B)      

#### argumenty do funkcji i korelacje
# shift = np.arange(0,10,0.2)
nLastDays = np.arange(0.01, 0.02,0.001)
coef = np.zeros(len(nLastDays))
augWin = np.zeros([len(timeEQs), len(nLastDays)])

# for i in range(len(shift)):
for j in range(len(nLastDays)):
    print(f"Starting loop for: the day {nLastDays[j]}.")
    coef[j] = calcCorr(0, nLastDays[j])