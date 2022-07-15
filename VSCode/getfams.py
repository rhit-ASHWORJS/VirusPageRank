from tkinter.filedialog import askopenfilename
import csv

famlist = []
filename = askopenfilename()
with open(filename) as csvfile:
    datareader = csv.reader(csvfile, delimiter=',')
    for row in datareader:
        if(row[1] not in famlist) and ('SINGLETON' not in row[1]):
            famlist.append(row[1])

print(famlist)