import pandas as pd
import tkinter as tk
import numpy as np
import pyexcel as p
import sys
import math
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from ortools.sat.python import cp_model

#required before we can ask for input file
root = tk.Tk()
root.withdraw()

gui = Tk()
gui.geometry("400x150")
gui.title("UT Mat Order")
gui.columnconfigure(0, weight=3)
gui.columnconfigure(1, weight=1)

folderPath = StringVar()
filePath = StringVar()

class FolderSelect(Frame):
    def __init__(self,parent=None,folderDescription="",**kw):
        Frame.__init__(self,master=parent,**kw)
        self.folderPath = StringVar()
        self.lblName = Label(self, text=folderDescription)
        self.lblName.grid(row=0,column=0, sticky="ew", pady=1)
        self.entPath = Entry(self, textvariable=self.folderPath)
        self.entPath.grid(row=1,column=0, sticky="ew", pady=1)
        self.btnFind = ttk.Button(self, text="Select Folder",command=self.setFolderPath)
        self.btnFind.grid(row=1,column=1, pady=1)
    def setFolderPath(self):
        folder_selected = filedialog.askdirectory()
        self.folderPath.set(folder_selected)
        self.entPath.insert(0,folder_selected)
    @property
    def folder_path(self):
        self.entPath.update()
        return self.folderPath.get()
    
class FileSelect(Frame):
    def __init__(self,parent=None,folderDescription="",**kw):
        Frame.__init__(self,master=parent,**kw)
        self.filePath = StringVar()
        self.lblName = Label(self, text=folderDescription)
        self.lblName.grid(row=0,column=0, sticky="ew", pady=1)
        self.entPath = Entry(self, textvariable=self.filePath)
        self.entPath.grid(row=1,column=0, sticky="ew", pady=1)
        self.btnFind = ttk.Button(self, text="Select File",command=self.setFilePath)
        self.btnFind.grid(row=1,column=1, pady=1)
    def setFilePath(self):
        file_selected = filedialog.askopenfilename()
        self.filePath.set(file_selected)
        self.entPath.insert(0,file_selected)
    @property
    def file_path(self):
        self.entPath.update()
        return self.filePath.get()
        

def doStuff():
    global excel_file
    global output_directory
    excel_file = file1Select.file_path
    output_directory = directory1Select.folder_path
    root.quit()

def endProgram():
    sys.exit()


file1Select = FileSelect(gui,"Excel BOM File:")
file1Select.grid(row=0)

directory1Select = FolderSelect(gui,"Order File Output Folder:")
directory1Select.grid(row=1)

c = ttk.Button(gui, text="RUN", command=doStuff)
c.grid(row=4,column=0, pady=1)
e = ttk.Button(gui, text="EXIT", command=endProgram)
e.grid(row=4,column=1, pady=1)
gui.mainloop()

if excel_file[len(excel_file)-1] == "s":
        p.save_book_as(file_name=excel_file,
               dest_file_name=excel_file + "x")
        excel_file = excel_file + "x"

###
#  GET LIBRARY DATA AND ONCOR ID / PROJECTS DATA
###
df = pd.read_excel(excel_file, sheet_name="LINEAR LIBARY (ORDER)", header=[1], dtype_backend="pyarrow")
dfDict = pd.read_excel(excel_file, sheet_name="ONCOR I.D.'S", header=[0], dtype_backend="pyarrow")
dict = dfDict.set_index('ID').T.to_dict('series')

###
#  APPLY ID QUANTITIES AND PROJECT NAMES TO LIBRARY
###
outputList = []
for key in dict:
    tempDF = df.copy(deep=True)
    tempDF['ID QTY'] = tempDF['ID'].map(dict[key])
    tempDF['TOTAL QTY'] = tempDF['QTY PER ID'] * tempDF['ID QTY']
    tempDF = tempDF.drop(tempDF[tempDF['TOTAL QTY'] == 0].index)
    tempDF['TOTAL LENGTH TOTAL QTY'] = tempDF['PART LENGTH PER ID'] * tempDF['ID QTY']
    tempDF['PERCENT USED'] = tempDF['TOTAL LENGTH TOTAL QTY'] / tempDF['STOCK LENGTH']
    tempDF['ROUNDED'] = tempDF['PERCENT USED'].apply(lambda x:(math.ceil(x*4)/4))
    tempDF['PARTS PER STICK'] = tempDF.apply(lambda row:(math.floor(row['STOCK LENGTH']/row['LENGTH2'])), axis=1)
    tempDF['MINIMUM STOCK'] = tempDF['TOTAL QTY'] / tempDF['PARTS PER STICK']
    tempDF['DROP SIZE'] = 0
    tempDF['DROP QTY'] = 0
    maskDF = tempDF[tempDF['PARTS PER STICK'] == 1]
    tempDF = tempDF[tempDF['PARTS PER STICK'] != 1]
    maskDF['DROP SIZE'] = maskDF['STOCK LENGTH'] - maskDF['LENGTH2']
    maskDF['DROP QTY'] = maskDF['MINIMUM STOCK']
    tempDF['PROJECT'] = key
    maskDF['PROJECT'] = key
    outputList.append(tempDF)
    outputList.append(maskDF)
outputDF = pd.concat(outputList, ignore_index=True)
outputDF.to_excel(output_directory + "//" + " dfOrder.xlsx", sheet_name="Sheet 1")

###
#  PREP DATA AND NEST USING CP-SAT SOLVER
###
nestDF = outputDF.copy(deep=True)
nestDF['LENGTH2'] = nestDF['LENGTH2'].apply(lambda x: x*1000)
nestDF = nestDF.loc[nestDF.index.repeat(nestDF['TOTAL QTY'])].reset_index(drop=True)

NestWorksetDataFrame = []
CutTicketWorksetDataFrame = []

def create_data_model():
    data = {}
    # part lengths
    data['weights'] = dfMatType['LENGTH2'].astype(int).values.tolist()
    data['items'] = list(range(len(data['weights'])))
    data['bins'] = data['items']
    # stick size
    data['bin_capacity'] = dfMatType.iloc[0, 12] * 1000
    data['material'] = dfMatType.iloc[0, 1]
    data['stock-code'] = dfMatType.iloc[0, 0]
    return data

# angle nesting function
for group, dfMatType in nestDF.groupby(['STOCK CODE', 'STOCK LENGTH']):
    data = create_data_model()

    # Create the CP-SAT model.
    model = cp_model.CpModel()

    # Variables
    # x[i, j] = 1 if item i is packed in bin j.
    x = {}
    for i in data['items']:
        for j in data['bins']:
            x[(i, j)] = model.NewIntVar(0, 1, 'x_%i_%i' % (i, j))

    # y[j] = 1 if bin j is used.
    y = {}
    for j in data['bins']:
        y[j] = model.NewIntVar(0, 1, 'y[%i]' % j)

    # Constraints
    # Each item must be in exactly one bin.
    for i in data['items']:
        model.Add(sum(x[i, j] for j in data['bins']) == 1)

    # The amount packed in each bin cannot exceed its capacity.
    for j in data['bins']:
        model.Add(
            sum(x[(i, j)] * data['weights'][i] for i in data['items']) <= y[j] *
            data['bin_capacity'])

    # Objective: minimize the number of bins used.
    model.Minimize(sum(y[j] for j in data['bins']))

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = 120.0
    status = solver.Solve(model)

    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        num_bins = 0
        bin_usage = 0
        for j in data['bins']:
            if solver.Value(y[j]) == 1:
                bin_items = []
                bin_weight = 0
                for i in data['items']:
                    if solver.Value(x[i, j]) > 0:
                        bin_items.append(i)
                        # stick usage
                        bin_weight += data['weights'][i]
                        CutTicketDictionary = {'PART': dfMatType.iloc[i,2], 'QTY': 1, 'STOCK CODE': data['stock-code'], 'MAT. DESCRIPTION': data['material'], 'LENGTH': dfMatType.iloc[i,4], 'ID': dfMatType.iloc[i,7], 'NESTED LENGTH': (data['weights'][i])/1000, 'STICK': j}
                        # List of parts to dataframe
                        CutTicketDictionaryDataFrame = pd.DataFrame(data=CutTicketDictionary, index=[0])
                        # Add the parts to the overall list
                        CutTicketWorksetDataFrame.append(CutTicketDictionaryDataFrame)
                if bin_items:
                    # counting number of sticks pulled
                    num_bins += 1
        # make list of parts
        NestDictionary = {'STOCK CODE': data['stock-code'], 'ORDER': num_bins}
        # list to dataframe
        NestDictionaryDataFrame = pd.DataFrame(data=NestDictionary, index=[0])
        # add parts to overall list
        NestWorksetDataFrame.append(NestDictionaryDataFrame)
    else:
        # there's either a fatal problem, or there's too many "good" solutions
        print('Nesting problem does not have an optimal or feasible solution.')

#new excel file
PostNestDataFrame = pd.concat(NestWorksetDataFrame, ignore_index=True)
PostNestDataFrame.to_excel(output_directory + "//" + " NestOrder.xlsx", sheet_name="Sheet 1")

CutTicketDataFrame = pd.concat(CutTicketWorksetDataFrame, ignore_index=True)
#combining multiple quantities of the same part on the same stick
CutTicketDataFrame = CutTicketDataFrame.groupby(['PART', 'STOCK CODE', 'MAT. DESCRIPTION', 'LENGTH', 'NESTED LENGTH', 'ID', 'STICK'])['QTY'].sum(numeric_only=True).reset_index()
CutTicketDataFrame['SHOP NOTES'] = "CUT " + CutTicketDataFrame['QTY'].apply(str) + " PCS @ " + CutTicketDataFrame['LENGTH']
CutTicketDataFrame['RAW MAT QTY'] = None
CutTicketDataFrame['HEAT NUMBER'] = None
CutTicketDataFrame['LOCATION'] = None
CutTicketDataFrame['GRADE'] = CutTicketDataFrame['MAT. DESCRIPTION'].str.split("'").str[-1]
CutTicketDataFrame = CutTicketDataFrame[['PART', 'QTY', 'STOCK CODE', 'GRADE', 'MAT. DESCRIPTION', 'RAW MAT QTY', 'HEAT NUMBER', 'LOCATION', 'SHOP NOTES', 'LENGTH', 'ID', 'NESTED LENGTH', 'STICK']]
CutTicketDataFrame.to_excel(output_directory + "//" + " DEBUGCutTicket.xlsx", sheet_name="Sheet 1")
writer = pd.ExcelWriter(output_directory + "//" + " CutTicket.xlsx")
for group, MaterialGroup in CutTicketDataFrame.groupby(['STOCK CODE']):
    MaterialGroup = MaterialGroup.sort_values(by=['STICK'])
    MaterialGroup.to_excel(writer, sheet_name=MaterialGroup.iloc[0,2])
writer.close()

###
#  GET ORDER QUANTITY PER MATERIAL ACROSS ALL PROJECTS
###
dfOrder = outputDF.groupby(['STOCK CODE'])[['ROUNDED', 'MINIMUM STOCK']].sum(numeric_only=True).reset_index()
dfOrder['STICKS'] = dfOrder[['ROUNDED', 'MINIMUM STOCK']].max(axis=1)
dfOrder['+0.01'] = dfOrder['STICKS'] + 0.01
dfOrder['QUANTITY ORDER'] = np.ceil(dfOrder['+0.01'])
dfOrder = dfOrder.drop('ROUNDED', axis=1)
dfOrder = dfOrder.drop('STICKS', axis=1)
dfOrder = dfOrder.drop('MINIMUM STOCK', axis=1)
dfOrder = dfOrder.drop('+0.01', axis=1)
dfOrder.to_excel(output_directory + "//" + " UT Linear Order.xlsx", sheet_name="Sheet 1")

###
#  GET ALLOCATION QUANTITY PER PROJECT
###
writerAllocate = pd.ExcelWriter(output_directory + "//" + " UT Linear Allocation.xlsx")
dfAllocate = outputDF.groupby(['STOCK CODE', 'PROJECT'])['ROUNDED'].sum(numeric_only=True).reset_index()
for group, StationGroup in dfAllocate.groupby(['PROJECT']):
    StationGroup = StationGroup[['STOCK CODE', 'ROUNDED', 'PROJECT']]
    StationGroup.rename(columns = {'ROUNDED':'ALLOCATE'}, inplace=True)
    StationGroup.to_excel(writerAllocate, sheet_name=StationGroup.iloc[0,2])
writerAllocate.close()