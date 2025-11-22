# -*- coding: utf-8 -*-
""" 
Useful unit conversions for GeoDT
"""
#initializations
import numpy as np

def save_csv(out=[['null',0.0]],filename='save.csv',append=True):
    #output to file
    out = list(map(list,zip(*out)))
    head = out[0][0]
    for i in range(1,len(out[0])):
        head = head + ',' + out[0][i]
    try:
        data = '%i' %(out[1][0])
    except:
        data = '' + out[1][0]
    for i in range(1,len(out[1])):
        try:
            data = data + ',%.5e' %(out[1][i])
        except:
            data = data + ',' + out[1][i]
    if append:
        try:
            with open(filename,'r') as f:
                test = f.readline()
            f.close()
            if test != '':
                with open(filename,'a') as f:
                    f.write(data + '\n')
                f.close()
            else:
                with open(filename,'a') as f:
                    f.write(head + '\n')
                    f.write(data + '\n')
                f.close()
        except:
            with open(filename,'a') as f:
                f.write(head + '\n')
                f.write(data + '\n')
            f.close()
    else:
        with open(filename,'w') as f:
            f.write(head + '\n')
            f.write(data + '\n')
        f.close()
    
def save_set(out=[['null',0.0,0.0,0.0,'none']],filename='save.set'):
    #output to file in column format
    with open(filename,'w') as f:
        for i in range(0,len(out)):
            line = out[i][0]
            for j in range(1,len(out[0])):
                try:
                    line += ',' + '%.5e' %(out[i][j])
                except:
                    line += ',' + out[i][j]
            f.write(line + '\n')
    f.close()

def build_csv(header=[],data=[],filename='save.csv',append=True):
    print('data length: %i' %(len(data)))
    for i in range(0,len(data)):
        out = []
        for j in range(0,len(header)):
            out += [[header[j],data[header[j]][i]]]
            # if j == 0:
            #     print(data[header[j]][i])
        if i == 0:
            save_csv(out,filename,append)
        else:
            save_csv(out,filename,True)

def special_csv(header=[],data=[],source='',target='',filename='save.csv',append=True):
    for i in range(0,len(data)):
        out = []
        for j in range(0,len(header)):
            if (i == 1) and (header[j] == 'Strategy_SourceDirectory_path'):
                out += [[header[j],source]]
            elif (i == 1) and (header[j] == 'Strategy_TargetDirectory_path'):
                out += [[header[j],target]]
            else:
                out += [[header[j],data[header[j]][i]]]
        if i == 0:
            save_csv(out,filename,append)
        else:
            save_csv(out,filename,True)

def copy_csv(sourcefile='load.csv',targetfile='save.csv',append=True):
    #read sourcefile
    data = []
    try:
        with open(sourcefile,'r') as f1:
            data = f1.readlines()
            f1.close()
        if not(append):
            with open(targetfile,'w') as f2:
                f2.writelines(data)
                f2.close()
        else:
            with open(targetfile,'r') as f2:
                test = f2.readline()
                f2.close()
            with open(targetfile,'a') as f2:
                if test != '':
                    f2.writelines(data[1:])
                else:
                    f2.writelines(data[:])
                f2.close()
    except:
        print('failed to copy file %s to %s' %(sourcefile, targetfile)) 

def load_csv(filename='save.csv'):
    names = []
    data = []
    #load data as list of rows; note that genfromtxt doesn't work, so a custom function is needed
    with open(filename,'r') as f:
        file = f.readlines()
        f.close()
    #convert into arrays
    names = file[0].split(',')
    names[-1] = names[-1][:-1]
    for i in range(1,len(file)):
        data += [file[i][:-1].split(',')]
    #convert into structured array
    numbers = np.ones(len(names))
    dtype = []
    for j in range(0,len(names)):
        try:
            for i in range(0,len(data)):
                float(data[i][j])
            if (data[i][j].find('.') > -1) or (data[i][j] == "inf"):
                numbers[j] = 2
        except:
            numbers[j] = 0
        if numbers[j] == 2:
            dtype += [tuple([names[j],'<f8'])]
        elif numbers[j] == 1:
            dtype += [tuple([names[j],'<i8'])]
        else:
            dtype += [tuple([names[j],'<U600'])]
    for i in range(0,len(data)):
        data[i] = tuple(data[i])
    dtype = np.dtype(dtype)
    data = np.array(data,dtype=dtype)
    return names, data

def load_set(filename='save.set'):
    names = []
    data = []
    #load data as list of rows; note that genfromtxt doesn't work, so a custom function is needed
    with open(filename,'r') as f:
        file = f.readlines()
        f.close()
    #convert into arrays
    names = []
    data = []
    dtype = []
    for i in range(0,len(file)):
        hold = file[i][:-1].split(',')
        names += [hold[0]]
        data += [[hold[1],hold[2],hold[3],hold[4]]]
        dtype += [tuple([names[i],'<U600'])]
    # list(map(list,zip(*data)))
    data = [list(row) for row in zip(*data)]
    for i in range(0,len(data)):
        data[i] = tuple(data[i])
    dtype = np.dtype(dtype)
    data = np.array(data,dtype=dtype)
    return names, data

def convert(source='save.gts'):
    #load settings
    names, data = load_csv(source)
    out = []
    for i in range(0,len(names)):
        out += [[names[i], data[names[i]][0], data[names[i]][1], data[names[i]][2], data[names[i]][3]]]
    save_set(out,source[:-4]+'.set')