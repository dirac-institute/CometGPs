import os
import pandas as pd
import matplotlib.pyplot as plt

# ALCDEF data files can be downloaded from http://alcdef.org/ (just download the entire set of files, then unzip)
# ALCDEF data standard is described at http://alcdef.org/docs/ALCDEF_Standard.pdf

def read_alcdef(filename, datadir='.'):
    with open(os.path.join(datadir, filename)) as f:
        data = list(f)
    len(data)
    lc = []
    metadata = []
    delimiter = '|'
    for d in data:
        if d.startswith('STARTMETADATA'):
            mi = []
            li = []
        if d.startswith('ENDMETADATA'):
            metadata.append(mi)
            lc.append(li)
        if d.startswith('DATA'):
            tmpvals = d.rstrip('\n').lstrip('DATA=').split(delimiter)
            if tmpvals[2] == '':
                tmpvals[2] = -66
            if tmpvals[3] == '':
                tmpvals[3] = -66
            li.append([float(tmpvals[0]), float(tmpvals[1]), float(tmpvals[2]), float(tmpvals[3])])
        else:
            if d.startswith('STARTMETADATA') or d.startswith('ENDMETADATA') or d.startswith('ENDDATA'):
                continue
            mi.append(d.rstrip('\n'))
            if d.startswith('DELIMITER'):
                delimiter = d.split('=')[1]
                if delimiter.startswith('PIPE'):
                    delimiter = '|'
                elif delimiter.startswith('TAB'):
                    delimiter = '\t'
                elif delimiter.startswith('SPACE'):
                    delimiter = ' '
                elif delimiter.starstwith('COMMA'):
                    delimiter = ','
    nblocks = len(metadata)
    # Convert metadata segments into dictionaries
    for i in range(nblocks):
        md = {}
        for line in metadata[i]:
            vals = line.split('=')
            key = vals[0]
            val = vals[1]
            md[key] = val
        metadata[i] = md
    # Convert light curve segments into data frames
    for i in range(nblocks):
        lc[i] = pd.DataFrame(lc[i], columns=['JD', 'Mag', 'DeltaMag', 'Airmass'])
    return metadata, lc
    

def plot_lightcurve(metadata, lc, filename):
    colors = ['g', 'r', 'b', 'y']
    for i in range(len(metadata)):
        plt.figure(figsize=(12, 10))
        plt.errorbar(lc[i].JD - lc[i].JD[0], lc[i].Mag, yerr=lc[i].DeltaMag, linestyle='', marker='o', markersize=5, color=colors[i % 4])
        plt.title('%s: Start date %f' % (filename.rstrip('.txt'), lc[i].JD[0]))
        print('%s block %d: (midpoint %s, %s)' % (filename.rstrip('.txt'), i, metadata[i]['SESSIONDATE'], metadata[i]['SESSIONTIME']))
        print('Observed in %s, reporting %s' % (metadata[i]['FILTER'], metadata[i]['MAGBAND']))
        if metadata[i]['REDUCEDMAGS'] == 'AVERAGE':
            print('Magnitudes have been reduced to unity distance, corrected by %s (AVERAGE)' % metadata[i]['UCORMAG'])
        if metadata[i]['REDUCEDMAGS'] == 'POINT':
            print('Magnitudes have been reduced to unity distance, per measurement (POINT). Average is %s.' % metadata[i]['UCORMAG'])
        if metadata[i]['DIFFERMAGS']:
            print('Reported differential magnitudes')
            if metadata[i]['MAGADJUST'] != '0.0':
                print('A suggested correction of %s could place these onto a standard system' % metadata[i]['MAGADJUST'])
        else:
            print('Reported standard magnitudes, using %s refcat' % metadata[i]['STANDARD'])
        if metadata[i]['LTCAPP'] == 'AVERAGE':
            print('Times are light time corrected by %s (AVERAGE)' % (metadata[i]['LTCDAYS']))
        if metadata[i]['LTCAPP'] == 'POINT':
            print('Times are light time correct per measurement (POINT). Average is %s.' % (metadata[i]['LTCDAYS']))
    plt.figure(figsize=(24, 10))
    for i in range(len(metadata)):
        plt.errorbar(lc[i].JD, lc[i].Mag, yerr=lc[i].DeltaMag, linestyle='', marker='o', markersize=5, color=colors[i % 4])
