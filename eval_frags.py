# Written by Mukundan S mukundan.kollam@gmail.com

'''
Copyright (C) 2025 Mukundan S <mukundan.kollam@gmail.com>

You can use the contents of this file under the conditions of Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International deed (https://creativecommons.org/licenses/by-nc-sa/4.0/). 
A copy of this deed in provided in the LICENSE.md file.

You are free to copy, distribute, remix, transform and build upon this program as long as 1) it is for non commercial purpose, 2) Attributions are provided 3) Shared under the same license. 

Please cite https://doi.org/10.1101/2025.06.16.659877 if you use this program
'''

def load_csv_table_as_di(file):
    # read_file
    
    with open(file) as fp:
        lines = fp.readlines()
    di = {}  # {'NNNN':(relative efficenfy in percentage,percentage of correct reads aka fidility)}
    
    diagonal_pos = 1
    max_count = 0

    for line in lines[1:]:
        words = line.strip().split(',')
        if max_count< int(words[diagonal_pos]):
            max_count = int(words[diagonal_pos])
        diagonal_pos+=1

    diagonal_pos = 1

    for line in lines[1:]:# 0 is headder
        # its a summentric metrix so we only need to analyse one line. the same value can be entered fro the reverse complement
        words = line.strip().split(',')
        wc_pair = int(words[diagonal_pos])
        total = sum([int(x) for x in words[1:]])

        # print(total,max_count)

        relative_efficency =  (wc_pair/max_count)*100               
        fidility = (wc_pair/total)*100
        # print(e,words[0])

        di[words[0]] = (relative_efficency,fidility)

        diagonal_pos+=1 
    return di

def read_fasta(file):
    seqs = {}
    seq = ''
    headder = False
    with open(file) as inf:
        for line in inf:
            if line[0] == '>':
                if headder:
                    seqs[headder] = seq
                headder = line.strip()
                seq = ''
            else:
                seq+=line.strip()
    seqs[headder] = seq.upper()
    return seqs



import argparse

# parsing arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description='''This program takes in a list of 4 base pair overhang sequences in 5'-3' orientation and 
                                 evaluates the net score of efficency and net probability of correct assembly of those fragments. 
                                 The data used to calculate was published by Potapov et. al. 2018. 
                                 https://doi.org/10.1021/acssynbio.8b00333
                                 '''
                                 )
# # positional arguments 
parser.add_argument("overhangs", nargs='*', help="4 base pair sequences space seperated",type=str)
# optional argument1
parser.add_argument("-data_file", help="Supplimentarty file from Potapov et. al. 2018 to be used. They are in ./lib",type=str,default='./lib/FileS04_T4_18h_37C.csv')
parser.add_argument("-fasta_input", help="Fasta input from SplitSet",type=str,default=False)

a = parser.parse_args()

ovrhgs = a.overhangs

if a.fasta_input:
    seq_di = read_fasta(a.fasta_input)
    for x in seq_di:
        ovrhgs.append(x[-4:])

print('Fragments:',', '.join([str(x) for x in ovrhgs]))

score_di = load_csv_table_as_di(a.data_file)

net_eff,net_fid = 1,1
print('{:10s}    {:10s}    {:10s}'.format('Overhang','Efficency(%)','Fidility(%)'))
for ovrhg in ovrhgs:
    try:
        eff,fid = score_di[ovrhg]
    except KeyError:
        continue
    print('{:10s}     {:8.2f}     {:8.2f}'.format(ovrhg,round(eff),round(fid)))
    net_eff = net_eff*(eff/100)
    net_fid = net_fid*(fid/100)

print('\nNet efficency = {}'.format(net_eff*100))#,round(net_eff*100),'%')
print('Net fidility = {:8.2f}'.format(net_fid*100))#,round(net_fid*100),'%')
#print('\nNet efficency = ',round(net_eff*100),'%')

#print('Net fidility = ',round(net_fid*100),'%')
