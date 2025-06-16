# this program creates fragmants with overhangs which have optimal ligation efficencyes values

# run with python3 with -h option to see documentation

# 

'''
Copyright (C) 2025 Mukundan S <mukundan.kollam@gmail.com>

You can use the contents of this file under the conditions of Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International deed (https://creativecommons.org/licenses/by-nc-sa/4.0/). 
A copy of this deed in provided in the LICENSE.md file.

You are free to copy, distribute, remix, transform and build upon this program as long as 1) it is for non commercial purpose, 2) Attributions are provided 3) Shared under the same license. 
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

class Dyna_frag:
    
    def __init__(self,seq,min_len,max_len,overhang_table,n_frag=False,n_trace = 5, start_site = 0, eff_w = 1,fid_w = 1, exclusion_list=[],inti_score = 1,max_overhang_identity = 2, alien_overhangs = []) -> None:
        '''
        seq = DNA sequence (can have non ATGC residues, they are ignored and overlaps are not made)
        '''
        self.seq = seq
        self.min_len = min_len
        self.max_len = max_len
        self.start_site = start_site
        self.eff_w = eff_w
        self.fid_w = fid_w
        self.n_trace = n_trace
        self.exclusion_list = exclusion_list
        self.inti_score = inti_score # scores are in 100 so 1 shoulnt effect anything while giving true
        self.max_overhang_identity = max_overhang_identity
    
        self._5overhang = ''
        self._3overhang = ''

        self.alien_overhangs = alien_overhangs

        # load csv file as a dictionary
        # {'NNNN':(relative efficency in percentage,percentage of correct reads aka fidility)}

        self.score_di = load_csv_table_as_di(overhang_table)
        print('Overhang scores loaded from ',overhang_table)

        # make matrix
        self.make_matrix()
        print('Matrix populated')
        
        # do traceback
        self.traceback(n_frag)
        print('Trace back completed')

    def make_matrix(self):
        '''
        K = maximum number of fragments = len(seq)/min fragment len
        for i in range(K):
            for j in range(len(seq)):
                score(i,j) = max(processed(efficenfy score(i,j) and fidility score(i-1,j)))
                add trace to i-1,j that was max 
        traceback for top n combinations
        '''
        # maximum possible number of fragments 
        K = self.round_up(len(self.seq)/self.min_len)
                
        # create an matrix to iterate over        
        self.mat = [[False for i in range(len(self.seq))] for x in range(K+1)] # +1 for starting condition
        self.split_mat = [[False for i in range(len(self.seq))] for x in range(K+1)] # +1 for starting condition
        # fill intial conditions
        # ie set start positions as intial score
        self.mat[0][self.start_site] = self.inti_score
        self.split_mat[0][self.start_site] = (self.inti_score,self.inti_score)
        
        # dictionary to store traces
        self.trace_di = {} # {i_j:[i-1_j*,]}

        # compute scores
        for i in range(1,K+1):
            # find the range of possible j based on fragment length 
            min_j,max_j = self.get_j_range(i) 
            for j in range(len(self.seq)):
                # check if j is in the correct range, else skip
                if j < min_j or j > max_j or j in self.exclusion_list:
                    continue    

                # only j that is in window must reach here
                
                # score the overhang at the j
                ovrhg = self.seq[j:j+4]  ### modify here to add reverse orientaitons
                try: # skips non standard stuff and ends, also leves those position as False
                    eff,fid = self.score_di[ovrhg] 
                except KeyError as e:
                    continue  

                # score = self.eff_w*(eff/100) + self.fid_w*(fid*100)
                

                # find the j_ from i-1 state that gives maximum score identity penalty is also added here
                sum_scores = sorted(self.__get_score_list(i,j,eff,fid) , reverse= True, key=lambda x:x[0])
                if not sum_scores:
                    continue # if sum scores in empty it means that there were no valid traces for the J
                # print(sum_scores)
                score_tot,j_,eff_new,fid_new = sum_scores[0]
                # print(eff_new,fid_new)
                # j_ = sum_scores[0][1]
                # score_tot = sum_scores[0][0]
                # update matrix                
                self.mat[i][j] = score_tot
                self.split_mat[i][j] = (eff_new,fid_new)

                # update trace
                self.trace_di[self.make_key(i,j)] = self.make_key(i-1,j_)
            # ending the iteration through K if the range is

    def __get_score_list(self,i,j,eff,fid):
        score_li = []

        for j_ in range(len(self.seq)):
            frag_len = j-j_ 
            if frag_len >= self.min_len and frag_len <= self.max_len and self.mat[i-1][j_] and self.__overlap_pass(i,j,j_): 
                # print(j_)                 
                eff_tally,fid_tally = self.split_mat[i-1][j_]

                eff_new = eff_tally*(eff/100)
                fid_new = fid_tally*(fid/100)

                score = (eff_new**self.eff_w)*(fid_new**self.fid_w)
                score_li.append((score,j_,eff_new,fid_new))
        return score_li

    def __overlap_pass(self,i,j,j_):
        # finds the previous motifs and checks if the motifs overlap
        js = self.__trace(i-1,j_)
        current_overhang = self.seq[j:j+4]
        overhangs = [self.seq[x:x+4] for x in js]+[self.seq[-4:]] + self.alien_overhangs # addition is the last overhang
        identities = self.__find_identities(overhangs+[self.get_rc(x) for x in overhangs+[current_overhang]],current_overhang)
        # last addition checks if the reverse complement in the other side of the cut will interfere with any of the ligations
        
        for x in identities:
            if x > self.max_overhang_identity:
                # print(overhangs,current_overhang,i,j,j_)
                return False
        return True

    def __find_identities(self,overhangs,current_overhang):
        ns = []
        for ovrhg in overhangs:
            ns.append(sum([1 if current_overhang[x] == ovrhg[x] else 0 for x in range(len(current_overhang))]))
        return ns

    def get_j_range(self,i):
        # min_j = min_j of i-1 + min_len
        # max_j = max_j of i-1 + max_len
        scored_prev_set = [x for x in range(len(self.seq)) if self.mat[i-1][x]]
        # print(scored_prev_set)        
        if len(scored_prev_set) == 0:
            min_j = len(self.seq)+1 # all comparison in the parent function will fail and skip
            max_j = len(self.seq)+2
        else:
            min_j = min(scored_prev_set) + self.min_len
            max_j = max(scored_prev_set) + self.max_len  

        return min_j,max_j
    
    def traceback(self,n_frag):
        # used the generated matrix and trace dictionary to find n_trace number of traces
        # the traces will be rank ordered by trace score.

        # make a score list [(score,i,j)] for every j that is the terminal fragment i.e. j+max_len > len(seq)
        score_li = []
        if not n_frag:
            for i in range(len(self.mat)):
                for j in range(len(self.seq)):
                    if self.mat[i][j] and j+self.max_len > len(self.seq):
                        score_li.append((self.mat[i][j],self.split_mat[i][j],i,j))
        else:
            i = n_frag-1
            for j in range(len(self.seq)):
                if self.mat[i][j] and j+self.max_len > len(self.seq):
                    score_li.append((self.mat[i][j],self.split_mat[i][j],i,j))
        # sort score list
        # self.sorted_score_li = sorted(score_li,key=lambda x:(x[1],-x[0]))
        self.sorted_score_li = sorted(score_li,key=lambda x:(-x[0]))
        assert self.sorted_score_li, 'No starting point for traces. Usually due to infeasible fragment length and fragment number combination'
        self.traces, self.trace_scores = [], []
        self.trace_split_scores = []

        
        # print(self.sorted_score_li)
        # make_traces
        n = 0
        for score,split_scores,i,j in self.sorted_score_li[:self.n_trace]:  
            n+=1
            trace = self.__trace(i,j)
            print('Trace',n,' breaks:',' ,'.join([str(x) for x in trace]),'; Total score:',score)
            # print([self.get_score_for_j(j) for j in trace])
            self.traces.append(trace)
            self.trace_scores.append(score)
            self.trace_split_scores.append(split_scores)            

    def __trace(self,i,j):
        js = [j]
        while i:
            key = self.make_key(i,j)
            i,j = self.unmake_key(self.trace_di[key])
            js.append(j)
        return js

    # def find_number_of_identity(self,js,ovrhg):
    #     ns = [x for x in _ovrhg = self.seq[j:j+4]]

    def write_outfile(self,filename_prifix):
        import json
        '''
        Writes text output to the filename_prifix.txt
        matrix to filename_prifix_matrix.json
        traces to filename_prifix_traces.json
        '''  

        out_name = filename_prifix+'.txt'
        mat_name = filename_prifix+'_matrix.json'
        trace_name = filename_prifix+'_trace.json'      
        # describe inputs
        outlines = []
        a = outlines.append
        a('File contains the a formated output of the generated fragments')

        with open(mat_name,'w') as fp:
            json.dump(self.mat,fp)
        a('Matrix saved as: '+mat_name)
        print('Wrote matrix file:',mat_name)

        with open(trace_name,'w') as fp:
            json.dump(self.trace_di,fp)
        print('Wrote trace file:',trace_name)

        a('Trace saved as : '+trace_name)
        a('')
        a('Input sequence')
        a('==============')
        a(self.seq)
        a('')
        a('Lenght of input: '+str(len(self.seq)))    

        # describe parameters used
        a('Parameters Used')
        a('===============')
        a('\tLower limit of fragment length: {}'.format(self.min_len))
        a('\tUpper limit of fragment length: {}'.format(self.max_len))
        a('\tStart position index          : {}'.format(self.start_site))
        a('\tEfficenfy weight              : {}'.format(self.eff_w))
        a('\tFidility weight               : {}'.format(self.fid_w))
        a('\tNumber of output traces       : {}'.format(self.n_trace))
        a("\t5' overhang                   : {}".format(self._5overhang))
        a("\t3' overhang                   : {}".format(self._3overhang))
        a('\tInitiation score              : {}'.format(self.inti_score))
        a('\tExcluded positions            : {}'.format(' ,'.join([str(x) for x in self.exclusion_list])))
        a('\tAlien overhangs               : {}'.format(' ,'.join([str(x) for x in self.alien_overhangs])))
        a('')
        a('Tracebacks(Predictions)')
        a(''.join(['-' for x in range(len(self.seq))]))

        for i in range(len(self.traces)):
            a('Traceback '+str(i+1))
            a(''.join(['+' for x in outlines[-1]]))
            a('Score (higher is better)            : {:3.2f}'.format(self.trace_scores[i]))
            a('Fragmentation positions(starts at 0): '+' ,'.join([str(x) for x in self.traces[i]]))
            # compute net efficency and fidility
            eff_tot,fid_tot = self.trace_split_scores[i]
            a('Net ligation efficency              : {:5.2f} %'.format(eff_tot*100))
            a('Net Fidility (accuracy)             : {:5.2f} %'.format(fid_tot*100))
            a('Net success rate                    : {:5.2f} %'.format(eff_tot*fid_tot*100))
            a('')
            a('Numbering 100s:   '+''.join([str(x)[-3] if (len(str(x)) >=3 and str(x)[-1] == '0' and str(x)[-2] == '0') else ' ' for x in range(0,len(self.seq))]))
            a('Numbering  10s:   '+''.join([str(x)[-2] if (len(str(x)) >=2 and str(x)[-1] == '0') else ' ' for x in range(0,len(self.seq))]))
            a('Numbering   1s:   '+''.join([str(x)[-1] for x in range(0,len(self.seq))]))
            a("Sequence      :5'-"+self.seq+"-3'")
            a('Break points  :   '+''.join('^' if x in self.traces[i][:-1] else ' ' for x in range(len(self.seq)) ))      
            # compute the fragments from trace
            fragments,effs,fids,range_li = self.get_fragments(self.traces[i])
            a('') 
            a('Fragments')
            a('---------')
            for i,fragment in enumerate(fragments):
                full_frag = self._5overhang+fragment+self._3overhang
                a('Fragment {}'.format(i+1))
                a("Range {}-{}; length: {} ; length with additions: {}".format(range_li[i][0],range_li[i][1],len(fragment),len(full_frag)))
                a("5' and 3' efficencies: {}, {}; 5' and 3' fidilities: {}, {}".format(effs[i][0],effs[i][1],fids[i][0],fids[i][1]))
                gap_spots = list(range(len(self._5overhang),len(self._5overhang)+4))+list(range((len(full_frag)-len(self._3overhang)-4),(len(full_frag)-len(self._3overhang))))
                a('Overhangs sites   :   '+''.join('m' if x in gap_spots else ' ' for x in range(len(full_frag))))
                a("Sequence          :5'-{}-3'".format(full_frag))                
                a("Reverse complement:3'-{}-5'".format(self.get_rc(full_frag,reverse=False)))
                a('Non ATCG bases    :   '+''.join('*' if x not in ['A','T','G','C'] else ' ' for x in full_frag))
                a('')
            a('')
            a(''.join(['-' for x in range(len(self.seq))]))
        with open(out_name,'w') as fp:
            fp.writelines([x+'\n' for x in outlines]) 
        print('Wrote output file:',out_name)


    def get_fragments(self,js):
        js_rev = js[::-1]
        fragments = [self.seq[:js_rev[1]+4]]
        eff,fid = self.get_score_for_j(js_rev[1])
        efficencies = [('-',self.fl(eff))] # [[Eff of forward, Eff of reverse complent]]
        fidilities = [('-',self.fl(fid))] # [[fid of forward, fid of reverse complent]]   
        range_li = [(0,js_rev[1])]    
        
        for i in range(1,len(js_rev)):
            try:
                fragments.append(self.seq[js_rev[i]:js_rev[i+1]+4])
                eff,fid = self.get_score_for_j(js_rev[i])
                _eff,_fid = self.get_score_for_j(js_rev[i+1])
                efficencies.append((self.fl(eff),self.fl(_eff)))
                fidilities.append((self.fl(fid),self.fl(_fid)))
                range_li.append((js_rev[i],js_rev[i+1]))
            except IndexError:
                break
        fragments.append(self.seq[js_rev[i]:])        
        efficencies.append((self.fl(_eff),'-'))
        fidilities.append((self.fl(_fid),'-'))
        range_li.append((js_rev[i],len(self.seq)+1))

        return fragments,efficencies,fidilities,range_li
    
    def get_score_for_j(self,j):
        ovrhg = self.seq[j:j+4]
        eff,fid = self.score_di[ovrhg] 
        return eff,fid

    def get_rc(self,DNA,reverse=True):
        self.DNA_compelemnts= {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 'N': 'N'}
        if reverse:
            return ''.join([self.DNA_compelemnts[x] for x in DNA[::-1]])
        else:
            return ''.join([self.DNA_compelemnts[x] for x in DNA])
    
    # def get_n_identity(self,str1,str2): return sum([1 for x in range(len(str1)) if str1[x] == str2[x]])
    def make_key(self,i,j): return str(i)+'_'+str(j)
    def unmake_key(self,key): return int(key.split('_')[0]),int(key.split('_')[1])
    def round_up(self,number): return int(number) + (number % 1 > 0)
    def fl(self,a): return '{:5.2f}'.format(a)

import argparse

# parsing arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description='''Program takes in a DNA sequence and fragments it for golden gate cloning. 
                                 This is done such that the overhangs of each fragment (except the first and the last) have minimum error rate
                                 and maximum efficency. The data used to calculate was published by Potapov et. al. 2018. 
                                 https://doi.org/10.1021/acssynbio.8b00333
                                 Helpfull information: 3' attachment for BsaI = GGCTACGGTCTCC, 5' attachment for BsaI = GGAGACCGTAGCC
                                 '''
                                 )
# # positional arguments 
parser.add_argument("input", help="Input file with DNA sequence in fasta format. Only the first sequence is processed",type=str)
parser.add_argument("min_len", help="Minimum length of fragment",type=int)
parser.add_argument("max_len", help="Maximum length of fragment",type=int)
parser.add_argument("output", help="Prefix of outut files. If there is a directory in the prefix and it is not present the directory will be created",type=str)
# optional argument1
parser.add_argument("-n_frag", help="Number of fragments",type=int,default=False)
parser.add_argument("-data_file", help="Supplimentarty file from Potapov et. al. 2018 to be used. They are in ./lib",type=str,default='./lib/FileS04_T4_18h_37C.csv')
parser.add_argument("-start", help="Start site. This is the position of the first fragment",type=int,default=0)
parser.add_argument("-eff_w", help="Weight of efficency metric, Power of efficenfy metric", type=float,default=1)
parser.add_argument("-fid_w", help="Weight of fidility metric, Power of fidility metric", type=float,default=1)
parser.add_argument("-ovrhg_iden", help="Maximum identity possible between two overhanges", type=int,default=2)
parser.add_argument("-n_sol", help="Number of solutions", type=int,default=5)
parser.add_argument('-exclude', nargs='*', default=[], type=int, help='list of DNA indices where overhangs cannot be made')
parser.add_argument("-add3", help="Additional sequences to be added to 3' site. Usefull for adding restriction enzyme binding site",type=str,default='')
parser.add_argument("-add5", help="Additional sequences to be added to 5' site. Usefull for adding restriction enzyme binding site",type=str,default='')
parser.add_argument('-alien_overhangs', nargs='*', default=[], type=str, help='list of overhangs that are not from input that need to be considered for overlap check')

a = parser.parse_args()

assert (a.max_len - a.min_len) >= 4, 'Confused OOGGA booga: Difference between minimum and maximum length must be greater than 4.'

try:
    seq = list(read_fasta(a.input).values())[0]
    df = Dyna_frag(seq,a.min_len,a.max_len,a.data_file,n_frag=a.n_frag,n_trace = a.n_sol, start_site = a.start, eff_w = a.eff_w,fid_w = a.fid_w, exclusion_list=a.exclude,inti_score = 1,max_overhang_identity = a.ovrhg_iden, alien_overhangs=a.alien_overhangs)
    df._5overhang = a.add3
    df._3overhang = a.add5
    df.write_outfile(a.output)
except Exception as e:
    print('Confused OOGGA booga:',e)
