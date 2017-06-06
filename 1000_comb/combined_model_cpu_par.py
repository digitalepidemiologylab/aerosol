"""
Created on Mon Jan 11 16:21:26 2016

@author: timo ---> main simulation
@author: lazzari ---> loop over the index case & parallisation :)
"""

import networkx
import random
import math
import sys
import csv

#preparing for parallization
import multiprocessing as mp

print "we are using, max N thread on cpu=",mp.cpu_count()

#%% ####################
# Global variables #
####################

CONTACT_FILE = sys.argv[1] # school file!
  
II_FILE = sys.argv[2]  # "./analyses/i_to_i_ACR05.csv"
IR_FILE = sys.argv[3]  # "./analyses/i_to_r_ACR05.csv"
#INDEX_CASE = int(sys.argv[4])
SHEDDING_FACTOR   = 1.0 # float(sys.argv[5]) # shifting all quanta geenration rate 
ATTENDANCE_FACTOR = 0.25 #float(sys.argv[6])    # reducing the interaction after you get sick

FACTOR = SHEDDING_FACTOR * ATTENDANCE_FACTOR

RUNS = 1000 #int(sys.argv[7])

OUTPUT_SIGN = str(sys.argv[4])

TIMESEQUENCE = [True,False,
                True,False,
                True,False,
                True,False,
                True,False,
                False,False,
                False,False]

TIME_TO_TIMESTEP = 2.0

CUR_STEP = 1

ID1 = 0 # 1st column in csv-file is the infector ID
ID2 = 1 # 2nd column is the ID of the exposed individual
ID3 = 2 # 3rd column gives the weight

EXPOSURE_DATA = 2 # 3rd column contains the exposure data: i_to_i file!

SEED_SIZE = 1

INFECTIOUS_NOTCONFINED = 1

WEIGHT_FACTOR = 0.25

ALPHA = 0.5

cnt_contact = 0
cnt_aerosol = 0
cnt_both    = 0

#%% #############
# Functions #
#############

def read_pairwise_and_exposure(contact_net,contact_expo ):

    # filename points to the file that contains the following data:
    #          Infector ID, Exposed ID, Exposure [quanta]
    
    global ID1, ID2, ID3, EXPOSURE_DATA
    global WEIGHT_FACTOR
    
    G = networkx.DiGraph()
    
    input_file = open(contact_net,'r')
     
    for line in input_file:
        
        infector = int(line.split()[ID1])
        exposed = int(line.split()[ID2])
        weight_val = int(line.split()[ID3])
        
        if not infector in G.nodes():
            G.add_node(infector)
            G.node[infector]['status_time'] = None
        
        if not exposed in G.nodes():
            G.add_node(exposed)
            G.node[exposed]['status_time'] = None
        
        if weight_val > 0:
            G.add_edge(infector, exposed, duration=(float(weight_val)*
                                                  WEIGHT_FACTOR))
     
    # reading the exposure_file
    input_file = csv.reader(open(contact_expo,'rb'))
    
    for line in input_file:
        
        infector = int(line[ID1])
        exposed = int(line[ID2])
        exposure = float(line[EXPOSURE_DATA])
        
        if not infector in G.nodes():
            G.add_node(infector)
            G.node[infector]['status_time'] = None
        
        if not exposed in G.nodes():
            G.add_node(exposed)
            G.node[exposed]['status_time'] = None
        
        if exposure > 0.0:
            G.add_edge(infector, exposed, quanta=exposure)
    
    return G

#not necessary for the inpfection process
def read_ir_exposure_data(room_expo):
    # filename points to the file that contains the following data:
    #          Infector ID, Exposed ID, Exposure [quanta]
    
    global ID1, ID2, EXPOSURE_DATA
    
    G = networkx.DiGraph()
    nodes = []
    rooms = []
    
    input_file = csv.reader(open(room_expo,'rb'))
    
    for line in input_file:
        
        infector = int(line[ID1])
        room = int(line[ID2])
        exposure = float(line[EXPOSURE_DATA])
        
        if not infector in nodes:
            G.add_node(infector)
            nodes.append(infector)
        
        if not room in rooms:
            G.add_node(room)
            rooms.append(room)
        
        if exposure > 0.0:
            G.add_edge(infector, room, weight=exposure)
    
    return G, set(nodes), set(rooms)

def transform_exposure_data(G_iie, G_ire, factor):
    
    for i in G_iie:
        for j in G_iie[i]:
            if ('quanta' in G_iie[i][j]): # this is gian's line to solve bug!
                G_iie[i][j]['quanta'] = G_iie[i][j]['quanta'] * factor #FACTOR = SHEDDING_FACTOR * ATTENDANCE_FACTOR
        for r in G_ire[i]:
            if ('weight' in G_ire[i][r]): # this is gian's line to solve bug!
                G_ire[i][r]['weight'] = G_ire[i][r]['weight'] * factor

def infect_seed(seed_size,
                G_iie,
                index=False):
    
    sus = set()
    exp = set()
    inf = set()
    con = set()
    rec = set()
    
    individuals = set(G_iie.nodes())
    if index == False:
        infectors = set(random.sample(individuals,seed_size))
    else:
        infectors = set([index])
    
    for e in infectors:
        G_iie.node[e]['status_time'] = exposure_period()
        exp.add(e)
    
    sus = individuals - inf
    
    return sus, exp, inf, con, rec

def inf_prob_aerosol(exposure):
    # exposure is the amount of quanta inhaled
    return 1 - math.exp(-exposure) # eq4  

def inf_prob_contacts(timestamps): # prob of infection after time: timestamps - cmp. PNAS
    return 1 - ((1 - 0.003)**timestamps)  # timestamps= contact duration


def exposure_by_room(G_ire,
                     inf,
                     rooms,
                     schooltime=True):
    
    r_exp = {}
    
    for r in rooms:
        r_exp[r] = 0.0
        
    if schooltime == True:
        for i in inf:
            for r in G_ire[i]:
                r_exp[r] += G_ire[i][r]['weight']  
    
    return r_exp

def iteration(G_iie,
              sus,
              exp,
              inf,
              con,
              rec,
              schooltime=True):
    
    add = {'sus':set(),
           'exp':set(),
           'inf':set(),
           'con':set(),
           'rec':set()}
    
    remove = {'sus':set(),
              'exp':set(),
              'inf':set(),
              'con':set(),
              'rec':set()}
    
    add, remove = c_to_r(G_iie,
                         con,
                         add,
                         remove)
    
    add, remove = i_to_r(G_iie,
                         inf,
                         add,
                         remove)
    
    add, remove = i_to_c(G_iie,
                         inf,
                         add,
                         remove)
    
    add, remove = e_to_i(G_iie,
                         exp,
                         add,
                         remove)
    
    if schooltime==True:
        add, remove = s_to_e(G_iie,
                             sus,
                             inf,
                             add,
                             remove)
    
    sus, exp, inf, con, rec = update_healthlists(add,
                                                 remove,
                                                 sus,
                                                 exp,
                                                 inf,
                                                 con,
                                                 rec)
    
    return sus, exp, inf, con, rec

def update_healthlists(add,
                       remove,
                       sus,
                       exp,
                       inf,
                       con,
                       rec):
    
    sus = sus - remove['sus']
    exp = exp - remove['exp']
    inf = inf - remove['inf']
    con = con - remove['con']
    rec = rec - remove['rec']
    
    sus = sus | add['sus']
    exp = exp | add['exp']
    inf = inf | add['inf']
    con = con | add['con']
    rec = rec | add['rec']
    
    return sus, exp, inf, con, rec

def c_to_r(G_iie,
           con,
           add,
           remove):
    
    for c in con:
        if G_iie.node[c]['status_time'] == 1:
            add['rec'].add(c)
            remove['con'].add(c)
        elif G_iie.node[c]['status_time'] > 1:
            G_iie.node[c]['status_time'] -= 1
        else:
            print "Error"
    
    return add, remove

def i_to_r(G_iie,
           inf,
           add,
           remove):
    
    return add, remove

def i_to_c(G_iie,
           inf,
           add,
           remove):
    
    for i in inf:
        if G_iie.node[i]['status_time'] == 1:
            remove['inf'].add(i)
            con_period = confinement_period()
            if con_period == 0:
                add['rec'].add(i)
            else:
                add['con'].add(i)
                G_iie.node[i]['status_time'] = con_period
        elif G_iie.node[i]['status_time'] > 1:
            G_iie.node[i]['status_time'] -= 1
        else:
            print "Error"
    
    return add, remove

def confinement_period(): ### Correct?
    
    for t in range(0,12):
        p = 1.0 - (0.95 ** float(t+1))
        if random.random() <= p:
            break
    
    return t

def e_to_i(G_iie,
           exp,
           add,
           remove):
    
    global INFECTIOUS_NOTCONFINED
    
    for e in exp:
        if G_iie.node[e]['status_time'] == 1:
            add['inf'].add(e)
            remove['exp'].add(e)
            G_iie.node[e]['status_time'] == INFECTIOUS_NOTCONFINED
        elif G_iie.node[e]['status_time'] > 1:
            G_iie.node[e]['status_time'] -= 1
        else:
            print "Error"
    
    return add, remove

def s_to_e(G_iie,
           sus,
           inf,
           add,
           remove):

# counting infected by different routes
    global cnt_contact
    global cnt_aerosol
    global cnt_both
    global ALPHA

    for s in sus:
        
        exposure = 0.0

        already_infect = False

        for i in inf:
            if s in G_iie[i]:
                if ('quanta' in G_iie[i][s]): # this is gian's line to solve bug!
                    exposure += G_iie[i][s]['quanta']
        infection_probability = inf_prob_aerosol(exposure) * ALPHA

        if random.random() < infection_probability:
            cnt_aerosol += 1
            already_infect = True
            remove['sus'].add(s)
            add['exp'].add(s)
            G_iie.node[s]['status_time'] = exposure_period()

        timestamps = 0.0

        for i in inf:
            if s in G_iie.neighbors(i):
                if ('duration' in G_iie[i][s]): # this is gian's line to solve bug!
                    timestamps += G_iie[i][s]['duration']

        infection_probability = inf_prob_contacts(timestamps) * (1-ALPHA)
        if random.random() < infection_probability:
            if already_infect == False:
                cnt_contact += 1
                remove['sus'].add(s)
                add['exp'].add(s)
                G_iie.node[s]['status_time'] = exposure_period()
            else:
                cnt_both += 1
                cnt_aerosol -= 1

    return add, remove


def exposure_period():
    
    SCALE_PAR = 1.10
    SHAPE_PAR = 2.21
    OFFSET = 0.5
    
    t = int(round((random.weibullvariate(SCALE_PAR, SHAPE_PAR) + OFFSET)
                  * TIME_TO_TIMESTEP))
    
    return t


#%% ##############
# INITIALIZE #
##############

G_iie = read_pairwise_and_exposure(CONTACT_FILE, II_FILE)
print "Data of",len(G_iie.nodes())," nodes read in."
G_ire, individuals, rooms = read_ir_exposure_data(IR_FILE)
transform_exposure_data(G_iie, G_ire, FACTOR)


#%% ##############
# SIMULATION on parallel pool of worker :)  --> go to the end of this section! 
##############


def par_sim(node):

    INDEX_CASE = node  # set the index case 
    #print "simulation @ index case-",INDEX_CASE," started"

    #    get global vars
    global RUNS, SHEDDING_FACTOR, ATTENDANCE_FACTOR, OUTPUT_SIGN
    global G_iie,G_ire

    #    set the 'local**' vars - ** = for each index_case 
    i_tot = []
    time_to_peak = []
    max_time = []
    room_exp = []
    indiv_ill = []
    
    for run in range(1,RUNS+1):
        
#        print "---> RUN", run
        
        t_cnt = 0
        i_max = (0, 0)
        room_exposures = {}
        for r in rooms:
            room_exposures[r] = 0.0
        
        sus, exp, inf, con, rec = infect_seed(1,
                                              G_iie,
                                              INDEX_CASE)
        
        for week in range(0,10):
            
            if len(exp)+len(inf) == 0:
                break
            
            if week == 0:
                seq1 = random.randint(0,13)
            else:
                seq1 = 0
            seq2 = 14
            
            for t in TIMESEQUENCE[seq1:seq2]:
                
                t_cnt += 1
                
                if len(inf) > i_max[0]:
                    i_max = (len(inf), t_cnt)
                
                r_exp = exposure_by_room(G_ire,
                                         inf,
                                         rooms,
                                         t)
                
                for r in rooms:
                    room_exposures[r] += r_exp[r]
                
                sus, exp, inf, con, rec = iteration(G_iie,
                                                    sus,
                                                    exp,
                                                    inf,
                                                    con,
                                                    rec,
                                                    t)
                
                if len(exp)+len(inf) == 0:
                    break
                
        time_to_peak.append(i_max[1])
        max_time.append(t_cnt)
        i_tot.append(len(con)+len(rec))
        room_exp.append(room_exposures)
        indiv_ill.append(con | rec)

    #################
    # STORE RESULTS #
    #################
    
    # filename = './scratch/res_sgn%0s_ind%0i_fac%0i_fac%0i_epi.csv' %(OUTPUT_SIGN, # timo's code
    filename = './res_sgn%0s_ind%0i_fac%0i_fac%0i_epi.csv' %(OUTPUT_SIGN,
                                                         INDEX_CASE,
                                                         int(SHEDDING_FACTOR*100),
                                                         int(ATTENDANCE_FACTOR*100))
    o = csv.writer(open(filename, 'wb'))
    for i in range(0, len(i_tot)):
        o.writerow([INDEX_CASE,
                    SHEDDING_FACTOR,
                    ATTENDANCE_FACTOR,
                    i_tot[i],
                    time_to_peak[i],
                    max_time[i]])
    
    # filename = './scratch/res_sgn%0s_ind%0i_fac%0i_fac%0i_rooms.csv' %(OUTPUT_SIGN, # timo s code
#    filename = './res_sgn%0s_ind%0i_fac%0i_fac%0i_rooms.csv' %(OUTPUT_SIGN, 
#                                                         INDEX_CASE,
#                                                         int(SHEDDING_FACTOR*100),
#                                                         int(ATTENDANCE_FACTOR*100))
#    
#    # amount of inhaled quanta (sum of exposure for all students) for each room (cols),
#    # for each run (rows), at the end of the run!
#    
#    o = csv.writer(open(filename, 'wb'))
#    room_list = sorted(rooms)
#    o.writerow(room_list)
#    for line in room_exp:
#        output_line = []
#        for r in room_list:
#            output_line.append('%2.10f'%round(line[r],10))
#        o.writerow(output_line)
#    
#    # filename = './scratch/res_sgn%0s_ind%0i_fac%0i_fac%0i_indivs.csv' %(OUTPUT_SIGN,# timo s code
#    filename = './res_sgn%0s_ind%0i_fac%0i_fac%0i_indivs.csv' %(OUTPUT_SIGN,
#                                                         INDEX_CASE,
#                                                         int(SHEDDING_FACTOR*100),
#                                                         int(ATTENDANCE_FACTOR*100))
#    
#    # binary values for each individuals (cols), for each run (rows)
#    o = csv.writer(open(filename, 'wb'))
#    indiv_list = sorted(individuals)
#    o.writerow(indiv_list)
#    for line in indiv_ill:
#        output_line = []
#        for i in indiv_list:
#            if i in line:
#                output_line.append(1)
#            else:
#                output_line.append(0)
#        o.writerow(output_line)

#%% ACUTAL RUN OF SIMULATION :)

import multiprocessing

#pool_size  = 4
pool = multiprocessing.Pool()

pool.map(par_sim, G_iie.nodes())

#%%  
