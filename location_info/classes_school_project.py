class Health_Status(object):
    # enumerates all possible health statuses
    
    susceptible=0
    exposed=1
    infectious=2
    immune=3
    confined=4

class Individual(object):
    
    def __init__(self, id, health_status=None, shedding=None, breathing=None):
        # id is the individual's id
        # shedding in shedded quanta per hour
        # breathing in liter per minute; stored as m3 per hour
        
        L_TO_M3 = 1000.0
        MIN_TO_H = 60.0
        RESP_MIN_VOL_ADULT_L = 8.0
        DEFAULT_HEALTH = Health_Status.susceptible
        DEFAULT_BREATHING = (RESP_MIN_VOL_ADULT_L * MIN_TO_H) / L_TO_M3
        DEFAULT_SHEDDING = 0.0
        
        self.id = id
        
        if not health_status:
            self.health_status = DEFAULT_HEALTH
        else:
            self.health_status = health_status
        
        if not shedding:
            self.shedding = DEFAULT_SHEDDING
        else:
            self.shedding = shedding
            
        if not breathing:
            self.breathing = DEFAULT_BREATHING
        else:
            self.breathing = (breathing * MIN_TO_H) / L_TO_M3 

class Room(object):
    
    def __init__(self, id, volume, air_change_rate=None):
        # id is the room's id
        # volume is the room's volume in m3
        # air_change_rate is air changes per hour
        
        DEFAULT_AIR_CHANGE_RATE = 0.0
        
        self.id = id
        self.volume = volume
        
        if not air_change_rate:
            self.air_change_rate = DEFAULT_AIR_CHANGE_RATE
        else:
            self.air_change_rate = air_change_rate

class Occupancy(object):
    
    def __init__(self, timestamp, individual, room):
        # timestamp stands for the moment when the occupancy occurred (integer)
        # individual is an Individual object
        # room is a Room object
        
        self.timestamp = timestamp
        self.individual = individual
        self.room = room

class Sequence(object):
    
    def __init__(self):
        self.room_dict = {}
        self.individual_dict = {}
        self.t_rooms = {}
        self.t_individuals = {}
        self.concentrations = {}
        self.exposures = {}
    
    def add_occupancy(self, occupancy):
        t = occupancy.timestamp
        i = occupancy.individual
        r = occupancy.room
        
        if not (r.id in self.room_dict):
            self.room_dict[r.id] = r
        
        if not (i.id in self.individual_dict):
            self.individual_dict[i.id] = i
        
        if not (t in self.t_rooms):
            self.t_rooms[t] = {}
            self.t_rooms[t][r] = [occupancy]
        else:
            if not (r in self.t_rooms[t]):
                self.t_rooms[t][r] = [occupancy]
            else:
                if not (occupancy in self.t_rooms[t][r]):
                    self.t_rooms[t][r].append(occupancy)
        
        if not (t in self.t_individuals):
            self.t_individuals[t] = {}
            self.t_individuals[t][i] = occupancy
        else:
            if not (i in self.t_individuals[t]):
                self.t_individuals[t][i] = occupancy
            else:
                print "Individual",i,"at time",t,"posses already an occupancy."
                self.t_individuals[t][i] = occupancy
    
    def complete_room_dict(self, room):
        # room is a Room object
        
        if not room.id in self.room_dict:
            self.room_dict[room.id] = room
            message = "Room "+str(room.id)+" is never occupied"
            return message
        else:
            return None
    
    def complete_individual_dict(self, individual):
        # individual is an Individual object
        
        if not individual.id in self.individual_dict:
            self.individual_dict[individual.id] = individual
            message = "Individual "+str(individual.id)+ " is never present."
            return message
        else:
            return None
    
    def init_concentration(self, room, initial_concentration=None):
        # room is a Room object
        # initial concentration in quanta per m3
        
        DEFAULT_CONCENTRATION=0.0
        T0 = 0
        r = room.id
        
        if not initial_concentration:
            initial_concentration = DEFAULT_CONCENTRATION
        
        if not (r in self.concentrations):
            self.concentrations[r] = {}
            self.concentrations[r][T0] = initial_concentration
        else:
            print "Concentration timeline for room",r,"already exists."
            self.concentrations[r] = {}
            self.concentrations[r][T0] = initial_concentration
    
    def init_exposure(self, individual):
        # individual is an Individual object
        
        T0 = 0
        i = individual.id
        
        if not (i in self.exposures):
            self.exposures[i] = {}
            self.exposures[i][T0] = 0.0
        else:
            print "Exposure timeline for individual",i,"already exists."
            self.exposures[i] = {}
            self.exposures[i][T0] = 0.0
    
    def iterate_concentration(self, cur_t, time_stepwidth):
        # cur_t last calculated timestamp
        # time_stepwidth in minutes
        
        H_TO_MIN = (1.0/60.0)
        next_t = cur_t + 1
        
        for r in self.concentrations:
            room = self.room_dict[r]
            inflow = 0.0
            
            if next_t in self.t_rooms:
                if room in self.t_rooms[next_t]:
                    for o in self.t_rooms[next_t][room]:
                        if o.individual.shedding > 0.0:
                            inflow += (o.individual.shedding * H_TO_MIN *
                                      time_stepwidth)
            
            outflow = (self.concentrations[r][cur_t] *
                       self.room_dict[r].volume *
                       self.room_dict[r].air_change_rate *
                       H_TO_MIN * time_stepwidth)
            
            self.concentrations[r][next_t] = (self.concentrations[r][cur_t]
                                     + (inflow / self.room_dict[r].volume)
                                     - (outflow / self.room_dict[r].volume))
            
            if self.concentrations[r][next_t] < 1.0E-10:
                self.concentrations[r][next_t] = 0.0
    
    def iterate_exposure(self, cur_t, time_stepwidth):
        # cur_t last calculated timestamp
        # time_stepwidth in minutes
        
        H_TO_MIN = (1.0/60.0)
        next_t = cur_t + 1
        
        for i in self.exposures:
            individual = self.individual_dict[i]
            
            if not next_t in self.t_individuals:
                self.exposures[i][next_t] = 0.0
            elif not individual in self.t_individuals[next_t]:
                self.exposures[i][next_t] = 0.0
            else:
                o = self.t_individuals[next_t][individual]
                concentration = self.concentrations[o.room.id][next_t]
                self.exposures[i][next_t] = (concentration *
                                             self.individual_dict[i].breathing *
                                             H_TO_MIN * time_stepwidth)
    
    def exposure(self, individual, timestep):
        # individual is an Individual object
        # timestep is the timestep for which the instantaneous exposure shall
        #          be reported
        
        i = individual.id
        if timestep in self.exposures[i]:
            return self.exposures[i][timestep]
        else:
            return None
    
    def concentration(self, room, timestep):
        # room is a Room object
        # timestep is the timestep for which the instantaneous exposure shall
        #          be reported
        
        r = room.id
        if timestep in self.concentrations[r]:
            return self.concentrations[r][timestep]
        else:
            return None
    
    def tot_exposure(self, individual):
        # individual is an Individual object
        
        i = individual.id
        
        return sum(self.exposures[i].values())
    
    def peak_concentration(self, room):
        # room is a Room object
        
        r = room.id
        
        return max(self.concentrations[r].values())
    
    def exposure_caused(self, room):
        # room is a Room object
        # Returns the total exposure caused by Room on all Individual
        
        r = room
        tot_exposure = 0.0
        
        for t in self.t_rooms:
            if r in self.t_rooms[t]:
                for o in self.t_rooms[t][r]:
                    i = o.individual.id
                    tot_exposure += self.exposures[i][t]
        
        return tot_exposure

    def exposure_caused_to_individual(self, room, individual):
        # room is a Room object
        # individual is an Individual object
        # Returns the exposure caused by Room on Individual
        
        r = room
        i = individual
        tot_exposure = 0.0
        
        for t in self.t_rooms:
            if r in self.t_rooms[t]:
                for o in self.t_rooms[t][r]:
                    if i == o.individual:
                        tot_exposure += self.exposures[i.id][t]
        
        return tot_exposure
    
    def people_in_room(self, room):
        ####
        
        r = room
        individual_set = set()
        timestamps = 0
        
        for t in self.t_rooms:
            if r in self.t_rooms[t]:
                for o in self.t_rooms[t][r]:
                    individual_set.add(o.individual.id)
                    timestamps += 1
        return len(individual_set), timestamps
    
    def sum_of_timestamps(self, room):
        ####
        
        r = room
        t_cnt = 0
        
        for t in self.t_rooms:
            if r in self.t_rooms[t]:
                t_cnt += len(self.t_rooms[t][r])
        return t_cnt
    
    def max_time(self):
        t_max_i = max(self.t_individuals.keys())
        t_max_r = max(self.t_rooms.keys())
        
        if t_max_i == t_max_r:
            return t_max_i
        else:
            print "Error - different maximal timestamps for rooms and indivs"
            raise ValueError