# Data structure of aerosol exposure files

## i_to_i_ACR$$_£££££_schedule.csv

- $$: 05 means 0.5 air changes per hour (poor ventilation); 30 means 3.0 air changes (good ventilation), etc.

- £££££: simulated  exposure changes if parts of the building are better ventilated; this describes how these parts were identified and what the percentage of the indoor air is that is better ventilated. Examples:

  - opt05: 5% of indoor air has 3.0 ACR instead of 0.5 such that the effect on exposure is max
  - sch25: 25% of indoor air is improved, choice of rooms was based on accumulated occupancy according to the official school roster
  - siz15: 15% of indoor air improved, acc. to school roster inversely weighted by room size
   
- content structure:

  - infector id: id of individual who is assumed to be the infector
  - receiver id: id of individual who might or might not be exposed to quanta from infector
  - exposure: accumulated quanta exposure of receiver due to infector 


## i_to_r_ACR$$_£££££_schedule.csv

filename structure as above

- content structure: infector id, room id, exposure

	- infector id: id of individual who is assumed to be the infector
	- room id: as label indicates
	- exposure: accumulated quanta exposure that infector caused in that room
