# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 13:55:19 2026

@author: haapalt4
"""
from netCDF4 import Dataset
import math
import numpy as np
import time
from scipy.stats import gamma
import os
from datetime import timedelta, datetime


"""Functions for event analysis."""
def summer_flag(all_time_data):
    """This function is for creating a summer mask to allow study of the summer / growing season."""
    summer = []
    length = len(all_time_data)
    a = [timedelta(x) + datetime(1970,1,1) for x in all_time_data]
    for i in range(length):
        kk = a[i].month
        if kk < 9 and kk > 4:
            summer.append(1)
        else:
            summer.append(0)
    return a, summer

def analyze_droughts(SPI_list, threshold_start_drought, threshold_end_drought, summer_mask, summer):
    """This function determines the dry events ("droughts") occurring during the analysis period."""
    list_length = len(SPI_list) #length for loop calculation
    drought_events = [] #set lists for logging events
    summer_drought_events_dates = []
    number_of_droughts = 0
    number_of_summer_droughts = 0
    drought = False #initial drought status is false.
    
    """Maximum SPI calculations"""
    maximum_drought_SPI = min(SPI_list)
    summer_spi_previous = 0 
    for a in range(list_length):
        if summer_mask[a] == 1:
            summer_spi_current = SPI_list[a]
            if summer_spi_current < summer_spi_previous:
                   summer_spi_previous = summer_spi_current
    maximum_drought_SPI_summer = summer_spi_previous
    drought_onoff = [] #including this variable creates a time series where dry events are shown with boolean value 1 during the analysis period
    for i in range(list_length): #go through every item on the spi list. i is a date index
        if drought == False: #if there is currently no drought,
            if SPI_list[i] <= threshold_start_drought: #compare whether the current SPI list value is below the threshold
                drought_date_indices = [] #set the list for the pair of date indices, corresponding to the beginning and ending indices of the drought
                drought = True #drought status is true. required to distinguish multiday droughts
                drought_date_indices.append(i) #append beginning date index
                if i == list_length -1: #if the last item has a drought-qualifying SPI, then it is a single day drought. In an annual dataset, this means that a drought beginning on Dec 31st is counted as a single day drought occurring at the end of the year, even if the actual event continues into the new year.
                    drought_date_indices.append(i)
                    number_of_droughts += 1 #count increases
                    drought_events.append(drought_date_indices) #index pair appended
        else: # if drought is true, progress through the list. 
            if SPI_list[i] > threshold_end_drought: #check the SPI value while it is happening though. If it goes above the threshold,
                drought = False #drought is over.
                number_of_droughts += 1 #one drought event has happened.
                drought_date_indices.append(i) #append end date index. it is now in the shape of [i_begin, i_end]
                drought_events.append(drought_date_indices) #append the index pair into events
            else: #check if final index has drought. Drought gets cut short as the year ends even if it would continue in real life. This is to identify the actual number of events occurring during a year, as yearly event valuea are used later.
                if i == list_length -1:
                    drought_date_indices.append(i)
                    number_of_droughts += 1
                    drought_events.append(drought_date_indices)
        if drought == True:
            drought_onoff.append(1)
        else:
            drought_onoff.append(0)
    """Total length calculations"""
    lengths_of_droughts = [] #set list for lengths of droughts, determined by subtracting the end index from the beginning index
    lengths_of_summer_droughts = [] #subset list for droughts that occur during summer
    """Droughts"""
    for j in range(number_of_droughts): #go through the number,        
        lengths_of_droughts.append(drought_events[j][1]-drought_events[j][0]) #and append the lengths by subtracting the beginning from the end index
        #print(drought_events[j][0],drought_events[j][1])
        if summer_mask[drought_events[j][0]] > 0 or summer_mask[drought_events[j][1]] > 0: #check if date indices match with summer mask in either index!
            number_of_summer_droughts += 1 #drought occurred entirely or partially during summer
            lengths_of_summer_droughts.append(drought_events[j][1]-drought_events[j][0]) #length appended in similar fashion to three rows above
            summer_drought_events_dates.append([summer[drought_events[j][0]],summer[drought_events[j][1]]])
    """Average drought length calculations"""
    if number_of_droughts == 0: #if droughts did not occur (unlikely, but possible)
        average_drought_length = 0 #average length set to zero to avoid division by zero
    else: #droughts did occur. calculate average as usual
        average_drought_length = sum(lengths_of_droughts) / number_of_droughts
    if number_of_summer_droughts == 0:
        average_summer_drought_length = 0
    else:
        average_summer_drought_length = sum(lengths_of_summer_droughts) / number_of_summer_droughts
    
    return drought_events, number_of_droughts, average_drought_length, lengths_of_droughts, \
        number_of_summer_droughts, average_summer_drought_length, lengths_of_summer_droughts, \
        maximum_drought_SPI, maximum_drought_SPI_summer, summer_drought_events_dates, drought_onoff

def analyze_wetnesses(SPI_list, threshold_start_wetness, threshold_end_wetness, summer_mask, summer):
    """The same code, but mirrored for wet events."""
    list_length = len(SPI_list) #length for loop calculations
    wetness_events = []
    summer_wetness_events_dates = []
    number_of_wetnesses = 0
    number_of_summer_wetnesses = 0
    wetness = False #initial wetness status too
    """Maximum SPI calculations"""
    maximum_wetness_SPI = max(SPI_list)
    summer_spi_previous = 0
    for a in range(list_length):
        if summer_mask[a] == 1:
            summer_spi_current = SPI_list[a]
            if summer_spi_current > summer_spi_previous:
                   summer_spi_previous = summer_spi_current
    maximum_wetness_SPI_summer = summer_spi_previous
    wetness_onoff = []
    for i in range(list_length):
        if wetness == False:
            if SPI_list[i] >= threshold_start_wetness:
                wetness_date_indices = []
                wetness = True
                wetness_date_indices.append(i)
                if i == list_length -1:
                    wetness_date_indices.append(i)
                    number_of_wetnesses += 1
                    wetness_events.append(wetness_date_indices)
        else:
            if SPI_list[i] < threshold_end_wetness:
                wetness = False
                number_of_wetnesses += 1
                wetness_date_indices.append(i)
                wetness_events.append(wetness_date_indices)
            else:
                if i == list_length -1:
                    wetness_date_indices.append(i)
                    number_of_wetnesses += 1
                    wetness_events.append(wetness_date_indices)
        if wetness == True:
            wetness_onoff.append(1)
        else:
            wetness_onoff.append(0)
    lengths_of_wetnesses = []
    lengths_of_summer_wetnesses = []
    for k in range(number_of_wetnesses):
        lengths_of_wetnesses.append(wetness_events[k][1]-wetness_events[k][0])
        if summer_mask[wetness_events[k][0]] > 0 or summer_mask[wetness_events[k][1]] > 0:
            number_of_summer_wetnesses += 1
            lengths_of_summer_wetnesses.append(wetness_events[k][1]-wetness_events[k][0])
            summer_wetness_events_dates.append([summer[wetness_events[k][0]],summer[wetness_events[k][1]]])
    if number_of_wetnesses == 0:
        average_wetness_length = 0
    else:
        average_wetness_length = sum(lengths_of_wetnesses) / number_of_wetnesses
    if number_of_summer_wetnesses == 0:
        average_summer_wetness_length = 0
    else:
        average_summer_wetness_length = sum(lengths_of_summer_wetnesses) / number_of_summer_wetnesses
    return wetness_events, number_of_wetnesses, average_wetness_length, lengths_of_wetnesses, number_of_summer_wetnesses, average_summer_wetness_length, lengths_of_summer_wetnesses, maximum_wetness_SPI, maximum_wetness_SPI_summer, summer_wetness_events_dates, wetness_onoff

def trend_of_summer_drought_lengths(summer_drought_events_dates, lengths_of_summer_droughts, analysis_length):
    yearly_summer_drought_lengths = []
    yearly_summer_drought_number = []
    current_year = summer_drought_events_dates[0][0].year #current year taken from datetime formatting
    total_length_of_current_year_summer_droughts = 0 #no initial total length
    number_of_current_year_summer_droughts = 0 #or number
    for i in range(len(summer_drought_events_dates)): #go through the summer drought events
        year_start = summer_drought_events_dates[i][0].year #starting year of the summer drought
        year_end = summer_drought_events_dates[i][1].year #ending year of the summer drought
        if year_start == current_year or year_end == current_year: #check if the current summer drought occurs during the currently stored year
            number_of_current_year_summer_droughts += 1 #add 1 to the number, the currently drought occurred during the current year
            total_length_of_current_year_summer_droughts += lengths_of_summer_droughts[i] #add the length of this drought to the current year total length
        else: #if the summer drought did not occur during the current year...
            if number_of_current_year_summer_droughts == 0: #if there are no droughts occurring the current year,
                mean_length_of_current_year_summer_droughts = 0 #the mean length is zero - there were no summer droughts this year. Example: discretisationstep = 10, cell (109,33), year 1962
                yearly_summer_drought_lengths.append(mean_length_of_current_year_summer_droughts) #append the mean
            else: #no more summer droughts occur during the current year
                yearly_summer_drought_lengths.append(total_length_of_current_year_summer_droughts / number_of_current_year_summer_droughts) #append the mean length, calculated properly
            yearly_summer_drought_number.append(number_of_current_year_summer_droughts) #store the number of summer droughts into a variable
            total_length_of_current_year_summer_droughts = 0 #reset total length variable
            number_of_current_year_summer_droughts = 0 #reset number variable
            current_year += 1 #proceed to next year
            if year_start == current_year or year_end == current_year: #handles summer droughts that skip a year(?)
                number_of_current_year_summer_droughts += 1 #drought occurrence recorded
                total_length_of_current_year_summer_droughts += lengths_of_summer_droughts[i] #length of drought recorded
        #print(i, current_year, year_start, year_end, number_of_current_year_summer_droughts)
    """Last year appends"""
    if number_of_current_year_summer_droughts == 0: #if there are no summer droughts during the last year,
        mean_length_of_current_year_summer_droughts = 0 #mean length is zero
        yearly_summer_drought_lengths.append(mean_length_of_current_year_summer_droughts) #and it gets appended
    else: #if there are summer droughts during the last year,
        yearly_summer_drought_lengths.append(total_length_of_current_year_summer_droughts / number_of_current_year_summer_droughts) #append the properly calculated mean length
    yearly_summer_drought_number.append(number_of_current_year_summer_droughts) #and store the number of summer droughts into a variable 
    total_length_of_current_year_summer_droughts = 0 #reset for next cell
    number_of_current_year_summer_droughts = 0 #reset for next cell
    return yearly_summer_drought_lengths, yearly_summer_drought_number

def trend_of_summer_wetness_lengths(summer_wetness_events_dates, lengths_of_summer_wetnesses, analysis_length):
    yearly_summer_wetness_lengths = [] #sisältää kuivuuksien keskipituuden kesällä
    yearly_summer_wetness_number = [] #kesäkuivuuksien määrä
    current_year = summer_wetness_events_dates[0][0].year #current year taken from datetime formatting
    total_length_of_current_year_summer_wetnesses = 0 #no initial total length
    number_of_current_year_summer_wetnesses = 0 #or number
    for i in range(len(summer_wetness_events_dates)): #go through the summer drought events
        year_start = summer_wetness_events_dates[i][0].year #starting year of the summer drought
        year_end = summer_wetness_events_dates[i][1].year #ending year of the summer drought
        if year_start == current_year or year_end == current_year: #check if the current summer drought occurs during the currently stored year
            number_of_current_year_summer_wetnesses += 1 #add 1 to the number, the currently drought occurred during the current year
            total_length_of_current_year_summer_wetnesses += lengths_of_summer_wetnesses[i] #add the length of this drought to the current year total length
        else: #if the summer drought did not occur during the current year...
            if number_of_current_year_summer_wetnesses == 0: #if there are no droughts occurring the current year,
                mean_length_of_current_year_summer_wetnesses = 0 #the mean length is zero - there were no summer droughts this year. Example: discretisationstep = 10, cell (109,33), year 1962
                yearly_summer_wetness_lengths.append(mean_length_of_current_year_summer_wetnesses) #append the mean
            else: #no more summer droughts occur during the current year
                yearly_summer_wetness_lengths.append(total_length_of_current_year_summer_wetnesses / number_of_current_year_summer_wetnesses) #append the mean length, calculated properly
            yearly_summer_wetness_number.append(number_of_current_year_summer_wetnesses) #store the number of summer droughts into a variable
            total_length_of_current_year_summer_wetnesses = 0 #reset total length variable
            number_of_current_year_summer_wetnesses = 0 #reset number variable
            current_year += 1 #proceed to next year
            if year_start == current_year or year_end == current_year: #handles summer droughts that skip a year(?)
                number_of_current_year_summer_wetnesses += 1 #drought occurrence recorded
                total_length_of_current_year_summer_wetnesses += lengths_of_summer_wetnesses[i] #length of drought recorded
        #print(i, current_year, year_start, year_end, number_of_current_year_summer_droughts)
    """Last year appends"""
    if number_of_current_year_summer_wetnesses == 0: #if there are no summer droughts during the last year,
        mean_length_of_current_year_summer_wetnesses = 0 #mean length is zero
        yearly_summer_wetness_lengths.append(mean_length_of_current_year_summer_wetnesses) #and it gets appended
    else: #if there are summer droughts during the last year,
        yearly_summer_wetness_lengths.append(total_length_of_current_year_summer_wetnesses / number_of_current_year_summer_wetnesses) #append the properly calculated mean length
    yearly_summer_wetness_number.append(number_of_current_year_summer_wetnesses) #and store the number of summer droughts into a variable 
    total_length_of_current_year_summer_wetnesses = 0 #reset for next cell
    number_of_current_year_summer_wetnesses = 0 #reset for next cell
    return yearly_summer_wetness_lengths, yearly_summer_wetness_number
    
def trend_string_to_trend_number(trend_string):
    if trend_string == "increasing":
        trend_number = 3
    elif trend_string == "decreasing":
        trend_number = 1
    else:
        trend_number = 2
    return trend_number

discretisationstep = int(input("Please enter the size of the grid cell used ('1' for one square kilometer, '10' for 10 km * 10 km, '100' for 100 km * 100 km.)\n"))
accumulation_period_input = int(input("Please enter the accumulation period ('1' for SPI-1, '3' for SPI-3, '6' for SPI-6, '12' for SPI-12)\n"))
accumulation_periods_dictionary = {1:30, 3:90, 6:180, 12:360} 
accumulation_period = accumulation_periods_dictionary[accumulation_period_input] #Set the amount of days used for each accumulation period.

timer_start = time.time()
raw_data_list = os.listdir(path='raw_data')
file_dictionary = {}
for i in range(len(raw_data_list)):
    file_dictionary[i] = 'raw_data/' + raw_data_list[i]
file = Dataset(file_dictionary.get(0), mode='r')
files = [0, len(raw_data_list)-1]
analysis_length = files[1] - files[0] + 1
lons = file.variables['Lon'][:] #longitudes and latitudes to determine the spatial extent of the data
lats = file.variables['Lat'][:]
prtemp = file.variables['RRday'][:] 
file.close()

length_of_data = 0
for x in range(files[0], files[1]+1):
    file = Dataset(file_dictionary.get(x), mode='r')
    length_of_data += len(file.variables['RRday'][:])

"""Grid size in 1kmx1km"""
imax = len(lats)
jmax = len(lons)
"""Modify grid size according to discretisation step"""
imaxpick = int(imax / discretisationstep)  # Lat index
jmaxpick = int(jmax / discretisationstep)  # Lon index

spi_map = np.zeros((imaxpick,jmaxpick,length_of_data))

"""SPI constants for calculation (Lloyd-Hughes)"""
#gamma_variables = [2.515517, 0.802853, 0.010328, 1.432788, 0.189269, 0.001308]
c_0 = 2.515517
c_1 = 0.802853
c_2 = 0.010328
d_1 = 1.432788
d_2 = 0.189269
d_3 = 0.001308

temp_name = 'Time'
precipitationdatas = []
progress_milestones = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]


for i in range(imaxpick):
    for j in range(jmaxpick):
        all_precipitation_data = []
        all_time_data = []
        long_SPI_list = []
        icurrent = i*discretisationstep
        jcurrent = j*discretisationstep
        """Check whether the current cell has data via comparison against the fill value of the data array."""
        if prtemp[0,icurrent,jcurrent] != -3.4e+38: #check for non-mask values
            for n in range(files[0], files[1]+1): #Cell identified. Go through data year by year.
                file = Dataset(file_dictionary.get(n), mode='r') #lue tiedosto.
                precipitation_data = file.variables['RRday'][:,icurrent,jcurrent]
                time_data = file.variables[temp_name][:]
                all_precipitation_data = np.append(all_precipitation_data, precipitation_data)
                all_time_data = np.append(all_time_data, time_data)
            total_sum_precipitations = np.zeros(len(all_precipitation_data))
            precipitationdatas.append(all_precipitation_data)
            average_sum_precipitation = np.mean(all_precipitation_data)*accumulation_period
            for m in range(len(all_precipitation_data)):
                if m < (accumulation_period - 1):
                    total_sum_precipitations[m] = average_sum_precipitation
                else:
                    total_sum_precipitations[m] = sum(all_precipitation_data[m-(accumulation_period - 1):m+1])               
            """Input dataset created. Calculation of additional variables for SPI calculation"""
            count = len(total_sum_precipitations)
            spi_series_average = sum(total_sum_precipitations) / count
            spi_series_ln_average = math.log(spi_series_average)
            spi_series_ln_precipitations = []
            """Occurrence of zero values in the dataset."""
            count_zero = 0
            for o in range(len(total_sum_precipitations)):
                if total_sum_precipitations[o] == 0:
                    count_zero +=1
                    spi_series_ln_precipitations.append(0)
                else:
                    spi_series_ln_precipitations.append(math.log(total_sum_precipitations[o]))
            spi_series_ln_precipitations_sum = sum(spi_series_ln_precipitations)
            variable_a = spi_series_ln_average - (spi_series_ln_precipitations_sum/count)
            variable_alpha = (1/(4*variable_a))*(1+math.sqrt(1+(4/3)*variable_a))
            variable_beta = spi_series_average / variable_alpha
            """Cumulative distribution function required."""
            """cdf(x,a,loc=o,scale=1) t. scipystats"""
            #a = shape = alpha, loc = shift = location, scale = scale = beta
            """G(x)"""
            gamma_distro_function = gamma.cdf(total_sum_precipitations, variable_alpha, loc=0, scale=variable_beta)
            variable_q = count_zero / count
            """H(x) = q+(1-q)*G(x)"""
            cumulative_probability_function = variable_q + (1-variable_q)*gamma_distro_function
            h = 0
            SPI_list = []
            for h in range(len(cumulative_probability_function)):
                if cumulative_probability_function[h] > 0.5:
                    """IF 0.5 < H(x) < 1"""
                    variable_t = math.sqrt(math.log(1/((1-(cumulative_probability_function[h]))**2)))
                    SPI_value = variable_t - ((c_0 + c_1*variable_t + c_2*variable_t**2)/ \
                                              (1+d_1*variable_t+d_2*variable_t**2+d_3*variable_t**3))
                    SPI_list.append(SPI_value)
                else:
                    """IF 0 < H(x) <= 0.5"""                            
                    variable_t = math.sqrt(math.log(1/(cumulative_probability_function[h])**2))
                    SPI_value = -1* (variable_t - ((c_0 + c_1*variable_t + c_2*variable_t**2)/ \
                                                   (1+d_1*variable_t+d_2*variable_t**2+d_3*variable_t**3)))
                    SPI_list.append(SPI_value)
                long_SPI_list.append(SPI_value)
            spi_map[i,j] = long_SPI_list
            summer, summer_mask = summer_flag(all_time_data)

    progress = round(i/imaxpick,1)
    if progress in progress_milestones:
        print(int(progress*100), "%")
        progress_milestones.remove(progress)

timer_end = time.time()
duration = timer_end - timer_start
print(round(duration,1), "seconds elapsed.")


"""Saving the SPI map and key indicators into analysis_folder as csv files for later reading"""

def create_analysis_folder():
    analysis_folder = os.getcwd() + "/analysis_folder" 
    try:
        os.makedirs(analysis_folder)
    except OSError as iex:
        print(f"Creation of directories failed! {iex}")
    return analysis_folder

analysis_folder = create_analysis_folder()

spi_map_ravelled = np.ravel(spi_map) #turn 3dmap to 1d for saving
np.savetxt(analysis_folder + "/spi_map_ravelled.csv", spi_map_ravelled, delimiter=",")
np.savetxt(analysis_folder + '/key_discstep_datalengths.csv', [discretisationstep, len(raw_data_list), length_of_data], delimiter=",")

summer_array = np.array(summer)
np.savetxt(analysis_folder + '/summer.csv', summer_array, delimiter = ",", fmt='%s')
np.savetxt(analysis_folder + '/summer_mask.csv', summer_mask, delimiter = ",")





