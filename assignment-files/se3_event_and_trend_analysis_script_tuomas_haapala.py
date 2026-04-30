# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 13:55:19 2026

@author: haapalt4
"""

import numpy as np
import time
import pymannkendall as mk
import os
from datetime import  datetime

"""SE2 was a rough transcription of my research code."""
"""An updated SE2 is included alongside SE3."""
"""The updates SE2 creates a map with the SPI series, and SE3 analyzes them with user set thresholds"""

"""gets"""
#map of Finland where each cell contains the list of daily SPI values based on user input accumulation period
"""outputs"""
#SPI trends based on threshold values input by the user

"""Submission to include"""
"""Python code for training a statistical model"""
#mann-kendall based trend tests over time
"""Python code to extract relevant information from model output"""
#trend statistics are recorded.
"""Python code to calculate model fit coefficients"""
#p-values.

"""Functions for event analysis."""
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

"""Trends"""                
"""Can't be calculated if there's only one event or no events."""
"""Function to account these situations and also calculate the proper trends."""
def mk_trend_maps_from_yearly_summer_events(event_value):
    if event_value == 0: #covers when yearly... = 0
        trend = 2
        slope = 0
        p_value = 1
    if type(event_value) is list and len(event_value) == 1: #covers when yearly ... = [x]
        trend = 2
        slope = 0
        p_value = 1   
    if type(event_value) is list and len(event_value) >= 2: # covers when yearly ... = [x,y,...,n]
        trend = trend_string_to_trend_number(mk.original_test(event_value).trend)
        slope = mk.original_test(event_value).slope
        p_value = mk.original_test(event_value).p
    return trend, slope, p_value

def hr_trend_maps_from_yearly_summer_events(event_value):
    if event_value == 0: #covers when yearly... = 0
        trend = 2
        slope = 0
        p_value = 1
    if type(event_value) is list and len(event_value) == 1: #covers when yearly ... = [x]
        trend = 2
        slope = 0
        p_value = 1   
    if type(event_value) is list and len(event_value) >= 2: # covers when yearly ... = [x,y,...,n]
        trend = trend_string_to_trend_number(mk.hamed_rao_modification_test(event_value).trend)
        slope = mk.hamed_rao_modification_test(event_value).slope
        p_value = mk.hamed_rao_modification_test(event_value).p
    return trend, slope, p_value

def yw_trend_maps_from_yearly_summer_events(event_value):
    if event_value == 0: #covers when yearly... = 0
        trend = 2
        slope = 0
        p_value = 1
    if type(event_value) is list and len(event_value) == 1: #covers when yearly ... = [x]
        trend = 2
        slope = 0
        p_value = 1   
    if type(event_value) is list and len(event_value) >= 2: # covers when yearly ... = [x,y,...,n]
        trend = trend_string_to_trend_number(mk.yue_wang_modification_test(event_value).trend)
        slope = mk.yue_wang_modification_test(event_value).slope
        p_value = mk.yue_wang_modification_test(event_value).p
    return trend, slope, p_value

currentdir = os.getcwd()+ '/analysis_folder/'
spi_map = np.loadtxt(currentdir + 'spi_map_ravelled.csv') #Note: an SPI map with discretisationstep 10 has a file size of ~700 MB
summer = np.loadtxt(currentdir + 'summer.csv', dtype=datetime)
summer = np.delete(summer, 1, axis = 1)
summer_mask = np.loadtxt(currentdir + 'summer_mask.csv')

def flatten(xss): #I love you stackexchange
    return [x for xs in xss for x in xs]
summer = flatten(np.ndarray.tolist(summer))
summer = [datetime.strptime(x, '%Y-%m-%d') for x in summer]

"""Reshape parameters depend on discretisation step."""
keyvalues = np.loadtxt(currentdir + 'key_discstep_datalengths.csv')
discretisationstep = int(keyvalues[0])
analysis_length_years = int(keyvalues[1])
analysis_length_days = int(keyvalues[2])

discretisationstep_reshape_dict = {1:[1147,661],10:[114,66],100:[11,6]}
spi_map = spi_map.reshape(discretisationstep_reshape_dict[discretisationstep][0], discretisationstep_reshape_dict[discretisationstep][1], analysis_length_days)
map_shape = list(spi_map.shape)
del map_shape[-1]
"""Event thresholds, unit = standard deviations. Negative depicts dryness, positive depicts wetness."""
print("Please enter the threshold values for determining SPI events of specific severity.")
threshold_start_drought = int(input("Dry event start threshold: "))
threshold_end_drought = int(input("Dry event end threshold: "))
threshold_start_wetness = int(input("Wet event start threshold: "))
threshold_end_wetness = int(input("Wet event end threshold: "))

timer_start = time.time()

"""SPI constants for calculation (Lloyd-Hughes and Saunders, 2002)"""
#gamma_variables = [2.515517, 0.802853, 0.010328, 1.432788, 0.189269, 0.001308] #c0,c1,c2,d1,d2,d2
c_0 = 2.515517
c_1 = 0.802853
c_2 = 0.010328
d_1 = 1.432788
d_2 = 0.189269
d_3 = 0.001308

progress_milestones = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]

"""Independent maps"""
"""nmb for event count, 
len for event length, 
s for summer events, 
and D and W for droughts and wetnesses."""
map_nmb_D = np.zeros(map_shape) #number of dry events
map_nmb_W = np.zeros(map_shape) #number of wet events
map_len_D = np.zeros(map_shape) #average length of dry events
map_len_W = np.zeros(map_shape) #average length of wet events
map_nmb_Ds = np.zeros(map_shape) #number of dry events during summer
map_nmb_Ws = np.zeros(map_shape) #number of wet events during summer
map_len_Ds = np.zeros(map_shape) #average length of dry events during summer
map_len_Ws = np.zeros(map_shape) #average length of wet events during summer

"""Mann-Kendall trend maps"""
"""Trend, slope and p-value (signifigance) recorded"""
map_len_Ds_trend_mk = np.zeros(map_shape)
map_nmb_Ds_trend_mk = np.zeros(map_shape) 
map_len_Ds_slope_mk = np.zeros(map_shape)
map_nmb_Ds_slope_mk = np.zeros(map_shape)
map_len_Ds_p_value_mk = np.zeros(map_shape)
map_nmb_Ds_p_value_mk = np.zeros(map_shape)
map_len_Ws_trend_mk = np.zeros(map_shape) 
map_nmb_Ws_trend_mk = np.zeros(map_shape) 
map_len_Ws_slope_mk = np.zeros(map_shape)
map_nmb_Ws_slope_mk = np.zeros(map_shape)
map_len_Ws_p_value_mk = np.zeros(map_shape)
map_nmb_Ws_p_value_mk = np.zeros(map_shape)

"""Hamed and Rao modified Mann-Kendall trend maps"""
map_len_Ds_trend_hr = np.zeros(map_shape)
map_nmb_Ds_trend_hr = np.zeros(map_shape) 
map_len_Ds_slope_hr = np.zeros(map_shape)
map_nmb_Ds_slope_hr = np.zeros(map_shape)
map_len_Ds_p_value_hr = np.zeros(map_shape)
map_nmb_Ds_p_value_hr = np.zeros(map_shape)
map_len_Ws_trend_hr = np.zeros(map_shape) 
map_nmb_Ws_trend_hr = np.zeros(map_shape) 
map_len_Ws_slope_hr = np.zeros(map_shape)
map_nmb_Ws_slope_hr = np.zeros(map_shape)
map_len_Ws_p_value_hr = np.zeros(map_shape)
map_nmb_Ws_p_value_hr = np.zeros(map_shape)

"""Yue and Wang modified Mann-Kendall trend maps"""
map_len_Ds_trend_yw = np.zeros(map_shape)
map_nmb_Ds_trend_yw = np.zeros(map_shape) 
map_len_Ds_slope_yw = np.zeros(map_shape)
map_nmb_Ds_slope_yw = np.zeros(map_shape)
map_len_Ds_p_value_yw = np.zeros(map_shape)
map_nmb_Ds_p_value_yw = np.zeros(map_shape)
map_len_Ws_trend_yw = np.zeros(map_shape) 
map_nmb_Ws_trend_yw = np.zeros(map_shape) 
map_len_Ws_slope_yw = np.zeros(map_shape)
map_nmb_Ws_slope_yw = np.zeros(map_shape)
map_len_Ws_p_value_yw = np.zeros(map_shape)
map_nmb_Ws_p_value_yw = np.zeros(map_shape)


for i in range(map_shape[0]):
    for j in range(map_shape[1]):
        """Check whether the current cell has data via comparison against the fill value of the data array."""
        if sum(spi_map[i,j,:]) != 0:
            long_SPI_list = spi_map[i,j,:]
            drought_events, number_of_droughts, average_drought_length, lengths_of_droughts, number_of_summer_droughts, average_summer_drought_length, lengths_of_summer_droughts, max_drought_spi, max_drought_spi_summer, summer_drought_events_dates, drought_onoff = analyze_droughts(long_SPI_list, threshold_start_drought, threshold_end_drought, summer_mask, summer)
            wetness_events, number_of_wetnesses, average_wetness_length, lengths_of_wetnesses, number_of_summer_wetnesses, average_summer_wetness_length, lengths_of_summer_wetnesses, max_wetness_spi, max_wetness_spi_summer, summer_wetness_events_dates, wetness_onoff = analyze_wetnesses(long_SPI_list, threshold_start_wetness, threshold_end_wetness, summer_mask, summer)
            
            """It is possible that no events or only one event occurs. Trend calculation doesn't work in these cases."""
            if len(summer_drought_events_dates) != 0: #possible that there are no summer droughts.
                yearly_summer_drought_lengths, yearly_summer_drought_number = trend_of_summer_drought_lengths(summer_drought_events_dates, lengths_of_summer_droughts, analysis_length_years)
            else:
                yearly_summer_drought_lengths = 0
                yearly_summer_drought_number = 0
            if len(summer_wetness_events_dates) != 0: #possible that there are no summer wetnesses.
                yearly_summer_wetness_lengths, yearly_summer_wetness_number = trend_of_summer_wetness_lengths(summer_wetness_events_dates, lengths_of_summer_wetnesses, analysis_length_years)
            else:
                yearly_summer_wetness_lengths = 0
                yearly_summer_wetness_number = 0
                        
            """Maps independent of trend tests"""
            map_nmb_D[i,j] += number_of_droughts
            map_nmb_W[i,j] += number_of_wetnesses
            map_len_D[i,j] += average_drought_length
            map_len_W[i,j] += average_wetness_length
            map_nmb_Ds[i,j] += number_of_summer_droughts
            map_nmb_Ws[i,j] += number_of_summer_wetnesses
            map_len_Ds[i,j] += average_summer_drought_length
            map_len_Ws[i,j] += average_summer_wetness_length
           
            """Mann-Kendall trend tests"""
            map_len_Ds_trend_mk[i,j], map_len_Ds_slope_mk[i,j], map_len_Ds_p_value_mk[i,j] = mk_trend_maps_from_yearly_summer_events(yearly_summer_drought_lengths)
            map_nmb_Ds_trend_mk[i,j], map_nmb_Ds_slope_mk[i,j], map_nmb_Ds_p_value_mk[i,j] = mk_trend_maps_from_yearly_summer_events(yearly_summer_drought_number)
            map_len_Ws_trend_mk[i,j], map_len_Ws_slope_mk[i,j], map_len_Ws_p_value_mk[i,j] = mk_trend_maps_from_yearly_summer_events(yearly_summer_wetness_lengths)
            map_nmb_Ws_trend_mk[i,j], map_nmb_Ws_slope_mk[i,j], map_nmb_Ws_p_value_mk[i,j] = mk_trend_maps_from_yearly_summer_events(yearly_summer_wetness_number)
            
            """Hamed and Rao modified Mann-Kendall trend tests"""
            map_len_Ds_trend_hr[i,j], map_len_Ds_slope_hr[i,j], map_len_Ds_p_value_hr[i,j] = hr_trend_maps_from_yearly_summer_events(yearly_summer_drought_lengths)
            map_nmb_Ds_trend_hr[i,j], map_nmb_Ds_slope_hr[i,j], map_nmb_Ds_p_value_hr[i,j] = hr_trend_maps_from_yearly_summer_events(yearly_summer_drought_number)
            map_len_Ws_trend_hr[i,j], map_len_Ws_slope_hr[i,j], map_len_Ws_p_value_hr[i,j] = hr_trend_maps_from_yearly_summer_events(yearly_summer_wetness_lengths)
            map_nmb_Ws_trend_hr[i,j], map_nmb_Ws_slope_hr[i,j], map_nmb_Ws_p_value_hr[i,j] = hr_trend_maps_from_yearly_summer_events(yearly_summer_wetness_number)
            
            """Yue and Wang modified Mann-Kendall trend tests"""
            map_len_Ds_trend_yw[i,j], map_len_Ds_slope_yw[i,j], map_len_Ds_p_value_yw[i,j] = yw_trend_maps_from_yearly_summer_events(yearly_summer_drought_lengths)
            map_nmb_Ds_trend_yw[i,j], map_nmb_Ds_slope_yw[i,j], map_nmb_Ds_p_value_yw[i,j] = yw_trend_maps_from_yearly_summer_events(yearly_summer_drought_number)
            map_len_Ws_trend_yw[i,j], map_len_Ws_slope_yw[i,j], map_len_Ws_p_value_yw[i,j] = yw_trend_maps_from_yearly_summer_events(yearly_summer_wetness_lengths)
            map_nmb_Ws_trend_yw[i,j], map_nmb_Ws_slope_yw[i,j], map_nmb_Ws_p_value_yw[i,j] = yw_trend_maps_from_yearly_summer_events(yearly_summer_wetness_number)


    progress = round(i/map_shape[0],1)
    if progress in progress_milestones:
        print(int(progress*100), "%")
        progress_milestones.remove(progress)
timer_end = time.time()
duration = timer_end - timer_start
print(round(duration,1), "seconds elapsed.")


"""Creating a map folder for output storage."""
def create_map_folder():
    map_folder = os.getcwd() + "/map_folder" 
    try:
        os.makedirs(map_folder)
    except OSError as iex:
        print(f"Creation of folder failed! {iex}")
    return map_folder

map_folder = create_map_folder()

def maps_to_csv(map_name, map_item):
    np.savetxt(map_folder + '/' + map_name, map_item, delimiter = ",", fmt='%s')
    
independent_maps = [map_nmb_D, map_nmb_W, map_len_D, map_len_W, map_nmb_Ds, map_nmb_Ws, map_len_Ds, map_len_Ws]
independent_map_names = ["nmb_D", "nmb_W", "len_D", "len_W", "nmb_Ds",  "nmb_Ws", "len_Ds", "len_Ws"]

for i in range(len(independent_map_names)):
    maps_to_csv(independent_map_names[i], independent_maps[i])
    
mk_trend_maps = [map_nmb_Ds_trend_mk, map_len_Ds_trend_mk, map_nmb_Ws_trend_mk, map_len_Ws_trend_mk, map_nmb_Ds_slope_mk, map_len_Ds_slope_mk, map_nmb_Ws_slope_mk, map_len_Ws_slope_mk]
mk_trend_map_names = ["nmb_Ds_trend_mk", "len_Ds_trend_mk", "nmb_Ws_trend_mk", "len_Ws_trend_mk", "nmb_Ds_slope_mk", "len_Ds_slope_mk", "nmb_Ws_slope_mk", "len_Ws_slope_mk"]

for i in range(len(mk_trend_map_names)):
    maps_to_csv(mk_trend_map_names[i], mk_trend_maps[i])
    
hr_trend_maps = [map_nmb_Ds_trend_hr, map_len_Ds_trend_hr, map_nmb_Ws_trend_hr, map_len_Ws_trend_hr, map_nmb_Ds_slope_hr, map_len_Ds_slope_hr, map_nmb_Ws_slope_hr, map_len_Ws_slope_hr]
hr_trend_map_names = ["nmb_Ds_trend_hr", "len_Ds_trend_hr", "nmb_Ws_trend_hr", "len_Ws_trend_hr", "nmb_Ds_slope_hr", "len_Ds_slope_hr", "nmb_Ws_slope_hr", "len_Ws_slope_hr"]

for i in range(len(hr_trend_map_names)):
    maps_to_csv(hr_trend_map_names[i], hr_trend_maps[i])
    
yw_trend_maps = [map_nmb_Ds_trend_yw, map_len_Ds_trend_yw, map_nmb_Ws_trend_yw, map_len_Ws_trend_yw, map_nmb_Ds_slope_yw, map_len_Ds_slope_yw, map_nmb_Ws_slope_yw, map_len_Ws_slope_yw]
yw_trend_map_names = ["nmb_Ds_trend_yw", "len_Ds_trend_yw", "nmb_Ws_trend_yw", "len_Ws_trend_yw", "nmb_Ds_slope_yw", "len_Ds_slope_yw", "nmb_Ws_slope_yw", "len_Ws_slope_yw"]

for i in range(len(yw_trend_map_names)):
    maps_to_csv(yw_trend_map_names[i], yw_trend_maps[i])
    
print("Maps saved into analysis folder as .csv files.")





