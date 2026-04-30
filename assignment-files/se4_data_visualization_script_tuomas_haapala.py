# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 15:39:04 2026

@author: haapalt4
"""
import numpy as np
import matplotlib.pyplot as plt
import os

"""SPI event and SPI event trend map visualization"""

def create_image_folder():
    image_folder = os.getcwd() + "/image_folder" 
    try:
        os.makedirs(image_folder)
    except OSError as iex:
        print(f"Creation of folder failed! {iex}")
    return image_folder

def create_trend_plot(index, map_data, colormap, title, caption):
    fig, ax = plt.subplots()
    map_image = plt.imshow(map_data, cmap = plt.get_cmap(colormap, np.max(map_data) - np.min(map_data) + 1))
    plt.tick_params(bottom=False, left = False, labelbottom=False, labelleft = False)
    colorbar = plt.colorbar(map_image, ticks=[0,1,2,3])
    colorbar.ax.set_yticklabels(["NaN","Negative","No trend","Positive"], fontsize=8)
    plt.title(title, loc = "Left", wrap = True)
    plt.text(0, map_height*1.15, caption, fontsize = 8.5, wrap=True)
    plt.savefig(image_folder + '/' + str(map_list[index]), dpi=600, bbox_inches="tight")
    plt.show()  
    return fig

def create_event_plot(index, map_data, colormap, title, caption):
    fig, ax = plt.subplots()
    plt.imshow(map_data, colormap)
    plt.tick_params(bottom=False, left = False, labelbottom=False, labelleft = False)
    plt.colorbar()
    plt.title(title, loc = "Left", wrap = True)
    plt.text(0, map_height*1.15, caption, fontsize = 8.5, wrap=True)
    plt.savefig(image_folder + '/' + str(map_list[index]), dpi=600, bbox_inches="tight")
    plt.show()
    return fig

map_folder = os.getcwd() + '\\map_folder\\'
image_folder = create_image_folder()
map_list = os.listdir(map_folder) #all available files - csv arrays
map_height = len((np.loadtxt(map_folder + map_list[0], delimiter = ",")))

"""Length of dry events"""
len_Ds = create_event_plot(1,np.loadtxt(map_folder + map_list[1], delimiter=','), "coolwarm", "Length of summer dry events", "The average length of dry events [days] occurring during summer that surpass set SPI thresholds in the analysis period.")
len_Ds_trend_mk = create_trend_plot(6,np.loadtxt(map_folder + map_list[6], delimiter=','), "coolwarm", "Length of summer dry events,\n(Mann-Kendall)", "Trend of average dry event length modelled with the Mann-Kendall trend test.")
len_Ds_trend_hr = create_trend_plot(5,np.loadtxt(map_folder + map_list[5], delimiter=','), "coolwarm", "Length of summer dry events,\n(Hamed and Rao)", "Trend of average dry event length modelled with the Hamed and Rao modified Mann-Kendall trend test.")
len_Ds_trend_yw = create_trend_plot(7,np.loadtxt(map_folder + map_list[7], delimiter=','), "coolwarm", "Length of summer dry events,\n(Yue and Wang)", "Trend of average dry event length modelled with the Yue and Wang modified Mann-Kendall trend test.")

"""Number of dry events"""
nmb_Ds = create_event_plot(17,np.loadtxt(map_folder + map_list[17], delimiter=','), "coolwarm", "Number of summer dry events", "The average number of dry events occurring during summer that surpass set SPI thresholds in the analysis period.")
nmb_Ds_trend_mk = create_trend_plot(22,np.loadtxt(map_folder + map_list[22], delimiter=','), "coolwarm", "Number of summer dry events,\n(Mann-Kendall)", "Trend of average dry event number modelled with the Mann-Kendall trend test.")
nmb_Ds_trend_hr = create_trend_plot(23,np.loadtxt(map_folder + map_list[23], delimiter=','), "coolwarm", "Number of summer dry events,\n(Hamed and Rao)", "Trend of average dry event number modelled with the Hamed and Rao modified Mann-Kendall trend test.")
nmb_Ds_trend_yw = create_trend_plot(21,np.loadtxt(map_folder + map_list[21], delimiter=','), "coolwarm", "Number of summer dry events,\n(Yue and Wang)", "Trend of average dry event number modelled with the Yue and Wang modified Mann-Kendall trend test.")

"""Length of wet events"""
len_Ws = create_event_plot(9,np.loadtxt(map_folder + map_list[9], delimiter=','), "coolwarm_r", "Length of summer wet events", "The average length of wet events [days] occurring during summer that surpass set SPI thresholds in the analysis period.")
len_Ws_trend_mk = create_trend_plot(14,np.loadtxt(map_folder + map_list[14], delimiter=','), "coolwarm_r", "Length of summer wet events,\n(Mann-Kendall)", "Trend of average wet event length modelled with the Mann-Kendall trend test.")
len_Ws_trend_hr = create_trend_plot(13,np.loadtxt(map_folder + map_list[13], delimiter=','), "coolwarm_r", "Length of summer wet events,\n(Hamed and Rao)", "Trend of average wet event length modelled with the Hamed and Rao modified Mann-Kendall trend test.")
len_Ws_trend_yw = create_trend_plot(15,np.loadtxt(map_folder + map_list[15], delimiter=','), "coolwarm_r", "Length of summer wet events,\n(Yue and Wang)", "Trend of average wet event length modelled with the Yue and Wang modified Mann-Kendall trend test.")

"""Number of wet events"""
nmb_Ws = create_event_plot(25,np.loadtxt(map_folder + map_list[25], delimiter=','), "coolwarm_r", "Number of summer wet events", "The average number of wet events occurring during summer that surpass set SPI thresholds in the analysis period.")
nmb_Ws_trend_mk = create_trend_plot(30,np.loadtxt(map_folder + map_list[30], delimiter=','), "coolwarm_r", "Number of summer wet events,\n(Mann-Kendall)", "Trend of average wet event number modelled with the Mann-Kendall trend test.")
nmb_Ws_trend_hr = create_trend_plot(29,np.loadtxt(map_folder + map_list[29], delimiter=','), "coolwarm_r", "Number of summer wet events,\n(Hamed and Rao)", "Trend of average wet event number modelled with the Hamed and Rao modified Mann-Kendall trend test.")
nmb_Ws_trend_yw = create_trend_plot(31,np.loadtxt(map_folder + map_list[31], delimiter=','), "coolwarm_r", "Number of summer wet events,\n(Yue and Wang)", "Trend of average wet event number modelled with the Yue and Wang modified Mann-Kendall trend test.")

 