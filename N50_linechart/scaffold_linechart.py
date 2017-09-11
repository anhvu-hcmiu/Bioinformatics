#!/usr/bin/env python3



#############
# LIBRARIES #
#############

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys



###################
# MINOR FUNCTIONS #
###################

def make_node_length_dict(scaffold_fasta):
    #Read scaffolds.fasta
    scaffold_content = open(scaffold_fasta)
    scaffold_dict = {}

    #Add length and node to scaffold_dict from scaffolds.fasta
    for line in scaffold_content:
        if line.startswith(">") == True:
            line = line[1:].strip().split("_") #Remove ">"
            node = ' '.join(line[0:2])
            length = int(line[3])
            scaffold_dict[node] = length

    return scaffold_dict



def calc_total_length(node_length_dict):
    list_all_length = [length for (node,length) in node_length_dict.items()]
    total_length = sum(list_all_length)
    return total_length



def calc_average_length(node_length_dict):
    list_all_length = [length for (node, length) in node_length_dict.items()]
    average_length = sum(list_all_length) / len(list_all_length)
    return average_length



def calc_max_length(node_length_dict):
    list_all_length = [length for (node,length) in node_length_dict.items()]
    return max(list_all_length)



def calc_min_length(node_length_dict):
    list_all_length = [length for (node, length) in node_length_dict.items()]
    return min(list_all_length)



###################
# MAJOR FUNCTIONS #
###################

################## Accumulative chart ##################

def make_accumulative_node_length_axis(node_length_dict, total_length):
    #Whole node-length dict
    whole_node_length_dict = node_length_dict
    
    #Convert node-length dictionary into descending length list of tuples
    descending_list = sorted(
        whole_node_length_dict.items(),
        key=lambda x: x[1]
    )[::-1]

    #Make length accumulative list
    n = 0
    accumulative_length_list = []
    accumulative_node_list = []
    length_list = []

    for (node,length) in descending_list:
        n = n + length
        length_list.append(length)
        accumulative_length_list.append(n)
        accumulative_node_list.append(node)
        if n >= (total_length / 2):
            break
    
    #Make a csv file containing 3 columns (otpional)
    df = pd.DataFrame()
    df['Node'] =  accumulative_node_list
    df['Length'] = length_list
    df['Accumulative length'] = accumulative_length_list
    df.to_csv("N50_table.csv", sep='\t', encoding='utf-8')
    
    return [accumulative_node_list,accumulative_length_list]



def accumulative_linechart(node_length_axis):
    fig,ax = plt.subplots(figsize=(20, 12))
    ax.set_xlabel("Node", fontsize=15)
    ax.set_ylabel("Accumulative length", fontsize=15)
    ax.set_title("N50 Line Graph", fontsize=25)
    x = node_length_axis[0]
    y = node_length_axis[1]
    sns.pointplot(x,y,color='green')
    fig.savefig("N50linchart.png")

########################################################################
########################################################################

################## Whole accumulative chart ##################

def make_node_length_axis(node_length_dict ,total_length):
    #Whole node-length dict
    whole_node_length_dict = node_length_dict
    
    #Convert node-length dictionary into descending length list of tuples
    descending_list = sorted(
        whole_node_length_dict.items(),
        key=lambda x: x[1]
    )[::-1]

    #Make length accumulative list
    n = 0
    accumulative_length_list = []
    accumulative_node_list = []
    length_list = []
    
    for (node,length) in descending_list:
        n = n + length
        length_list.append(length)
        accumulative_length_list.append(n)
        accumulative_node_list.append(node)

    '''
    #Make a csv file containing 3 columns (otpional)
    df = pd.DataFrame()
    df['Node'] =  accumulative_node_list
    df['Length'] = length_list
    df['Accumulative length'] = accumulative_length_list
    df.to_csv("whole_accumulative_table.csv", sep='\t', encoding='utf-8')
    '''
    return [accumulative_node_list,accumulative_length_list]



def whole_linechart(node_length_axis_1, node_length_axis_2):
    fig,ax = plt.subplots(figsize=(20, 12))
    ax.set_xlabel("Node", fontsize=15)
    ax.set_ylabel("Accumulative length", fontsize=15)
    ax.set_title("Scaffold Line Graph", fontsize=25)

    sns.pointplot(node_length_axis_1[0], node_length_axis_1[1],
                  color='green', label="Scaffold 1")
    sns.pointplot(node_length_axis_2[0], node_length_axis_2[1],
                  color='blue', label="Scaffold 2")

    fig.savefig("whole_accumulative_linechart.png")

########################################################################
########################################################################

def whole_pipeline(scaffold_input_1,scaffold_input_2,scaffold_output):
    result = open(scaffold_output,'w')
    
    #First step: make a dictionary of node and length
    node_length_dict_1 = make_node_length_dict(scaffold_input_1)
    node_length_dict_2 = make_node_length_dict(scaffold_input_2)
    
    #Second step: calculate the total length of scaffold 1 and 2
    result.write("The total length of {} is: ".format(scaffold_input_1))
    result.write(str(calc_total_length(node_length_dict_1)))
    result.write("\n")

    result.write("The total length of {} is: ".format(scaffold_input_2))
    result.write(str(calc_total_length(node_length_dict_2)))
    result.write("\n\n")

    #Third step: calculate the average length of scaffold 1 and 2
    result.write(str(calc_average_length(node_length_dict_1)))
    result.write("\n")

    result.write(str(calc_average_length(node_length_dict_2)))
    result.write("\n\n")

    #Fourth step: calculate the maximal length of scaffold 1 and 2
    result.write(str(calc_max_length(node_length_dict_1)))
    result.write("\n")

    result.write(str(calc_max_length(node_length_dict_2)))
    result.write("\n\n")
    
    #Fifth step: calculate the minimal length of scaffold 1 and 2
    result.write(str(calc_min_length(node_length_dict_1)))
    result.write("\n")

    result.write(str(calc_min_length(node_length_dict_2)))
    result.write("\n")
    
    #Sixth step: make axes of accumulative node and accumulative length
    axis_1 = make_node_length_axis(node_length_dict_1,
                                   calc_total_length(node_length_dict_1))
    axis_2 = make_node_length_axis(node_length_dict_2,
                                   calc_total_length(node_length_dict_2))
    
    #Final step: draw line chart for whole scaffolds
    whole_linechart(axis_1,axis_2)



if __name__ == '__main__':
    whole_pipeline(*sys.argv[1:])
