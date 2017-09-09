import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys


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
    list_all_length = [length for (node,length) in node_length_dict.items()]
    average_length = sum(list_all_length)/len(list_all_length)
    return average_length

def calc_max_length(node_length_dict):
    list_all_length = [length for (node,length) in node_length_dict.items()]
    return max(list_all_length)

def calc_min_length(node_length_dict):
    list_all_length = [length for (node,length) in node_length_dict.items()]
    return min(list_all_length)

def make_accumulative_node_length_axis(node_length_dict,total_length):
    #Whole node-length dict
    whole_node_length_dict = node_length_dict
    
    #Convert node-length dictionary into descending length list of tuples
    descending_list = sorted(whole_node_length_dict.items(), key=lambda x: x[1])[::-1]

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
    
def linechart(node_length_axis):
    fig,ax = plt.subplots(figsize=(20,12))
    ax.set_xlabel("Node",fontsize=15)
    ax.set_ylabel("Accumulative length",fontsize=15)
    ax.set_title("N50 Line Graph",fontsize=25)
    x = node_length_axis[0]
    y = node_length_axis[1]
    sns.pointplot(x,y,color='green')
    fig.savefig("N50linchart.png")

def whole_pipeline(scaffold_input,scaffold_output):
    result = open(scaffold_output,'w')
    
    #First step: make a dictionary of node and length
    node_length_dict = make_node_length_dict(scaffold_input)
    
    #Second step: calculate the total length of scaffold
    result.write("The total length of {} is: ".format(scaffold_input))
    result.write(str(calc_total_length(node_length_dict)))
    result.write("\n")

    #Third step: calculate the average length of scaffold
    result.write(str(calc_average_length(node_length_dict)))
    result.write("\n")
    
    #Fourth step: calculate the maximal length of scaffold
    result.write(str(calc_max_length(node_length_dict)))
    result.write("\n")
    
    #Fifth step: calculate the minimal length of scaffold
    result.write(str(calc_min_length(node_length_dict)))
    result.write("\n")
    
    #Sixth step: make axes of accumulative node and accumulative length
    axes = make_accumulative_node_length_axis(node_length_dict,calc_total_length(node_length_dict))
    
    #Final step: draw line chart for N50
    linechart(axes)
    
whole_pipeline(*sys.argv[1:])
