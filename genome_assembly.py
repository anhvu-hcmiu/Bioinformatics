import json
import os
import urllib.request

def csv_to_json_converter(path_input,path_output):
    csv_input = open(path_input)
    json_output = []
    header = next(csv_input).strip().split(",") #Name,Link,Status
    for line in csv_input:
        line = line.strip().split(",")
        json_output.append(dict(list(zip(header,line))))
    
    with open(path_output, 'w') as result:
        result.write(json.dumps(json_output))

def sequence_download(path_input,download_path,path_output):
    json_input = json.load(open(path_input))
    for element in json_input:
        # Add fastq_file to json
        link = element['Link'].split(";")
        forward_link = "ftp://" + link[0]
        reverse_link = "ftp://" + link[1]
        forward_name = link[0].split("/")[-1]
        reverse_name = link[1].split("/")[-1]
        forward_path = download_path + forward_name
        reverse_path = download_path + reverse_name
        element['fastq_file'] = "{};{}".format(forward_path,reverse_path)
        
        # Download fastq
        print("Downloading {}".format(forward_name))
        urllib.request.urlretrieve(forward_link,forward_path)
        print("Downloading {}".format(reverse_name))
        urllib.request.urlretrieve(reverse_link,reverse_path)

    # Create a new json file
    with open(path_output, 'w') as result:
        result.write(json.dumps(json_input))

def trim_data(path_input,path_output):
    json_input = json.load(open(path_input))
    for element in json_input:
        fastq_list = element['fastq_file'].split(";")
        forward = fastq_list[0]
        reverse = fastq_list[1]
        
        fastq_path = '/'.join(forward.split("/")[:-1])
        forward_name = forward.split("/")[-1]
        reverse_name = reverse.split("/")[-1]
        
        paired_forward = fastq_path + '/paired_' + forward_name
        paired_reverse = fastq_path + '/paired_' + reverse_name
        unpaired_forward = fastq_path + '/unpaired_' + forward_name
        unpaired_reverse = fastq_path + '/unpaired_' + reverse_name
        
        trim_execute = "java -jar /home/anhvu/bioinformatics_tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE {} {} {} {} {} {} ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36".format(forward,reverse,paired_forward,unpaired_forward,paired_reverse,unpaired_reverse)
        os.system(trim_execute)
        
        element['trimmed_data'] = "{};{}".format(paired_forward,paired_reverse)
        
    with open(path_output, 'w') as result:
        result.write(json.dumps(json_input))
        
if __name__ == '__main__':
    main()


#### FIRST STEP
#csv_to_json_converter(path_input,path_output)
#csv_to_json_converter("/home/anhvu/vu_data/Tri_new/ef.csv","/home/anhvu/vu_data/Tri_new/ef.json")

#### SECOND STEP
#sequence_download(path_input,download_path,path_output)
#sequence_download("/home/anhvu/vu_data/Tri_new/ef.json","/home/anhvu/vu_data/Tri_new/","/home/anhvu/vu_data/Tri_new/new_ef.json")

#### THIRD STEP
trim_data("/home/anhvu/vu_data/Tri_new/new_ef.json","/home/anhvu/vu_data/Tri_new/trimmed_ef.json")
