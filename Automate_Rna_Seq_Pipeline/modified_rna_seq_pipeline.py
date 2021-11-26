import os
from bs4 import BeautifulSoup
import subprocess
from tkinter import *

command_paths = {}
genome_gtf_path = ""
genome_fa_path = ""
genome_index = ""
bap = 0
main_folder_path = ""
text1 = ""
probar = ""
def main(se_C,se_T,pe_C,pe_T,main_folder_path_local,bap_local,text1_local,probar_local):
    global command_paths
    global genome_gtf_path
    global genome_fa_path
    global genome_index
    global bap
    global main_folder_path
    global text1
    global probar
    text1 = text1_local
    probar = probar_local
    main_folder_path = main_folder_path_local
    bap = bap_local
    command_paths = fetch_command_paths()
    genome_gtf_path, genome_fa_path = fetch_gtf_fa(main_folder_path)
    genome_index = bowtie(genome_fa_path)

    text1.insert(END,"OPERATION STARTED\n")
    probar.config(value=5)
    if len(se_C)!=0:
        text1.insert(END, "Working on single end files:\n")
        se_pipeline(se_C,se_T)
    if len(pe_C)!=0:
        text1.insert(END, "Working on paired end files:\n")
        pe_pipeline(pe_C,pe_T)

def se_pipeline(se_C,se_T):
    se_C_T = se_C + se_T
    for se_file in se_C_T:
        se_filename = se_file.split("/")[-1]
        folder_path = se_file.strip(".fastq")
        try:
            os.mkdir(folder_path)
        except:
            print("folder  already exist")
        mv_comm = "mv"+" "+se_file+" "+folder_path
        os.system(mv_comm)
        new_file_path = folder_path + "/" + se_filename
        fastqc_input = new_file_path
        fastqc_quality = "bad"
        cnt = 0
        while (fastqc_quality == "bad"):
            html_filepath = fastqc(fastqc_input,command_paths["fastqc"])
            fastqc_quality = quality_check(html_filepath)
            if fastqc_quality == "good":
                print("no need to cutadapt")
                
            elif fastqc_quality == "bad":
                cnt += 1
                cutadapt_out_file_path = cutadapt(fastqc_input, html_filepath, cnt, command_paths["cutadapt"])
                fastqc_input = cutadapt_out_file_path

        tophat_input_files = fastqc_input
        tophat_out_folder_path = tophat(tophat_input_files, genome_index, genome_gtf_path)
        cufflinks_out_folder_path = cufflinks(tophat_out_folder_path, genome_gtf_path,command_paths["cufflinks"])
        assembly_txt_filepath = assembly_txt(cufflinks_out_folder_path, main_folder_path)

    cuffmerge_out_folder_path = cuffmerge(assembly_txt_filepath, genome_gtf_path, genome_fa_path,command_paths["cuffmerge"])
    cuffdiff(cuffmerge_out_folder_path,bap,se_C,se_T,command_paths["cuffdiff"])

def generate_pairs(pe_C_T):
    text1.insert(END, "Generating distinct pairs for control and treated paired end files\n")

    pairs = []
    for i in range(len(pe_C_T) - 1):
        for j in range(i + 1, len(pe_C_T)):
            if pe_C_T[i].split("/")[-1].split("_")[0] == pe_C_T[j].split("/")[-1].split("_")[0]:
                pair = [pe_C_T[i], pe_C_T[j]]
                pairs.append(pair)
    return pairs

def pe_pipeline(pe_C,pe_T):
    pe_C_T = pe_C + pe_T
    pairs = generate_pairs(pe_C_T)
    for each_pair in pairs:
        tophat_input_files_list = []
        for each_pe_file in each_pair:
            pe_filename = each_pe_file.split("/")[-1]
            folder_name = pe_filename.split("_")[0]
            folder_path = each_pe_file.strip(pe_filename) + folder_name
            try:
                os.mkdir(folder_path)
            except:
                print("folder already exist")
            mv_comm = "mv" + " " + each_pe_file + " " + folder_path
            os.system(mv_comm)
            new_file_path = folder_path + "/" + pe_filename
            fastqc_input = new_file_path
            fastqc_quality = "bad"
            cnt = 0
            while (fastqc_quality == "bad"):
                html_filepath = fastqc(fastqc_input, command_paths["fastqc"])
                fastqc_quality = quality_check(html_filepath)
                if fastqc_quality == "good":
                    print("no need to cutadapt")
                    
                elif fastqc_quality == "bad":
                    cnt += 1
                    cutadapt_out_file_path = cutadapt(fastqc_input, html_filepath, cnt, command_paths["cutadapt"])
                    fastqc_input = cutadapt_out_file_path
            tophat_input_files_list.append(fastqc_input + " ")
        tophat_out_folder_path = tophat(tophat_input_files_list, genome_index, genome_gtf_path)
        cufflinks_out_folder_path = cufflinks(tophat_out_folder_path, genome_gtf_path,command_paths["cufflinks"])
        assembly_txt_filepath = assembly_txt(cufflinks_out_folder_path, main_folder_path)
    cuffmerge_out_folder_path = cuffmerge(assembly_txt_filepath, genome_gtf_path, genome_fa_path,command_paths["cuffmerge"])
    cuffdiff(cuffmerge_out_folder_path, bap, pe_C, pe_T,command_paths["cuffdiff"])

def fetch_gtf_fa(main_folder_path):
    text1.insert(END, "Fetching .gtf and .fa file\n")
    file_list = os.listdir(main_folder_path)
    for file in file_list:
        if file.endswith(".gtf"):
            genome_gtf_path = main_folder_path + "/" + file
        if file.endswith(".fa"):
            genome_fa_path = main_folder_path + "/" + file
    return genome_gtf_path, genome_fa_path

def bowtie(genome_fa_path):
    text1.insert(END, "Bowtie running on: "+genome_fa_path+"\n")
    probar.config(value=10)
    bowtie_out_path = genome_fa_path.strip(genome_fa_path.split("/")[-1])+"genome"
    bowtie_comm = "bowtie2-build " + genome_fa_path + " " + bowtie_out_path
    os.system(bowtie_comm)
    return bowtie_out_path

def fastqc(fastqc_input,fastqc_comm_path):
    text1.insert(END, "Fastqc running on: "+fastqc_input+"\n")
    probar.config(value=20)
    fastqc_comm = fastqc_comm_path + fastqc_input
    os.system(fastqc_comm)
    html_filepath = fastqc_input.strip(".fastq")+"_fastqc.html"
    return html_filepath

def quality_check(html_filepath):
    text1.insert(END, "Checking for quality of the generated .html file: "+html_filepath+"\n")
    html_filehandle = open(html_filepath)
    soup = BeautifulSoup(html_filehandle, 'lxml')
    summary = soup.find("div", class_="summary")
    for list in summary.find_all("li"):
        if list.a.text == "Overrepresented sequences":
            if list.img["alt"] == "[PASS]"or list.img["alt"] == "[WARNING]":
                fastqc_quality = "good"
                text1.insert(END, "No need to cutadapt\n")
                return fastqc_quality
            elif list.img["alt"] == "[FAIL]":
                fastqc_quality = "bad"
                text1.insert(END, "Cutadapt required\n")
                return fastqc_quality

def cutadapt(fastqc_input,html_filepath,cnt,cutadapt_comm_path):
    text1.insert(END, "Cutadapt running on: "+fastqc_input+"\n")
    probar.config(value=30)
    over_rep_seq = cutadapt_fetch(html_filepath)
    cutadapt_out_file_path = fastqc_input.strip(".fastq")+"_trimmed_"+cnt+".fastq"
    cutadapt_comm = cutadapt_comm_path + over_rep_seq + " -q30,30 -m20 -o " + cutadapt_out_file_path + " " + fastqc_input
    os.system(cutadapt_comm)
    return cutadapt_out_file_path

def cutadapt_fetch(html_filepath):
        soup = BeautifulSoup(open(html_filepath), "lxml")
        module_list = soup.find_all("div", class_="module")
        for module in module_list:
            if module.h2.text == "Overrepresented sequences":
                tbody = module.tbody
                tr_list = tbody.find_all("tr")
                over_rep_seq = ""
                for tr in tr_list:
                    seq = str(tr).split("<td>")[1].split("</td>")[0]
                    seq = "-b " + seq + " "
                    over_rep_seq += seq
        return over_rep_seq

def tophat(tophat_input_files,genome_index, genome_gtf_path):
    if type(tophat_input_files) == str:
        text1.insert(END, "Running tophat on: "+tophat_input_files+"\n")
        probar.config(value=40)
        tophat_input_files_str = tophat_input_files
        tophat_out_folder_path = tophat_input_files.strip(tophat_input_files.split("/")[-1]) + "tophat_" + tophat_input_files.split("/")[-1].strip(".fastq")
    elif type(tophat_input_files) == list:
        tophat_input_files_str = tophat_input_files[0] + " " + tophat_input_files[1]
        probar.config(value=40)
        text1.insert(END, "Running tophat on: " + tophat_input_files_str + "\n")
        tophat_out_folder_path = tophat_input_files[0].strip(tophat_input_files[0].split("/")[-1]) + "tophat_" + tophat_input_files[0].split("/")[-1].split("_")[0]
    tophat_comm = "tophat2 -o " + tophat_out_folder_path + " -G " + genome_gtf_path + " " + genome_index + " " + tophat_input_files_str
    os.system(tophat_comm)
    return tophat_out_folder_path

def cufflinks(tophat_out_folder_path, genome_gtf_path,cufflinks_comm_path):
    accepted_hits_bam_file_path = tophat_out_folder_path + "/" + "accepted_hits.bam"
    text1.insert(END, "Running cufflinks on: "+accepted_hits_bam_file_path+"\n")
    probar.config(value=60)
    cufflinks_out_folder_path = tophat_out_folder_path.strip(tophat_out_folder_path.split("/")[-1])+"cufflinks_"+tophat_out_folder_path.split("/")[-1].split("_")[-1]
    cufflinks_comm = cufflinks_comm_path + "-o" + " " + cufflinks_out_folder_path +" " + "-G" + " " + genome_gtf_path + " " + accepted_hits_bam_file_path
    os.system(cufflinks_comm)
    return cufflinks_out_folder_path

def assembly_txt(cufflinks_out_folder_path,main_folder_path):
    text1.insert(END, "Generating assembly.txt file\n")
    transcripts_gtf_file_path = cufflinks_out_folder_path + "/" + "transcripts.gtf"
    assembly_txt_filepath = main_folder_path + "/" + "assembly.txt"
    assembly_txt_file_handle = open(assembly_txt_filepath,"a")
    assembly_txt_file_handle.write(transcripts_gtf_file_path + "\n")
    assembly_txt_file_handle.close()
    return assembly_txt_filepath

def cuffmerge(assembly_txt_filepath,genome_gtf_path,genome_fa_path,cuffmerge_comm_path):
    text1.insert(END, "Running cuffmerge on: "+assembly_txt_filepath+"\n")
    probar.config(value=80)
    cuffmerge_out_folder_path = assembly_txt_filepath.strip(assembly_txt_filepath.split("/")[-1]) +  "cuffmerge_output"
    cuffmerge_comm = cuffmerge_comm_path + "-o" + " " + cuffmerge_out_folder_path + " -g " + genome_gtf_path + " -s " + genome_fa_path + " " + assembly_txt_filepath
    os.system(cuffmerge_comm)
    return cuffmerge_out_folder_path

def cuffdiff(cuffmerge_out_folder_path,bap,C_files_list,T_files_list,cuffdiff_comm_path):
    cuffdiff_out_folder_path = cuffmerge_out_folder_path.strip(cuffmerge_out_folder_path.split("/")[-1]) + "cuffdiff_output"
    merged_gtf_file_path = cuffmerge_out_folder_path + "/" + "merged.gtf"
    text1.insert(END, "Running cufdiff on: "+merged_gtf_file_path+"\n")
    all_C_accepted_hits_bam_set = set()
    all_T_accepted_hits_bam_set = set()

    for C_file in C_files_list:
        C_filename = C_file.split("/")[-1]
        if "_" not in C_filename:
            C_filename_first_half = C_filename.strip(".fastq")
            C_tophat_folder_path = C_file.strip(".fastq")+"/"+"tophat_"+C_filename_first_half
            C_accepted_hits_bam = C_tophat_folder_path + "/" + "accepted_hits.bam"
            all_C_accepted_hits_bam_set.add(C_accepted_hits_bam)
        elif "_" in C_filename:
            C_filename_first_half = C_filename.split("_")[0]
            C_tophat_folder_path = C_file.strip(C_file.split("/")[-1]) + C_filename_first_half+"/"+"tophat_"+C_filename_first_half
            C_accepted_hits_bam = C_tophat_folder_path + "/" + "accepted_hits.bam"
            all_C_accepted_hits_bam_set.add(C_accepted_hits_bam)
        all_C_accepted_hits_bam = ",".join(all_C_accepted_hits_bam_set)
    
    for T_file in T_files_list:
        T_filename = T_file.split("/")[-1]
        if "_" not in T_filename:
            T_filename_first_half = T_filename.strip(".fastq")
            T_tophat_folder_path = T_file.strip(".fastq")+"/"+"tophat_"+T_filename_first_half
            T_accepted_hits_bam = T_tophat_folder_path + "/" + "accepted_hits.bam"
            all_T_accepted_hits_bam_set.add(T_accepted_hits_bam)
        elif "_" in T_filename:
            T_filename_first_half = T_filename.split("_")[0]
            T_tophat_folder_path = T_file.strip(T_file.split("/")[-1]) + T_filename_first_half+"/"+"tophat_"+T_filename_first_half
            T_accepted_hits_bam = T_tophat_folder_path + "/" + "accepted_hits.bam"
            all_T_accepted_hits_bam_set.add(T_accepted_hits_bam)
        all_T_accepted_hits_bam = ",".join(all_T_accepted_hits_bam_set)
        
    if bap == 1:
        bap_path = all_C_accepted_hits_bam + " " + all_T_accepted_hits_bam
        bap_choice = "control,treated"
    elif bap == 0:
        bap_path = all_T_accepted_hits_bam + " " + all_C_accepted_hits_bam
        bap_choice = "treated,control"

    cuffdiff_comm = cuffdiff_comm_path + "-o" + " " + cuffdiff_out_folder_path + " " + "-L" + " " + bap_choice + " " + "-u" + " " + merged_gtf_file_path + " " + bap_path
    os.system(cuffdiff_comm)
    probar.config(value=95)

def fetch_command_paths():
    command_paths = {}
    fastqc_path_comm = "which fastqc".split(" ")
    temp = subprocess.Popen(fastqc_path_comm, stdout=subprocess.PIPE)
    fastqc_path_tup = temp.communicate()
    fastqc_comm_path = str(fastqc_path_tup[0], 'UTF8').strip("\n") + " "
    command_paths["fastqc"] = fastqc_comm_path

    cutadapt_path_comm = "which cutadapt".split(" ")
    temp = subprocess.Popen(cutadapt_path_comm, stdout=subprocess.PIPE)
    cutadapt_path_tup = temp.communicate()
    cutadapt_comm_path = str(cutadapt_path_tup[0], 'UTF8').strip("\n") + " "
    command_paths["cutadapt"] = cutadapt_comm_path

    cufflinks_path_comm = "which cufflinks".split(" ")
    temp = subprocess.Popen(cufflinks_path_comm, stdout=subprocess.PIPE)
    cufflinks_path_tup = temp.communicate()
    cufflinks_comm_path = str(cufflinks_path_tup[0], 'UTF8').strip("\n") + " "
    command_paths["cufflinks"] = cufflinks_comm_path

    cuffmerge_path_comm = "which cuffmerge".split(" ")
    temp = subprocess.Popen(cuffmerge_path_comm, stdout=subprocess.PIPE)
    cuffmerge_path_tup = temp.communicate()
    cuffmerge_comm_path = str(cuffmerge_path_tup[0], 'UTF8').strip("\n") + " "
    command_paths["cuffmerge"] = cuffmerge_comm_path

    cuffdiff_path_comm = "which cuffdiff".split(" ")
    temp = subprocess.Popen(cuffdiff_path_comm, stdout=subprocess.PIPE)
    cuffdiff_path_tup = temp.communicate()
    cuffdiff_comm_path = str(cuffdiff_path_tup[0], 'UTF8').strip("\n") + " "
    command_paths["cuffdiff"] = cuffdiff_comm_path

    return (command_paths)

#se_C = ["/home/ashish/Desktop/test1/c1.fastq"]
#se_T = ["/home/ashish/Desktop/test1/c2.fastq"]
#pe_C = []
#pe_T = []
#se_C=[]
#se_T=[]
#pe_C = ["/home/ashish/Desktop/test_rna_seq/SRR9665765_1.fastq",
 #       "/home/ashish/Desktop/test_rna_seq/SRR9665765_2.fastq"]
#pe_T = ["/home/ashish/Desktop/test_rna_seq/SRR9665767_1.fastq",
#        "/home/ashish/Desktop/test_rna_seq/SRR9665767_2.fastq"]
#main_folder_path = "/home/ashish/Desktop/pe_input_test1"
#bap = 1
#main(se_C,se_T,pe_C,pe_T,main_folder_path,bap)