from tkinter import *
from tkinter import messagebox
from tkinter import filedialog
from tkinter.ttk import *
from modified_rna_seq_pipeline import *


window=Tk()

window.geometry('350x400')
se_C = ()
se_T = ()
pe_C = ()
pe_T = ()
directory = ""

def browse_file():
    global se_C,se_T,pe_C,pe_T,directory

    directory = filedialog.askdirectory(initialdir="/home/ashish/Desktop",title='Please select a directory')
    #print(directory)
    choice = 1
    while(choice<=4):
        if choice==1:
            se_C = filedialog.askopenfilenames(initialdir=directory, title="Select single end control files")
            print(se_C)
            se_C_text.delete(1.0, END)
            for item in se_C:
                item = item.split("/")[-1]
                se_C_text.insert(END, item+" ")
        elif choice==2:
            se_T = filedialog.askopenfilenames(initialdir=directory, title="Select single end treated files")
            print(se_T)
            se_T_text.delete(1.0, END)
            for item in se_T:
                item = item.split("/")[-1]
                se_T_text.insert(END, item+" ")
        elif choice==3:
            pe_C = filedialog.askopenfilenames(initialdir=directory, title="Select paired end control files")
            print(pe_C)
            pe_C_text.delete(1.0, END)
            for item in pe_C:
                item = item.split("/")[-1]
                pe_C_text.insert(END, item+" ")
        elif choice==4:
            pe_T = filedialog.askopenfilenames(initialdir=directory, title="Select paired end treated files")
            print(pe_T)
            pe_T_text.delete(1.0, END)
            for item in pe_T:
                item = item.split("/")[-1]
                pe_T_text.insert(END, item+" ")
        choice += 1

def Go(se_C, se_T, pe_C, pe_T, directory):
    bap = (var.get())
    if bap == 0:
        messagebox.showerror("Warning!!","Biological analysis plan not selected")
    else:
        main(se_C, se_T, pe_C, pe_T, directory, bap,text1,probar)
        probar.config(value=100)
        messagebox.showinfo("Pipeline completes successfully", "output files generated in: " + directory)

#####header_and_footer_start######
header_label=Label(window,text="Rna Seq Pipeline",bg="black",fg="white")
header_label.pack(fill="x")
footer_label=Label(window,text="version 1.0",bg="black",fg="white")
footer_label.pack(side="bottom",fill="x")
#####header_and_footer_end######

#######Terminal_TEXT_FRAME_START######
text_frame=Frame(window,borderwidth=5,relief="sunken")
text_frame.pack(side="right",fill="y",padx=5)

#adding scrollbar to our text widget in our text_frame
scrollbar=Scrollbar(text_frame)
scrollbar.pack(side="right",fill="y")
#adding text widget to our text_frame
text1=Text(text_frame,width=27,yscrollcommand=scrollbar.set)
text1.pack(side="right",fill="y")
text1.insert(END,"Welcome To RNA Seq Pipeline"+"\n")

#connecting the scrollbar with our text wdget inside text_frame
scrollbar.config(command=text1.yview)
#######Terminal_TEXT_FRAME_END######


######image_frame_start######
image_frame=Frame(window,borderwidth=5,relief="sunken")
image_frame.pack(side="top",fill="x")
i=PhotoImage(file="/home/ashish/Desktop/Rna_seq_pipeline/dna.png")
image_label=Label(image_frame,image = i)
image_label.pack(anchor="w")
#######image_frame_end#######


######frame1_start######
frame1=Frame(window,relief="sunken")
frame1.pack(side="top",pady=2)

label1=Label(frame1,text="Click here to browse Input Folder")
label1.grid(row=0,column=0)

browse_button=Button(frame1,text="Browse",bg="black",fg="white",command = browse_file)
browse_button.grid(row=1,column=0)

probar= Progressbar(frame1,orient = HORIZONTAL,length=150,mode="determinate",value=0)
probar.grid(row=2,column=0)
######frame1_end######

######frame3_start######
frame3=Frame(window,borderwidth=5,relief="sunken")
frame3.pack(side="left")
label3=Label(frame3,text="Single end control files")
label3.grid(row=0,column=0,pady=2)
se_C_text=Text(frame3,height=1,width=10)
se_C_text.grid(row=0,column=1,padx=5,pady=2)

label4=Label(frame3,text="Single end treated files")
label4.grid(row=1,column=0,pady=2)
se_T_text=Text(frame3,height=1,width=10)
se_T_text.grid(row=1,column=1,padx=5,pady=2)

label5=Label(frame3,text="Paired end control files")
label5.grid(row=2,column=0,pady=2)
pe_C_text=Text(frame3,height=1,width=10)
pe_C_text.grid(row=2,column=1,padx=5,pady=2)

label6=Label(frame3,text="Paired end treated files")
label6.grid(row=3,column=0,pady=2)
pe_T_text=Text(frame3,height=1,width=10)
pe_T_text.grid(row=3,column=1,padx=5,pady=2)


######frame3_end#######

######frame2_start######
frame2=Frame(window,borderwidth=5,relief="sunken")
frame2.pack(side="right",fill="x",pady=2)
label2=Label(frame2,text="Select biogical analysis plan :")
label2.pack(anchor="w")
var=IntVar()
rb_dict={"C vs T":1,"T vs C":2}
for text,value in rb_dict.items():
    Radiobutton(frame2,text=text,variable=var,value=value,bg="white",fg="black",width=6).pack(side="top",anchor="w")

go_button = Button(frame2,text="Go",bg="white",fg="black", borderwidth=2,command = lambda :Go(se_C, se_T, pe_C, pe_T, directory))
go_button.pack(side="right", padx=2, pady=2)

######frame2_end#######


window.mainloop()