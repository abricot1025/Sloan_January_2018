# Run this file if you want to work with the graphical user interface
# No modification made by Clement, so same version since Summer 2017

import Tkinter as tk
import tkFileDialog as filedialog
#import Tkinter.filedialog as filedialog
import path_generator as pg

class Application(tk.Frame):

    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        self.pack()
        self.createWidgets()

    def createWidgets(self):
        self.if_animation = tk.IntVar()
        self.if_waveforms = tk.IntVar()

        self.text_upload_config = tk.Label(self, text="Configuration file").grid(row=0,sticky='e',pady = 20)
        self.upload_config_file = tk.Button(self, text="Load...", command=self.load_config_file).grid(row=0, column=1, sticky='w')

        self.text_upload_target = tk.Label(self, text="Target file").grid(row=1,sticky='e',pady = 20)
        self.upload_target_file = tk.Button(self, text="Load...", command=self.load_target_file).grid(row=1, column=1, sticky='w')

        #self.text_number_circles = tk.Label(self, text="Number of Circles").grid(row=2,sticky='e',pady = 20)
        #self.upload_target_file = tk.Entry(self).grid(row=2, column=1, sticky='w')

        self.if_animation_box = tk.Checkbutton(self,variable=self.if_animation, text="I want to see the positioners moving").grid(row=2, pady = 10,padx = 10,sticky='w')

        self.if_waveforms_box = tk.Checkbutton(self,variable=self.if_waveforms, text="I want to generate the motor waveforms").grid(row=3, pady = 10,padx = 10,sticky='w')

        self.result_dir_button = tk.Button(self, text='result folder', command=self.askdirectory).grid(row =4,sticky='w',padx = 10)

        self.go_button = tk.Button(self,text='GO',command=self.go_function ).grid(row=5,pady = 30)

    def load_config_file(self):

        self.config_file_name = filedialog.askopenfilename(filetypes=(("cfg files", "*.cfg"),
                                           ("HTML files", "*.html;*.htm"),
                                           ("All files", "*.*") ))
    def load_target_file(self):

        self.target_file_name = filedialog.askopenfilename(filetypes=(("text files", "*.txt"),
                                           ("HTML files", "*.html;*.htm"),
                                           ("All files", "*.*") ))

    def askdirectory(self):

        self.result_dir = filedialog.askdirectory()

    def go_function(self):
        try:
            pg.path_generator(self.config_file_name,self.target_file_name,self.result_dir,self.if_animation.get(),self.if_waveforms.get())
        # files are not given
        except Exception as e:
            print(e,'''Files are not properly selected.
            Check if you have chosen the configuration and target files.
            Check again the result directory.''')

if __name__=='__main__':
    root = tk.Tk()
    app = Application(master=root)
    app.mainloop()

