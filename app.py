"""
Main application class for Sequence Alignment Visualization Tool
"""

import tkinter as tk
from startpage import StartPage
from needleman_wunsch import PageOne
from smith_waterman import PageTwo

class AlgorithmApp(tk.Tk):
    """
    Main application class that contains the frame container and
    handles navigation between different algorithm pages.
    """

    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        self.title("Sequence Alignment Visualization Tool")
        self.geometry("1000x700")

        # The container is where we'll stack frames on top of each other
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}
        for F in (StartPage, PageOne, PageTwo):
            page_name = F.__name__
            frame = F(parent=container, controller=self)
            self.frames[page_name] = frame

            # Put all frames in the same location
            # The one on top will be visible
            frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame("StartPage")

    def show_frame(self, page_name):
        """Show a frame for the given page name"""
        frame = self.frames[page_name]
        frame.tkraise()
