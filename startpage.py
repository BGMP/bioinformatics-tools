"""
StartPage for Sequence Alignment Visualization Tool
"""

import tkinter as tk
from tkinter import Frame, Label, Button

class StartPage(tk.Frame):
    """
    Start page with buttons to select between Needleman-Wunsch and Smith-Waterman algorithms
    """

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller

        # Create main frame with padding
        main_frame = Frame(self, padx=20, pady=20)
        main_frame.pack(expand=True, fill="both")

        # Create title label
        title_label = Label(main_frame, text="Sequence Alignment Visualization Tool",
                          font=("Helvetica", 18, "bold"))
        title_label.pack(pady=20)

        # Create description
        desc_text = "An educational tool for understanding bioinformatics sequence alignment algorithms"
        desc_label = Label(main_frame, text=desc_text, font=("Helvetica", 12))
        desc_label.pack(pady=10)

        # Create frame for algorithm selection
        select_frame = Frame(main_frame)
        select_frame.pack(pady=20)

        # Label for algorithm selection
        select_label = Label(select_frame, text="Select an algorithm to visualize:",
                           font=("Helvetica", 12, "bold"))
        select_label.grid(row=0, column=0, columnspan=2, pady=10)

        # Create buttons for each algorithm
        button1 = Button(select_frame, text="Needleman-Wunsch Algorithm",
                       font=("Helvetica", 12),
                       width=25, height=2,
                       command=lambda: controller.show_frame("PageOne"))

        button2 = Button(select_frame, text="Smith-Waterman Algorithm",
                       font=("Helvetica", 12),
                       width=25, height=2,
                       command=lambda: controller.show_frame("PageTwo"))

        # Layout the buttons
        button1.grid(row=1, column=0, padx=10, pady=10)
        button2.grid(row=1, column=1, padx=10, pady=10)

        # Create frame for algorithm information
        info_frame = Frame(main_frame, relief="ridge", bd=1)
        info_frame.pack(pady=20, fill="x")

        # Add information about the algorithms
        info_title = Label(info_frame, text="Algorithm Information",
                         font=("Helvetica", 12, "bold"))
        info_title.pack(pady=(10, 5))

        nw_title = Label(info_frame, text="Needleman-Wunsch",
                        font=("Helvetica", 10, "bold"))
        nw_title.pack(pady=(5, 0))

        nw_info = Label(info_frame, text="Global sequence alignment algorithm that optimizes the alignment\n"
                                       "across the entire length of two sequences. Used when comparing\n"
                                       "sequences that are expected to be similar over their entire length.")
        nw_info.pack(pady=(0, 10))

        sw_title = Label(info_frame, text="Smith-Waterman",
                        font=("Helvetica", 10, "bold"))
        sw_title.pack(pady=(5, 0))

        sw_info = Label(info_frame, text="Local sequence alignment algorithm that identifies regions of similarity\n"
                                       "within sequences. Used for finding conserved patterns or motifs that\n"
                                       "may be shared between otherwise unrelated sequences.")
        sw_info.pack(pady=(0, 10))

        # Add example sequences section
        example_frame = Frame(main_frame)
        example_frame.pack(pady=10, fill="x")

        example_title = Label(example_frame, text="Example Sequences",
                            font=("Helvetica", 11, "bold"))
        example_title.pack(anchor="w")

        dna_label = Label(example_frame, text="DNA: ACGTACGT and ACGTACGA")
        dna_label.pack(anchor="w")

        protein_label = Label(example_frame, text="Protein: HEAGAWGHEE and PAWHEAE")
        protein_label.pack(anchor="w")

        # Credits
        credits = Label(main_frame, text="Sequence Alignment Visualization Tool",
                      font=("Helvetica", 8))
        credits.pack(side="bottom", pady=10)
