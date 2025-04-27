"""
StartPage for Sequence Alignment Visualization Tool
"""

import tkinter as tk


class StartPage(tk.Frame):
    """
    Start page with buttons to select between Needleman-Wunsch and Smith-Waterman algorithms
    """

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller

        # Create title label
        title_label = tk.Label(self, text="Sequence Alignment Visualization Tool",
                               font=("Helvetica", 18, "bold"))
        title_label.pack(pady=20)

        # Create description
        desc_text = "Select an algorithm to visualize sequence alignment:"
        desc_label = tk.Label(self, text=desc_text, font=("Helvetica", 12))
        desc_label.pack(pady=10)

        # Create frame for buttons
        button_frame = tk.Frame(self)
        button_frame.pack(pady=20)

        # Create buttons for each algorithm
        button1 = tk.Button(button_frame, text="Needleman-Wunsch Algorithm",
                            font=("Helvetica", 12),
                            width=25, height=2,
                            command=lambda: controller.show_frame("PageOne"))

        button2 = tk.Button(button_frame, text="Smith-Waterman Algorithm",
                            font=("Helvetica", 12),
                            width=25, height=2,
                            command=lambda: controller.show_frame("PageTwo"))

        # Add information labels about each algorithm
        nw_info = tk.Label(self, text="Needleman-Wunsch: Global sequence alignment",
                           font=("Helvetica", 10))
        sw_info = tk.Label(self, text="Smith-Waterman: Local sequence alignment",
                           font=("Helvetica", 10))

        # Layout the buttons and labels
        button1.grid(row=0, column=0, padx=10, pady=10)
        button2.grid(row=0, column=1, padx=10, pady=10)
        nw_info.pack(pady=(0, 5))
        sw_info.pack(pady=(0, 5))

        # Credits
        credits = tk.Label(self, text="Sequence Alignment Visualization Tool",
                           font=("Helvetica", 8))
        credits.pack(side="bottom", pady=10)
