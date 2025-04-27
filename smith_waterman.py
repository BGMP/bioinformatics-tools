"""
Smith-Waterman algorithm implementation for sequence alignment visualization
"""

import tkinter as tk
from tkinter import Frame, Label, Entry, Button
import numpy as np
from numpy import unravel_index


class PageTwo(tk.Frame):
    """
    Page implementing the Smith-Waterman algorithm for local sequence alignment
    """

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller

        # Initialize parameters
        self.gap_penalty = -1
        self.match_award = 1
        self.mismatch_penalty = -1
        self.seq1 = ""
        self.seq2 = ""
        self.n = 0
        self.m = 0
        self.score = []
        self.npscore = []
        self.max_score = 0
        self.index = 0
        self.index2 = 0
        self.i = 0
        self.j = 0
        self.align1 = ""
        self.align2 = ""
        self.a = 0  # For tracking max position
        self.b = 0  # For tracking max position

        # Create frames
        frame1 = Frame(self, relief="solid", bd=1)
        frame1.pack(side="left", fill="both", expand=True)
        frame2 = Frame(self)
        frame2.pack(side="right", fill="both", expand=True)

        # Store frame1 as an instance variable for later use
        self.frame1 = frame1

        # Create input UI elements
        label1 = Label(frame2, text="Sequence 1")
        entry1 = Entry(frame2)
        label2 = Label(frame2, text="Sequence 2")
        entry2 = Entry(frame2)
        label3 = Label(frame2, text="Gap/Mismatch Penalty")
        entry3 = Entry(frame2)

        # Position input UI elements
        label1.grid(row=0, column=0)
        entry1.grid(row=0, column=1)
        label2.grid(row=1, column=0)
        entry2.grid(row=1, column=1)
        label3.grid(row=2, column=0)
        entry3.grid(row=2, column=1)

        # Create label for displaying alignments
        label = Label(frame2)
        label.grid(row=4, column=1)
        self.label = label  # Store as instance variable

        # Create buttons
        button1 = Button(frame2, text="Execute", width=10, command=lambda: self.initialize(entry1, entry2, entry3))
        button1.grid(row=3, column=1)

        button2 = Button(frame2, text="<", width=3, command=self.left_button)
        button2.grid(row=3, column=2)

        button3 = Button(frame2, text=">", width=3, command=self.right_button)
        button3.grid(row=3, column=3)

        button4 = Button(frame2, text="<<", width=3, command=self.left_end_button)
        button4.grid(row=4, column=2)

        button5 = Button(frame2, text=">>", width=3, command=self.right_end_button)
        button5.grid(row=4, column=3)

        # Navigation button back to start page
        button = Button(frame2, text="Go to the start page",
                        command=lambda: controller.show_frame("StartPage"))
        button.grid(row=5, column=1)

    def zeros(self, rows, cols):
        """Create a matrix filled with zeros"""
        retval = []
        for x in range(rows):
            retval.append([])
            for y in range(cols):
                retval[-1].append(0)
        return retval

    def match_score(self, alpha, beta):
        """Calculate match score between two characters"""
        if alpha == beta:
            return self.match_award
        elif alpha == '-' or beta == '-':
            return self.gap_penalty
        else:
            return self.mismatch_penalty

    def smith_waterman(self, i, j):
        """Fill the matrix using Smith-Waterman algorithm"""
        match = self.score[i - 1][j - 1] + self.match_score(self.seq1[j - 1], self.seq2[i - 1])
        delete = self.score[i - 1][j] + self.gap_penalty
        insert = self.score[i][j - 1] + self.gap_penalty
        self.score[i][j] = max(0, match, delete, insert)  # Note the 0 to avoid negative scores

        self.show_matrix(self.score)

        # Highlight the reference cell of the score
        if self.score[i][j] == match and match > 0:
            e = Label(self.frame1, relief="solid", bd=1, fg="red")
            e.config(text=str(self.score[i - 1][j - 1]))
            e.grid(row=i, column=j, sticky="nsew")
        if self.score[i][j] == delete and delete > 0:
            e = Label(self.frame1, relief="solid", bd=1, fg="red")
            e.config(text=str(self.score[i - 1][j]))
            e.grid(row=i, column=j + 1, sticky="nsew")
        if self.score[i][j] == insert and insert > 0:
            e = Label(self.frame1, relief="solid", bd=1, fg="red")
            e.config(text=str(self.score[i][j - 1]))
            e.grid(row=i + 1, column=j, sticky="nsew")

    def traceback(self, i, j):
        """Traceback to find the optimal local alignment"""
        if self.score[self.i][self.j] == 0:  # If zero element has been reached
            self.label.config(text=self.align1[::-1] + '\n' + self.align2[::-1])
        elif self.i > 0 and self.j > 0:
            score_current = self.score[self.i][self.j]
            score_diagonal = self.score[self.i - 1][self.j - 1]
            score_up = self.score[self.i][self.j - 1]
            score_left = self.score[self.i - 1][self.j]

            if score_current == score_diagonal + self.match_score(self.seq1[self.j - 1], self.seq2[self.i - 1]):
                self.align1 += self.seq1[self.j - 1]
                self.align2 += self.seq2[self.i - 1]
                e = Label(self.frame1, relief="solid", bd=1, bg="green")
                e.config(text=str(score_diagonal))
                e.grid(row=self.i, column=self.j, sticky="nsew")
                self.i -= 1
                self.j -= 1
                self.label.config(text=self.align1[::-1] + '\n' + self.align2[::-1])
            elif score_current == score_up + self.gap_penalty:
                self.align1 += self.seq1[self.j - 1]
                self.align2 += '-'
                e = Label(self.frame1, relief="solid", bd=1, bg="green")
                e.config(text=str(score_up))
                e.grid(row=self.i + 1, column=self.j, sticky="nsew")
                self.j -= 1
                self.label.config(text=self.align1[::-1] + '\n' + self.align2[::-1])
            elif score_current == score_left + self.gap_penalty:
                self.align1 += '-'
                self.align2 += self.seq2[self.i - 1]
                e = Label(self.frame1, relief="solid", bd=1, bg="green")
                e.config(text=str(score_left))
                e.grid(row=self.i, column=self.j + 1, sticky="nsew")
                self.i -= 1
                self.label.config(text=self.align1[::-1] + '\n' + self.align2[::-1])
        return self.i, self.j

    def delete_traceback(self):
        """Delete trace for the result"""
        score_current = self.score[self.i][self.j]
        if self.i < self.a and self.j < self.b:  # For deleting traceback
            if self.align2[-1] == '-':
                self.align1 = self.align1[:len(self.align1) - 1]
                self.align2 = self.align2[:len(self.align2) - 1]
                e = Label(self.frame1, relief="solid", bd=1)
                e.config(text=str(score_current))
                e.grid(row=self.i + 1, column=self.j + 1, sticky="nsew")
                self.j += 1
                self.label.config(text=self.align1[::-1] + '\n' + self.align2[::-1])
            elif self.align1[-1] == '-':
                self.align1 = self.align1[:len(self.align1) - 1]
                self.align2 = self.align2[:len(self.align2) - 1]
                e = Label(self.frame1, relief="solid", bd=1)
                e.config(text=str(score_current))
                e.grid(row=self.i + 1, column=self.j + 1, sticky="nsew")
                self.i += 1
                self.label.config(text=self.align1[::-1] + '\n' + self.align2[::-1])
            else:
                self.align1 = self.align1[:len(self.align1) - 1]
                self.align2 = self.align2[:len(self.align2) - 1]
                e = Label(self.frame1, relief="solid", bd=1)
                e.config(text=str(score_current))
                e.grid(row=self.i + 1, column=self.j + 1, sticky="nsew")
                self.i += 1
                self.j += 1
                self.label.config(text=self.align1[::-1] + '\n' + self.align2[::-1])
        else:
            e = Label(self.frame1, relief="solid", bd=1)
            e.config(text=str(score_current))
            e.grid(row=self.i + 1, column=self.j + 1, sticky="nsew")
            self.index = len(self.seq1)
            self.index2 = len(self.seq2)
            f = Label(self.frame1, relief="solid", bd=1, bg="yellow")
            f.config(text=str(self.score[self.m][self.n]))
            f.grid(row=self.m + 1, column=self.n + 1, sticky="nsew")
        return self.index2, self.index

    def all_children(self, window):
        """Get all child widgets of a window"""
        _list = window.winfo_children()
        for item in _list:
            if item.winfo_children():
                _list.extend(item.winfo_children())
        return _list

    def initialize(self, entry1, entry2, entry3):
        """Initialize the algorithm with the input values"""
        self.index = 0
        self.index2 = 0
        widget_list = self.all_children(self.frame1)
        for item in widget_list:
            item.grid_forget()

        # Get input values
        self.seqA = entry1.get()
        self.seqB = entry2.get()
        self.seq1 = self.seqA.upper()
        self.seq2 = self.seqB.upper()
        self.gap_penalty = int(entry3.get())
        self.mismatch_penalty = int(entry3.get())
        self.n = len(self.seq1)
        self.m = len(self.seq2)
        self.i = self.m
        self.j = self.n
        self.score = self.zeros(self.m + 1, self.n + 1)
        self.align1 = ""  # delete previous result for a new alignment
        self.align2 = ""  # delete previous result for a new alignment

        # Display sequences in the UI
        for i in range(len(self.seq1)):  # location seq1
            e = Label(self.frame1)
            e.config(text=self.seq1[i])
            e.grid(row=0, column=i + 2)

        for i in range(len(self.seq2)):  # location seq2
            e = Label(self.frame1)
            e.config(text=self.seq2[i])
            e.grid(row=i + 2, column=0)

        # Initialize scoring matrix (all zeros for Smith-Waterman)
        for i in range(0, self.m + 1):
            self.score[i][0] = 0
        for j in range(0, self.n + 1):
            self.score[0][j] = 0

        self.label.config(text='')
        self.show_matrix(self.score)

    def right_button(self):
        """Handle right button click"""
        if self.index <= self.n and self.index2 <= self.m:
            self.index = self.index + 1
            if self.index % self.n == 1:  # Renew the column index
                self.index = 1
                self.index2 = self.index2 + 1
        self.button_event_right()

    def left_button(self):
        """Handle left button click"""
        if self.index >= 1 and self.index2 >= 1 and self.index2 != self.m + 1:  # Only for matrix element
            self.index = self.index - 1
            if self.index == 0 and self.index2 != 1:  # If at first column but not in first row
                self.index = self.n
                self.index2 = self.index2 - 1
        self.button_event_left()

    def button_event_right(self):
        """Handle right button event logic"""
        if self.index <= self.n and self.index2 <= self.m:
            self.smith_waterman(self.index2, self.index)
        elif self.index == 1 and self.index2 == self.m + 1:
            # Find maximum value in the matrix and start traceback from there
            self.npscore = np.array(self.score)
            (self.i, self.j) = unravel_index(self.npscore.argmax(), self.npscore.shape)
            (self.a, self.b) = (self.i, self.j)  # Store max position for later use
            e = Label(self.frame1, relief="solid", bd=1, bg="yellow")
            e.config(text=str(self.score[self.i][self.j]))
            e.grid(row=self.i + 1, column=self.j + 1, sticky="nsew")
            self.traceback(self.a, self.b)  # Trace from max value
            if self.i + 1 != self.m or self.j + 1 != self.n:
                f = Label(self.frame1, relief="solid", bd=1)
                f.config(text=str(self.score[self.m][self.n]))
                f.grid(row=self.m + 1, column=self.n + 1, sticky="nsew")
            self.index = 1
            self.index2 = self.m + 2
        else:
            self.traceback(self.i, self.j)

    def button_event_left(self):
        """Handle left button event logic"""
        if self.index <= self.n and self.index2 <= self.m:  # Only for matrix elements
            if self.index == 0 and self.index2 == 1:  # Do nothing if at first element
                self.index = 1
                self.index2 = 1
            else:
                self.smith_waterman(self.index2, self.index)
        else:
            self.delete_traceback()

    def tksleep(self, t):
        """Emulate time.sleep() while allowing Tkinter to update"""
        ms = int(t * 1000)
        # Use self.winfo_toplevel() instead of _get_default_root()
        root = self.winfo_toplevel()
        var = tk.IntVar(root)
        root.after(ms, lambda: var.set(1))
        root.wait_variable(var)

    def right_end_button(self):
        """Process the entire matrix and alignment in one go (forward)"""
        z = (self.m * self.n) + max(self.m, self.n)
        for k in range(z):
            self.right_button()
            self.tksleep(0.5)

    def left_end_button(self):
        """Process the entire matrix and alignment in one go (backward)"""
        z = (self.m * self.n) + max(self.m, self.n)
        if self.index != 0 and self.index2 != 0:
            for k in range(z):
                self.left_button()
                self.tksleep(0.5)

    def show_matrix(self, score):
        """Display the scoring matrix"""
        entry = {}
        # Create table of widgets
        for row in range(len(score)):
            for column in range(len(score[0])):
                index = (row, column)
                e = Label(self.frame1, relief="solid", bd=1)
                e.config(text=str(score[row][column]))
                e.grid(row=row + 1, column=column + 1, sticky="nsew")
                entry[index] = e
                # Highlight current cell
        e = Label(self.frame1, relief="solid", bd=1, bg="yellow")
        e.config(text=str(score[self.index2][self.index]))
        e.grid(row=self.index2 + 1, column=self.index + 1, sticky="nsew")
