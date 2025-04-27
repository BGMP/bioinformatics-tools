"""
Needleman-Wunsch algorithm implementation for sequence alignment visualization
"""

import tkinter as tk
from tkinter import Frame, Label, Entry, Button


class PageOne(tk.Frame):
    """
    Page implementing the Needleman-Wunsch algorithm for global sequence alignment
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
        self.index = 0
        self.index2 = 0
        self.i = 0
        self.j = 0
        self.align1 = ""
        self.align2 = ""

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

    def needleman_wunsch(self, i, j):
        """Fill the matrix using Needleman-Wunsch algorithm"""
        match = self.score[i - 1][j - 1] + self.match_score(self.seq1[j - 1], self.seq2[i - 1])
        delete = self.score[i - 1][j] + self.gap_penalty
        insert = self.score[i][j - 1] + self.gap_penalty
        self.score[i][j] = max(match, delete, insert)

        self.show_matrix(self.score)

        # Highlight the reference cell of the score
        if self.score[i][j] == match:
            e = Label(self.frame1, relief="solid", bd=1, fg="red")
            e.config(text=str(self.score[i - 1][j - 1]))
            e.grid(row=i, column=j, sticky="nsew")
        if self.score[i][j] == delete:
            e = Label(self.frame1, relief="solid", bd=1, fg="red")
            e.config(text=str(self.score[i - 1][j]))
            e.grid(row=i, column=j + 1, sticky="nsew")
        if self.score[i][j] == insert:
            e = Label(self.frame1, relief="solid", bd=1, fg="red")
            e.config(text=str(self.score[i][j - 1]))
            e.grid(row=i + 1, column=j, sticky="nsew")

    def traceback(self):
        """Traceback to find the optimal alignment"""
        score_current = self.score[self.i][self.j]
        score_diagonal = self.score[self.i - 1][self.j - 1]
        score_up = self.score[self.i][self.j - 1]
        score_left = self.score[self.i - 1][self.j]

        if self.i > 0 and self.j > 0:
            self.label.config(text=self.align1[::-1] + '\n' + self.align2[::-1])
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
        elif self.j > 0:
            score_current = self.score[self.i][self.j - 1]
            self.align1 += self.seq1[self.j - 1]
            self.align2 += '-'
            e = Label(self.frame1, relief="solid", bd=1, bg="green")
            e.config(text=str(score_current))
            e.grid(row=self.i + 1, column=self.j, sticky="nsew")
            self.j -= 1
            self.label.config(text=self.align1[::-1] + '\n' + self.align2[::-1])
        elif self.i > 0:
            score_current = self.score[self.i - 1][self.j]
            self.align1 += '-'
            self.align2 += self.seq2[self.i - 1]
            e = Label(self.frame1, relief="solid", bd=1, bg="green")
            e.config(text=str(score_current))
            e.grid(row=self.i, column=self.j + 1, sticky="nsew")
            self.i -= 1
            self.label.config(text=self.align1[::-1] + '\n' + self.align2[::-1])
        else:
            self.label.config(text=self.align1[::-1] + '\n' + self.align2[::-1])

    def delete_traceback(self):
        """Delete trace for the result"""
        score_current = self.score[self.i][self.j]
        if self.i < len(self.seq2) and self.j < len(self.seq1):
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
        elif self.j < len(self.seq1):
            self.align1 = self.align1[:len(self.align1) - 1]
            self.align2 = self.align2[:len(self.align2) - 1]
            e = Label(self.frame1, relief="solid", bd=1)
            e.config(text=str(score_current))
            e.grid(row=self.i + 1, column=self.j + 1, sticky="nsew")
            self.j += 1
            self.label.config(text=self.align1[::-1] + '\n' + self.align2[::-1])
        elif self.i < len(self.seq2):
            self.align1 = self.align1[:len(self.align1) - 1]
            self.align2 = self.align2[:len(self.align2) - 1]
            e = Label(self.frame1, relief="solid", bd=1)
            e.config(text=str(score_current))
            e.grid(row=self.i + 1, column=self.j + 1, sticky="nsew")
            self.i += 1
            self.label.config(text=self.align1[::-1] + '\n' + self.align2[::-1])
        else:
            self.label.config(text=self.align1[::-1] + '\n' + self.align2[::-1])
            self.index = len(self.seq1)
            self.index2 = len(self.seq2)
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
        self.count_step = 0
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

        # Initialize scoring matrix
        for i in range(0, self.m + 1):
            self.score[i][0] = self.gap_penalty * i
        for j in range(0, self.n + 1):
            self.score[0][j] = self.gap_penalty * j

        self.label.config(text='')
        self.show_matrix(self.score)

    def right_button(self):
        """Handle right button click"""
        if self.index <= self.n and self.index2 <= self.m:
            self.index = self.index + 1
            if self.index % self.n == 1:  # Wrap to next row
                self.index = 1
                self.index2 = self.index2 + 1
        self.button_event_right()

    def left_button(self):
        """Handle left button click"""
        if self.index >= 1 and self.index2 >= 1 and self.index2 != self.m + 1:  # Only until the last row of matrix
            self.index = self.index - 1
            if self.index == 0 and self.index2 != 1:  # If at beginning of row but not in the first row
                self.index = self.n
                self.index2 = self.index2 - 1
        self.button_event_left()

    def button_event_right(self):
        """Handle right button event logic"""
        if self.index <= self.n and self.index2 <= self.m:
            self.needleman_wunsch(self.index2, self.index)
        else:
            self.traceback()

    def button_event_left(self):
        """Handle left button event logic"""
        if self.index <= self.n and self.index2 <= self.m:  # Only for matrix elements
            if self.index == 0 and self.index2 == 1:  # Nothing if at first element
                self.index = 1
                self.index2 = 1
            else:
                self.needleman_wunsch(self.index2, self.index)
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
        z = (self.m * self.n) + self.m + self.n - 1
        for k in range(z):
            self.right_button()
            self.tksleep(0.5)

    def left_end_button(self):
        """Process the entire matrix and alignment in one go (backward)"""
        z = (self.m * self.n) + self.m + self.n - 1
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
