"""
Needleman-Wunsch algorithm implementation for sequence alignment visualization
"""

import tkinter as tk
from tkinter import Frame, Label, Entry, Button, StringVar, IntVar, Canvas, Scrollbar
from improved_algorithm import compute_needleman_wunsch, get_traceback_needleman_wunsch

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
        self.score = None
        self.computation_steps = []
        self.current_step_index = 0
        self.animation_in_progress = False
        self.animation_completed_matrix = False
        self.traceback_path = []
        self.cell_size = 30  # Size of each cell in the grid

        # Create frames
        self.left_frame = Frame(self)
        self.left_frame.pack(side="left", fill="both", expand=True)

        # Create canvas with scrollbars for matrix visualization
        self.canvas_frame = Frame(self.left_frame, bd=1, relief="solid")
        self.canvas_frame.pack(side="top", fill="both", expand=True)

        self.h_scrollbar = Scrollbar(self.canvas_frame, orient="horizontal")
        self.h_scrollbar.pack(side="bottom", fill="x")

        self.v_scrollbar = Scrollbar(self.canvas_frame)
        self.v_scrollbar.pack(side="right", fill="y")

        self.canvas = Canvas(self.canvas_frame,
                           xscrollcommand=self.h_scrollbar.set,
                           yscrollcommand=self.v_scrollbar.set)
        self.canvas.pack(side="left", fill="both", expand=True)

        self.h_scrollbar.config(command=self.canvas.xview)
        self.v_scrollbar.config(command=self.canvas.yview)

        # Right side panel for controls
        right_frame = Frame(self)
        right_frame.pack(side="right", fill="y")

        # Create title for the algorithm page
        algorithm_title = Label(right_frame, text="Needleman-Wunsch Algorithm", font=("Helvetica", 14, "bold"))
        algorithm_title.grid(row=0, column=0, columnspan=4, pady=10)

        # Create description
        desc = Label(right_frame, text="Global sequence alignment algorithm", font=("Helvetica", 10, "italic"))
        desc.grid(row=1, column=0, columnspan=4)

        # Create input UI elements
        label1 = Label(right_frame, text="Sequence 1")
        self.entry1_var = StringVar()
        entry1 = Entry(right_frame, textvariable=self.entry1_var, width=30)

        label2 = Label(right_frame, text="Sequence 2")
        self.entry2_var = StringVar()
        entry2 = Entry(right_frame, textvariable=self.entry2_var, width=30)

        # Enhanced scoring parameters section with separate controls
        scoring_label = Label(right_frame, text="Scoring Parameters:", font=("Helvetica", 10, "bold"))

        match_label = Label(right_frame, text="Match Reward")
        self.match_var = IntVar(value=1)
        match_entry = Entry(right_frame, textvariable=self.match_var, width=5)

        mismatch_label = Label(right_frame, text="Mismatch Penalty")
        self.mismatch_var = IntVar(value=-1)
        mismatch_entry = Entry(right_frame, textvariable=self.mismatch_var, width=5)

        gap_label = Label(right_frame, text="Gap Penalty")
        self.gap_var = IntVar(value=-1)
        gap_entry = Entry(right_frame, textvariable=self.gap_var, width=5)

        # Position input UI elements
        label1.grid(row=3, column=0, sticky="w", pady=5)
        entry1.grid(row=3, column=1, columnspan=3, sticky="we", pady=5)

        label2.grid(row=4, column=0, sticky="w", pady=5)
        entry2.grid(row=4, column=1, columnspan=3, sticky="we", pady=5)

        # Position scoring parameters
        scoring_label.grid(row=7, column=0, columnspan=4, sticky="w", pady=(10, 5))

        match_label.grid(row=8, column=0, sticky="w", pady=2)
        match_entry.grid(row=8, column=1, sticky="w", pady=2)

        mismatch_label.grid(row=9, column=0, sticky="w", pady=2)
        mismatch_entry.grid(row=9, column=1, sticky="w", pady=2)

        gap_label.grid(row=10, column=0, sticky="w", pady=2)
        gap_entry.grid(row=10, column=1, sticky="w", pady=2)

        # Add a scoring examples section
        examples_label = Label(right_frame, text="Example scoring schemes:", font=("Helvetica", 9))
        examples_label.grid(row=11, column=0, columnspan=2, sticky="w", pady=(10, 0))

        dna_button = Button(right_frame, text="DNA", width=8,
                          command=lambda: self.set_scoring_scheme(2, -1, -2))
        dna_button.grid(row=12, column=0, sticky="w", pady=2)

        protein_button = Button(right_frame, text="BLOSUM62", width=8,
                              command=lambda: self.set_scoring_scheme(1, -1, -1))
        protein_button.grid(row=12, column=1, sticky="w", pady=2)

        custom_button = Button(right_frame, text="Custom", width=8,
                             command=lambda: self.set_scoring_scheme(5, -4, -8))
        custom_button.grid(row=12, column=2, sticky="w", pady=2)

        # Add a reset button
        reset_button = Button(right_frame, text="Reset", width=8,
                            command=lambda: self.reset_form())
        reset_button.grid(row=12, column=3, sticky="w", pady=2)

        # Add explanation section
        explanation = Label(right_frame, text="Current scoring matrix: Match = +1, Mismatch = -1, Gap = -1",
                         font=("Helvetica", 9))
        explanation.grid(row=13, column=0, columnspan=4, sticky="w", pady=(10, 5))
        self.explanation = explanation  # Save for later updates

        # Create label for displaying alignments
        self.result_label = Label(right_frame, text="", font=("Courier", 10), justify="left")
        self.result_label.grid(row=14, column=0, columnspan=4, pady=10, sticky="w")

        # Create execution button
        execute_button = Button(right_frame, text="Execute", width=10,
                              command=self.initialize)
        execute_button.grid(row=15, column=0, columnspan=2, pady=10)

        # Step navigation buttons
        self.prev_button = Button(right_frame, text="<", width=3, command=self.previous_step, state="disabled")
        self.prev_button.grid(row=15, column=2, pady=10)

        self.next_button = Button(right_frame, text=">", width=3, command=self.next_step, state="disabled")
        self.next_button.grid(row=15, column=3, pady=10)

        # Progress label
        self.progress_label = Label(right_frame, text="", font=("Helvetica", 9))
        self.progress_label.grid(row=16, column=0, columnspan=4, pady=5)

        # Animation control buttons
        self.start_button = Button(right_frame, text="<<", width=3, command=self.go_to_start, state="disabled")
        self.start_button.grid(row=17, column=2, pady=5)

        self.end_button = Button(right_frame, text=">>", width=3, command=self.animate_to_end, state="disabled")
        self.end_button.grid(row=17, column=3, pady=5)

        # Animation speed control
        speed_label = Label(right_frame, text="Animation Speed:", font=("Helvetica", 9))
        speed_label.grid(row=18, column=0, sticky="w")

        self.speed_var = IntVar(value=2)  # Default medium speed
        speed_slow = Button(right_frame, text="Slow", width=6,
                           command=lambda: self.set_animation_speed(3))
        speed_slow.grid(row=18, column=1, sticky="w")

        speed_medium = Button(right_frame, text="Medium", width=6,
                             command=lambda: self.set_animation_speed(2))
        speed_medium.grid(row=18, column=2, sticky="w")

        speed_fast = Button(right_frame, text="Fast", width=6,
                           command=lambda: self.set_animation_speed(1))
        speed_fast.grid(row=18, column=3, sticky="w")

        # Navigation button back to start page
        button = Button(right_frame, text="Go to the start page",
                      command=lambda: controller.show_frame("StartPage"))
        button.grid(row=19, column=0, columnspan=4, pady=10)

        # Store entry variables for later use
        self.entry1 = entry1
        self.entry2 = entry2
        self.match_entry = match_entry
        self.mismatch_entry = mismatch_entry
        self.gap_entry = gap_entry

    def set_animation_speed(self, speed):
        """Set the animation speed (1=fast, 2=medium, 3=slow)"""
        self.speed_var.set(speed)

    def set_scoring_scheme(self, match, mismatch, gap):
        """Set a predefined scoring scheme"""
        self.match_var.set(match)
        self.mismatch_var.set(mismatch)
        self.gap_var.set(gap)
        self.update_explanation()

    def reset_form(self):
        """Reset the form to default values"""
        # Stop any ongoing animation
        if self.animation_in_progress:
            self.animation_in_progress = False
            self.canvas.after_cancel(self.animation_id)

        self.entry1_var.set("")
        self.entry2_var.set("")
        self.match_var.set(1)
        self.mismatch_var.set(-1)
        self.gap_var.set(-1)
        self.update_explanation()
        self.canvas.delete("all")
        self.result_label.config(text="")
        self.progress_label.config(text="")
        self.prev_button.config(state="disabled")
        self.next_button.config(state="disabled")
        self.start_button.config(state="disabled")
        self.end_button.config(state="disabled")

    def update_explanation(self):
        """Update the explanation text based on current scoring parameters"""
        text = f"Current scoring matrix: Match = +{self.match_var.get()}, "
        text += f"Mismatch = {self.mismatch_var.get()}, Gap = {self.gap_var.get()}"
        self.explanation.config(text=text)

    def initialize(self):
        """Initialize the algorithm with the input values"""
        # Stop any ongoing animation
        if hasattr(self, 'animation_id') and self.animation_in_progress:
            self.animation_in_progress = False
            self.canvas.after_cancel(self.animation_id)

        # Clear canvas
        self.canvas.delete("all")

        # Get input values
        self.seq1 = self.entry1_var.get().upper()
        self.seq2 = self.entry2_var.get().upper()

        if len(self.seq1) == 0 or len(self.seq2) == 0:
            self.result_label.config(text="Please enter both sequences.")
            return

        # Get scoring parameters from the UI
        self.match_award = self.match_var.get()
        self.mismatch_penalty = self.mismatch_var.get()
        self.gap_penalty = self.gap_var.get()

        self.n = len(self.seq1)
        self.m = len(self.seq2)

        # Reset state
        self.current_step_index = 0
        self.traceback_started = False
        self.highlighted_path = []
        self.animation_completed_matrix = False

        # Enable/disable navigation buttons
        self.prev_button.config(state="disabled")
        self.next_button.config(state="normal")
        self.start_button.config(state="disabled")
        self.end_button.config(state="normal")

        # Compute the matrix and get steps
        self.score, self.computation_steps = compute_needleman_wunsch(
            self.seq1, self.seq2, self.match_award, self.mismatch_penalty, self.gap_penalty
        )

        # Create the matrix visualization
        self.create_matrix_visualization()

        self.progress_label.config(text=f"Matrix initialized. Ready to start computation.")
        self.result_label.config(text="")

    def create_matrix_visualization(self):
        """Create the initial visualization of the matrix"""
        self.canvas.delete("all")

        # Calculate canvas size with extra space for shifted letters
        width = (self.n + 2) * self.cell_size + 1  # Extra space for the shift
        height = (self.m + 2) * self.cell_size

        self.canvas.config(scrollregion=(0, 0, width, height))

        # Draw sequence labels - shifted 100 pixels to the right
        for i in range(self.n):
            # Center the sequence 1 letters in the middle of their cells, shifted 100px right
            self.canvas.create_text(
                (i+2) * self.cell_size - self.cell_size/2 + 30, self.cell_size/2,
                text=self.seq1[i], font=("Helvetica", 10, "bold")
            )

        for i in range(self.m):
            # Center the sequence 2 letters in the middle of their cells, shifted 100px right
            self.canvas.create_text(
                self.cell_size/2, (i+2) * self.cell_size - self.cell_size/2 + 30,
                text=self.seq2[i], font=("Helvetica", 10, "bold")
            )

        # Draw the matrix cells - keeping the matrix in the same position
        for i in range(self.m+1):
            for j in range(self.n+1):
                x = (j+1) * self.cell_size
                y = (i+1) * self.cell_size

                # Draw cell background
                self.canvas.create_rectangle(
                    x, y, x + self.cell_size, y + self.cell_size,
                    fill="white", outline="black", tags=f"cell_{i}_{j}"
                )

                if i == 0 or j == 0:
                    value = 0
                    if i == 0 and j > 0:
                        value = self.gap_penalty * j
                    elif j == 0 and i > 0:
                        value = self.gap_penalty * i

                    # Draw the value
                    self.canvas.create_text(
                        x + self.cell_size/2, y + self.cell_size/2,
                        text=str(value), tags=f"text_{i}_{j}"
                    )

    def highlight_cell(self, i, j, color="#D6EAF8"):  # Light blue highlight
        """Highlight a specific cell in the matrix"""
        # Calculate cell position
        x = (j+1) * self.cell_size
        y = (i+1) * self.cell_size

        # Create or update the highlight rectangle
        tag = f"highlight_{i}_{j}"
        self.canvas.delete(tag)

        self.canvas.create_rectangle(
            x, y, x + self.cell_size, y + self.cell_size,
            fill=color, outline="black", tags=[tag, "highlight"]
        )

        # Make sure the text is on top
        self.canvas.tag_raise(f"text_{i}_{j}")

    def update_cell_value(self, i, j, value):
        """Update the value displayed in a cell"""
        # Calculate cell position
        x = (j+1) * self.cell_size
        y = (i+1) * self.cell_size

        # Delete old text
        self.canvas.delete(f"text_{i}_{j}")

        # Create new text
        self.canvas.create_text(
            x + self.cell_size/2, y + self.cell_size/2,
            text=str(value), tags=f"text_{i}_{j}"
        )

    def next_step(self):
        """Process the next step in the algorithm"""
        if self.traceback_started:
            # We're in traceback mode
            self.show_more_traceback()
            return

        if self.current_step_index < len(self.computation_steps):
            # Show next computation step
            step = self.computation_steps[self.current_step_index]
            i, j = step['i'], step['j']

            # Update the value and highlight the current cell
            self.update_cell_value(i, j, step['score'])
            self.highlight_cell(i, j)

            # Update progress
            self.progress_label.config(text=f"Step {self.current_step_index+1}/{len(self.computation_steps)}: "
                                         f"Computing cell ({i}, {j})")

            # Increment step counter
            self.current_step_index += 1

            # Enable previous button
            self.prev_button.config(state="normal")
            self.start_button.config(state="normal")

            # Check if we've completed the matrix
            if self.current_step_index >= len(self.computation_steps):
                self.animation_completed_matrix = True

                # Automatically start traceback if we were animating
                if self.animation_in_progress:
                    self.start_traceback()

        else:
            # Matrix computation complete, start traceback
            if not self.traceback_started:
                self.start_traceback()

    def previous_step(self):
        """Go back one step in the algorithm"""
        if self.traceback_started:
            # We're in traceback mode
            if len(self.highlighted_path) > 0:
                # Remove last highlighted cell from path
                i, j = self.highlighted_path.pop()

                # Remove highlight
                self.canvas.delete(f"highlight_{i}_{j}")

                # Update progress
                self.progress_label.config(text=f"Traceback: {len(self.highlighted_path)}/{len(self.traceback_path)} cells")

                # If we've removed all traceback cells, go back to matrix computation
                if len(self.highlighted_path) == 0:
                    self.traceback_started = False
                    self.next_button.config(state="normal")
                    self.end_button.config(state="normal")
                    self.result_label.config(text="")

                return

            # No more traceback to undo, so revert to matrix computation
            self.traceback_started = False
            self.canvas.delete("highlight")
            self.highlight_cell(self.m, self.n)  # Highlight bottom right cell
            self.next_button.config(state="normal")
            self.end_button.config(state="normal")
            return

        # Not in traceback mode, so go back a computation step
        if self.current_step_index > 0:
            self.current_step_index -= 1

            # Clear highlights
            self.canvas.delete("highlight")

            if self.current_step_index > 0:
                # Show previous step
                step = self.computation_steps[self.current_step_index-1]
                i, j = step['i'], step['j']

                # Highlight the cell
                self.highlight_cell(i, j)

                # Update progress
                self.progress_label.config(text=f"Step {self.current_step_index}/{len(self.computation_steps)}: "
                                         f"Showing cell ({i}, {j})")
            else:
                # At the beginning
                self.progress_label.config(text="Matrix initialized. Ready to start computation.")
                self.prev_button.config(state="disabled")
                self.start_button.config(state="disabled")

    def start_traceback(self):
        """Start the traceback process to find optimal alignment"""
        # Get the alignment and traceback path
        self.align1, self.align2, self.traceback_path = get_traceback_needleman_wunsch(
            self.seq1, self.seq2, self.score
        )

        # Display alignment
        self.result_label.config(text=f"{self.align1}\n{self.align2}")

        # Start with an empty path and add cells as we go
        self.highlighted_path = []

        # Mark that we've started traceback
        self.traceback_started = True

        # Update progress
        self.progress_label.config(text=f"Starting traceback from cell ({self.m}, {self.n})")

        # Start by highlighting the bottom right cell
        self.show_traceback_step(0)

    def show_more_traceback(self):
        """Show the next cell in the traceback path"""
        next_index = len(self.highlighted_path)
        if next_index < len(self.traceback_path):
            self.show_traceback_step(next_index)

    def show_traceback_step(self, index):
        """Show a specific step in the traceback"""
        if index < len(self.traceback_path):
            # Get the cell to highlight from the path - in REVERSE order
            # This shows the path from end to start (bottom-right to top-left)
            path_index = len(self.traceback_path) - 1 - index
            i, j = self.traceback_path[path_index]

            # Add to our highlighted path
            self.highlighted_path.append((i, j))

            # Highlight this cell - green for optimal path
            self.highlight_cell(i, j, "#ABEBC6")  # Light green

            # Update progress
            self.progress_label.config(text=f"Traceback: {len(self.highlighted_path)}/{len(self.traceback_path)} cells")

            # Disable next/end buttons when we're done
            if len(self.highlighted_path) == len(self.traceback_path):
                self.next_button.config(state="disabled")
                self.end_button.config(state="disabled")

    def go_to_start(self):
        """Go back to the beginning of the algorithm"""
        # Stop any ongoing animation
        if self.animation_in_progress:
            self.animation_in_progress = False
            self.canvas.after_cancel(self.animation_id)

        # Clear any highlights
        self.canvas.delete("highlight")

        # Reset state
        self.current_step_index = 0
        self.traceback_started = False
        self.highlighted_path = []
        self.animation_completed_matrix = False

        # Redraw the matrix
        self.create_matrix_visualization()

        # Update UI
        self.progress_label.config(text="Matrix initialized. Ready to start computation.")
        self.result_label.config(text="")

        # Update button states
        self.prev_button.config(state="disabled")
        self.next_button.config(state="normal")
        self.start_button.config(state="disabled")
        self.end_button.config(state="normal")

    def animate_to_end(self):
        """Animate all steps one by one until the end"""
        # Start animation if not already in progress
        if not hasattr(self, 'animation_in_progress') or not self.animation_in_progress:
            self.animation_in_progress = True
            self.animate_next_step()

    def animate_next_step(self):
        """Animate a single step and schedule the next one"""
        if not self.animation_in_progress:
            return

        # Determine if we have more steps
        has_more_steps = self.has_more_steps()

        if has_more_steps:
            # Perform the next step
            self.next_step()

            # Schedule the next animation step
            delay = self.get_animation_delay()
            self.animation_id = self.canvas.after(delay, self.animate_next_step)
        else:
            # Animation complete
            self.animation_in_progress = False

    def has_more_steps(self):
        """Check if there are more steps to animate"""
        if self.traceback_started:
            # Check if we still have traceback steps
            return len(self.highlighted_path) < len(self.traceback_path)
        elif self.animation_completed_matrix and not self.traceback_started:
            # Matrix is complete but traceback not started yet
            return True  # Need to start traceback
        else:
            # Check if we still have matrix computation steps
            return self.current_step_index < len(self.computation_steps)

    def get_animation_delay(self):
        """Get the delay between animation steps based on speed setting"""
        speed = self.speed_var.get()
        if speed == 1:  # Fast
            return 100
        elif speed == 2:  # Medium
            return 300
        else:  # Slow
            return 600
