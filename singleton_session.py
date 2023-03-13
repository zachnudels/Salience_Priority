from typing import List, Tuple
import os 

import numpy as np
from copy import deepcopy
from psychopy import visual
from psychopy import prefs
from psychopy import sound
prefs.hardware['audioLib'] = ['PTB']
from psychopy.tools.monitorunittools import deg2pix

from exptools2.core import PylinkEyetrackerSession
from singleton_trial import SingletonTrial
from utils import grid_coordinates, dva_per_pix, draw_instructions, construct_singleton_pairs

Coordinate = Tuple[int, int]


class SingletonSession(PylinkEyetrackerSession):
    def __init__(self, output_str, output_dir, settings_file, behav_file, eyetracker_on, subject_number, exp_num, debug):
        super().__init__(output_str=output_str, output_dir=output_dir, settings_file=settings_file,
                         eyetracker_on=eyetracker_on)

        self.results = []
        self.practice_trials = []
        self.trials = []
        self.trial_parameters = dict()
        self.trial_parameters["subject_number"] = subject_number
        self.RTs = []
        self.to_target_list = []
        self.exp_num = exp_num
        self.debug = debug
        self.screen = np.array([self.win.size[0], self.win.size[1]])

        if subject_number % 2 == 0:
            self.trial_parameters["target_orientation"] = self.settings["stimuli"]["singleton_orientation"][0]
            self.trial_parameters["distractor_orientation"] = self.settings["stimuli"]["singleton_orientation"][1]

        else:
            self.trial_parameters["target_orientation"] = self.settings["stimuli"]["singleton_orientation"][1]
            self.trial_parameters["distractor_orientation"] = self.settings["stimuli"]["singleton_orientation"][0]

        self.exp_str = f"exp{self.exp_num}"
        self.num_reps = self.settings["study"][self.exp_str]["num_reps"]
        self.bg_orientations = self.settings["stimuli"]["bg_orientation"]
        self.SOAs = [SOA*1e-3 for SOA in self.settings["study"][self.exp_str]["SOAs"]]
        print(self.SOAs)
        self.num_blocks = self.settings["study"][self.exp_str]["num_blocks"]

        
        self.behav_file = behav_file

        print(self.trial_parameters["target_orientation"], self.trial_parameters["distractor_orientation"])


        self.instructions = {
            "target_direction": ("left" if self.trial_parameters["target_orientation"] < 0 else "right"),
            "target_symbol": ("\\" if self.trial_parameters["target_orientation"] < 0 else "/"),
            "distractor_direction": ("left" if self.trial_parameters["distractor_orientation"] < 0 else "right"),
            "distractor_symbol": ("\\" if self.trial_parameters["distractor_orientation"] < 0 else "/"),
        }

                # TODO: Fixation point
        line_size_pix = deg2pix(self.settings["stimuli"]["fix_line_size_deg"], self.win.monitor)
        
        line_width_pix = deg2pix(self.settings["stimuli"]["fix_line_width_deg"], self.win.monitor)
        
        self.fixation = visual.ShapeStim(self.win,
                                         vertices=((0, -line_size_pix / 2), (0, line_size_pix / 2),
                                                   (0, 0),
                                                   (-line_size_pix / 2, 0), (line_size_pix / 2, 0)),
                                         lineWidth=line_width_pix,
                                         closeShape=False,
                                         lineColor=self.settings['stimuli']['fix_color'])
        
        self.max_dist = deg2pix(self.settings["stimuli"]["fix_check_rad"], self.win.monitor)
        
        ########################

    def run(self):
        self.create_trials()

        # if eyetracking then calibrate
        if self.eyetracker_on:
            self.calibrate_eyetracker()


        self.start_experiment()

        # draw instructions wait a few seconds
        this_instruction_string = (f"Instructions"
                                   f"\n\nThe task is to move your eyes to a line that is rotated to the "
                                   f"{self.instructions['target_direction']} ({self.instructions['target_symbol']})."
                                   f"\nTry not to make eye movements to the distractor line that is rotated to the "
                                   f"{self.instructions['distractor_direction']} "
                                   f"({self.instructions['distractor_symbol']})."
                                   f"\n Before each trial, press the -spacebar- to start"
                                   f"\n\n Press the -spacebar- to continue.")
        
        draw_instructions(self.win, this_instruction_string, keys='space')
        

        this_instruction_string = (f"\nTry to be as fast and accurate as possible!"
                                   f"\nMove as soon as you see the target (which is also when the fixation dot in the "
                                   f"middle disappears)"
                                   f"\n\n\nPlease press the -spacebar- to continue")

        draw_instructions(self.win, this_instruction_string, keys='space')

        this_instruction_string = (f"You will now start with a practice session."
                                   f"\n\nYou will get a warning if you select the distractor instead of the target"
                                   f"\nYou will also get a warning if you make more than one eye movement to reach the "
                                   f"target"
                                   f"\n\nPlease press the -spacebar- to begin.")

        draw_instructions(self.win, this_instruction_string, keys='space')

        if not self.debug:
            for trial in self.practice_trials:
                trial.run()

            this_instruction_string = (f"You will now start with the actual experiment."
                                    f"\n\nGoodluck."
                                    f"\n\nPlease press the -spacebar- to begin.")

            draw_instructions(self.win, this_instruction_string, keys='space')

            self.fixation.draw()
            self.win.flip()

        if self.eyetracker_on:
            self.start_recording_eyetracker()

        for trial in self.trials:
            trial.run()

        self.close()

    def create_trials(self):
        # TODO: Check this - I think Ines used vertical and psychopy uses horizontal to convert
        pixel_l = deg2pix(self.settings["stimuli"]["line_length"], self.win.monitor)

        pixel_w = deg2pix(self.settings["stimuli"]["line_width"], self.win.monitor)
        pixel_spacing = deg2pix(self.settings["stimuli"]["spacing"], self.win.monitor) + np.max([pixel_l, pixel_w])

        coordinates = grid_coordinates(self.settings["stimuli"]["x_count"],
                                       self.settings["stimuli"]["y_count"],
                                       pixel_spacing)

        possible_singleton_locations = construct_singleton_pairs(pixel_spacing,
                                                                 self.settings["stimuli"]["x_count"],
                                                                 self.settings["stimuli"]["y_count"])


        # # init variables for feedback
        #
        # check if they moved their eyes before targets are shown
        fixation_circle = visual.Circle(self.win, radius=50, pos=(0, 0))

        # Save memory - don't construct for each rep
        possible_grids = dict()
        possible_singletons = dict()
        target_circles = dict()
        dist_circles = dict()
        target_big_circles = dict()
        dist_big_circles = dict()

        for bg_orientation in self.bg_orientations:
            for target_coord, distractor_coord in possible_singleton_locations:
                grid_coords = [coordinate for coordinate in coordinates
                               if coordinate != target_coord and coordinate != distractor_coord]
                key = f"{bg_orientation}{target_coord}{distractor_coord}"

                possible_grids[key] = visual.ElementArrayStim(self.win,
                                                              nElements=len(grid_coords),
                                                              xys=grid_coords,
                                                              oris=bg_orientation,
                                                              units='pix',
                                                              autoLog=False,
                                                              sizes=(pixel_w, pixel_l),
                                                              elementMask=None,
                                                              elementTex=None,
                                                              interpolate=False)
                possible_singletons[key] = visual.ElementArrayStim(self.win,
                                                                   nElements=2,
                                                                   xys=[target_coord, distractor_coord],
                                                                   oris=[self.trial_parameters["target_orientation"],
                                                                         self.trial_parameters[
                                                                             "distractor_orientation"]],
                                                                   units='pix',
                                                                   autoLog=False,
                                                                   sizes=(pixel_w, pixel_l),
                                                                   elementMask=None,
                                                                   elementTex=None)
                target_circles[key] = visual.Circle(self.win, radius=50, pos=target_coord, lineColor='green')
                dist_circles[key] = visual.Circle(self.win, radius=50, pos=distractor_coord, lineColor='blue')
                target_big_circles[key] = visual.Circle(self.win, radius=107, pos=target_coord, lineColor='green')
                dist_big_circles[key] = visual.Circle(self.win, radius=107, pos=distractor_coord, lineColor='blue')

        num_trials = (self.num_reps *
                      len(self.bg_orientations) *
                      len(self.SOAs) *
                      len(possible_singleton_locations))
        
        print(num_trials)
        trials_per_block = num_trials / self.num_blocks
        trial_indices = np.arange(num_trials)
        # np.random.shuffle(trial_indices)  # Randomise trials
        trial_i = 0

        tone = sound.Sound('A')




        for rep in range(self.num_reps):
            for bg_orientation in self.bg_orientations:
                for SOA in self.SOAs:
                    for target_coord, distractor_coord in possible_singleton_locations:

                        trial_num = trial_indices[trial_i]
                        block_num = (trial_num // trials_per_block) + 1
                        trial_i += 1
                        key = f"{bg_orientation}{target_coord}{distractor_coord}"

                        _trial_parameters = deepcopy(self.trial_parameters)

                        # target salient or  not
                        if abs(self.trial_parameters["target_orientation"] - bg_orientation) > \
                                abs(self.trial_parameters["distractor_orientation"] - bg_orientation):
                            _trial_parameters["target_salience"] = "high"
                        else:
                            _trial_parameters["target_salience"] = "low"

                        _trial_parameters["target_pos"] = target_coord
                        _trial_parameters["distractor_pos"] = distractor_coord
                        _trial_parameters["bg_orientation"] = bg_orientation
                        _trial_parameters["SOA"] = SOA

                        _trial_parameters["practice"] = False

                        grid = possible_grids[key]
                        singletons = possible_singletons[key]

                        initialization_time = self.settings["trial_info"]["initialization_time"] 

                        if SOA < 0:
                            stimulus1 = grid
                            _trial_parameters["stimulus1_log"] = "stimuli_show"
                            stimulus2 = singletons
                            _trial_parameters["stimulus2_log"] = "target_display"
                            initialization_time += SOA # Adding a negative SOA will decrease time to ensure singletons presented after same duration
                        else:
                            stimulus1 = singletons
                            _trial_parameters["stimulus1_log"] = "target_display"
                            stimulus2 = grid
                            _trial_parameters["stimulus2_log"] = "stimuli_show"

                        phases = {
                            "initialization": initialization_time,
                            "stimulus1": np.abs(SOA),
                            "stimulus2": self.settings["trial_info"]["max_response_time"],
                            "ITI": self.settings["trial_info"]["ITI"] + np.random.uniform(-0.1, 0.1),
                        }

                        if trial_num % trials_per_block == 0:
                            phases["end_of_block"] = 1000

                        self.trials.append(SingletonTrial(self,
                                                          trial_nr=trial_num,
                                                          block_num=block_num,
                                                          phase_durations=tuple(
                                                              [float(duration) for duration in phases.values()]
                                                          ),
                                                          phase_names=tuple(phases.keys()),
                                                          parameters=_trial_parameters,
                                                          stimulus1=stimulus1,
                                                          stimulus2=stimulus2,
                                                          tone=tone,
                                                          debug=self.debug,
                                                          behavioural_file=self.behav_file,
                                                          fixation_circle=fixation_circle,
                                                          target_circle=target_circles[key],
                                                          distractor_circle=dist_circles[key],
                                                          target_circle_big=target_big_circles[key],
                                                          distractor_circle_big=dist_big_circles[key]))

        if not self.debug:
            self.practice_trials = np.random.choice(self.trials, self.settings["study"][self.exp_str]["practice_trials"])
            for i in range(len(self.practice_trials)):
                self.practice_trials[i].trial_num = i
                self.practice_trials[i].parameters["practice"] = True

                if i == len(self.practice_trials) - 1:
                    phases["end_of_block"] = 1000

            for i in range(len(self.practice_trials[0].phase_durations)):
                print(self.practice_trials[0].phase_names[i], self.practice_trials[0].phase_durations[i])

            # TODO: Change what makes it a practice: longer time somewhere?
            # practice_trial.phase_duration

        self.results = np.matrix(np.zeros([len(self.trials), 12]), dtype=object)


