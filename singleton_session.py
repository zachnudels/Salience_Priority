from typing import List, Tuple

import numpy as np
from copy import deepcopy
from psychopy import visual
from psychopy.tools.monitorunittools import deg2pix

from session import ExpSession
from singleton_trial import SingletonTrial, PracticeSingletonTrial
from utils import grid_coordinates, dva_per_pix

Coordinate = Tuple[int, int]


class SingletonSession(ExpSession):
    def __init__(self, output_str, output_dir, settings_file, eyetracker_on, subject_number, exp_num):
        super().__init__(output_str=output_str, output_dir=output_dir, settings_file=settings_file,
                         eyetracker_on=eyetracker_on)

        trial_parameters = dict()
        trial_parameters["subject_number"] = subject_number
        self.RTs = []

        if subject_number % 2 == 0:
            trial_parameters["target_orientation"] = self.settings["stimuli"]["targetOrientation"][0]
            trial_parameters["distractor_orientation"] = self.settings["stimuli"]["distractorOrientation"][0]

        else:
            trial_parameters["target_orientation"] = self.settings["stimuli"]["targetOrientation"][1]
            trial_parameters["distractor_orientation"] = self.settings["stimuli"]["distractorOrientation"][1]

        # target salient or  not
        if abs(trial_parameters["target_orientation"] - trial_parameters["bg_orientation"]) > \
                abs(trial_parameters["distractor_orientation"] - trial_parameters["bg_orientation"]):
            trial_parameters["target_salience"] = "high"
        else:
            trial_parameters["target_salience"] = "low"

        # TODO: when making instruction screen, use sign of target to help with [target] and [target_symbol]
        # TODO: As with distractor

        # TODO: Check this - I think Ines used vertical and psychopy uses horizontal to convert
        pixel_l = deg2pix(self.settings["stimuli"]["line_length"], self.win.monitor)

        pixel_w = deg2pix(self.settings["stimuli"]["line_width"], self.win.monitor)
        pixel_spacing = deg2pix(self.setting["stimuli"]["spacing"], self.win.monitor) + np.max([pixel_l, pixel_w])

        # x_ and y_ spacing are spacing + line width and length, respectively
        coordinates = grid_coordinates(self.settings["stimuli"]["x_count"],
                                       self.settings["stimuli"]["y_count"],
                                       pixel_spacing)

        possible_singleton_locations = self.construct_singleton_pairs()

        line_size_pix = self.settings['stimuli']['fix_line_size_deg'] / utils.dva_per_pix(
            height_cm=self.settings['monitor_extra']['height'],
            distance_cm=self.settings['monitor']['distance'],
            vert_res_pix=self.screen[-1])

        line_width_pix = self.settings['stimuli']['fix_line_width_deg'] / utils.dva_per_pix(
            height_cm=self.settings['monitor_extra']['height'],
            distance_cm=self.settings['monitor']['distance'],
            vert_res_pix=self.screen[-1])

        self.fixation = visual.ShapeStim(self.win,
                                         vertices=((0, -line_size_pix / 2), (0, line_size_pix / 2),
                                                   (0, 0),
                                                   (-line_size_pix / 2, 0), (line_size_pix / 2, 0)),
                                         lineWidth=line_width_pix,
                                         closeShape=False,
                                         lineColor=self.settings['stimuli']['fix_color'])

        ########################
        # # init variables for feedback
        #
        # check if they moved their eyes before targets are shown
        fixation_circle = visual.Circle(self.win, radius=50, pos=(0, 0))

        self.trials = []

        # Save memory - don't construct for each rep
        possible_grids = dict()
        possible_singletons = dict()
        target_circles = dict()
        dist_circles = dict()
        target_big_circles = dict()
        dist_big_circles = dict()

        for bg_orientation in self.settings["stimuli"]["bg_orientation"]:
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
                                                                   nElements=1,
                                                                   xys=[target_coord, distractor_coord],
                                                                   oris=[trial_parameters["target_orientation"],
                                                                         trial_parameters["distractor_orientation"]],
                                                                   units='pix',
                                                                   autoLog=False,
                                                                   sizes=(pixel_w, pixel_l),
                                                                   elementMask=None,
                                                                   elementTex=None)
                target_circles[key] = visual.Circle(self.win, radius=50, pos=target_coord, lineColor='green')
                dist_circles[key] = visual.Circle(self.win, radius=50, pos=distractor_coord, lineColor='blue')
                target_big_circles[key] = visual.Circle(self.win, radius=107, pos=target_coord, lineColor='green')
                dist_big_circles[key] = visual.Circle(self.win, radius=107, pos=distractor_coord, lineColor='blue')

        # TODO: Check with Mieke is all of the below should be shuffled (I think yes but then what should num_reps be?)
        trial_num = 0
        for rep in range(self.settings["study"][f"exp{exp_num}"]["num_reps"]):
            for bg_orientation in self.settings["stimuli"]["bg_orientation"]:
                for SOA in self.settings["study"][f"exp{exp_num}"]:
                    key = f"{bg_orientation}{target_coord}{distractor_coord}"

                    for target_coord, distractor_coord in possible_singleton_locations:
                        _trial_parameters = deepcopy(trial_parameters)
                        _trial_parameters["target_pos"] = target_coord
                        _trial_parameters["dist_pos"] = distractor_coord
                        _trial_parameters["bg_orientation"] = bg_orientation
                        _trial_parameters["SOA"] = SOA * 1e-3

                        grid = possible_grids[key]
                        singletons = possible_singletons[key]

                        if SOA < 0:
                            stimulus1 = grid
                            _trial_parameters["stimulus1_log"] = "stimuli_show"
                            stimulus2 = singletons
                            _trial_parameters["stimulus2_log"] = "target_display"
                        else:
                            stimulus1 = singletons
                            _trial_parameters["stimulus1_log"] = "target_display"
                            stimulus2 = grid
                            _trial_parameters["stimulus2_log"] = "stimuli_show"


                        # TODO: Deal with blocks here

                        phases = {
                            "initialization": self.settings["trial_info"]["initialization_time"],
                            "stimulus1": np.abs(SOA),
                            "stimulus2": self.settings["trial_info"]["max_response_time"],
                            "ITI": self.settings["trial_info"]["ITI"],
                        }

                        self.trials.append(SingletonTrial(self,
                                                          trial_nr=trial_num,
                                                          phase_durations=tuple(
                                                              [float(duration) for duration in phases.values()]
                                                          ),
                                                          phase_names=tuple(phases.keys()),
                                                          parameters=trial_parameters,
                                                          stimulus1=stimulus1,
                                                          stimulus2=stimulus2,
                                                          fixation_circle=fixation_circle,
                                                          target_circle=target_circles[key],
                                                          distractor_circle=dist_circles[key],
                                                          target_circle_big=target_big_circles[key],
                                                          distractor_circle_big=dist_big_circles[key]))
                        trial_num += 1

        np.random.shuffle(self.trials)
        self.results = np.matrix(np.zeros([len(self.trials), 12]), dtype=object)
        self.to_target_list = []

    # TODO
    # TODO: Check with Mieke if we want to only have singletons on the diagonals (paper says 6 possible locations)
    def construct_singleton_pairs(self) -> List[Tuple[Coordinate, Coordinate]]:
        """
        Should return every possible combination of target, distractor pairs
        Returns
        -------

        """
        return [((-1, 1), (1, -1)), ((-1, 1), (1, -1))]
