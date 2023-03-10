import numpy as np

from exptools2.core import Trial, Session
from typing import Dict, Optional, Tuple

import pickle
import psychtoolbox as ptb
from psychopy.visual.elementarray import ElementArrayStim
from psychopy.visual.circle import Circle
from psychopy import event, visual
from psychopy import prefs
prefs.hardware['audioLib'] = ['PTB']
from psychopy import sound

import utils


# Opening edf files:
# open powershell in the folder
# command: edf2asc <filename> -e(events)/s(samples)


class SingletonTrial(Trial):
    def __init__(self,
                 session: Session,
                 trial_nr: int,
                 block_num: int,
                 phase_durations: Tuple[float, ...],
                 phase_names: Tuple[str],
                 parameters: Dict,
                 stimulus1: ElementArrayStim,
                 stimulus2: ElementArrayStim,
                 fixation_circle: Circle,
                 target_circle: Circle,
                 distractor_circle: Circle,
                 target_circle_big: Circle,
                 distractor_circle_big: Circle):

        super().__init__(session, trial_nr, phase_durations, phase_names, verbose=False)

        self.RT = None
        self.success = None
        self.endpos = None
        self.startpos = None
        self.endtime = None
        self.t0 = None
        self.response = None
        self.target_pos = parameters['target_pos']
        self.distractor_pos = parameters['distractor_pos']
        self.SOA = parameters['SOA']
        self.bg_orientation = parameters['bg_orientation']
        self.target_orientation = parameters["target_orientation"]
        self.distractor_orientation = parameters["distractor_orientation"]

        self.stimulus1 = stimulus1
        self.stimulus2 = stimulus2

        self.trial_nr = trial_nr
        self.block_num = block_num

        self.parameters = parameters
        self.practice = self.parameters["practice"]

        self.a_sound = sound.Sound('A')

        # self.singletons = singletons
        # self.target_stim = target_stim
        # self.distractor_stim = distractor_stim

    def draw(self):
        # TODO: CHECK + What happens at block end??
        # check to see if it's time for a break

        if self.practice and self.trial_nr == self.session.settings["practice_trial_number"]:
            self.end_of_practise()

        if self.phase_names[int(self.phase)] == 'initialization':
            # self.session.tracker deals with everything here
            driftCheck = False
            while not driftCheck:
                driftCheck = self.session.tracker.drift_correction()
                self.session.win.flip()

            self.session.tracker.start_recording()
            self.session.tracker.status_msg(
                f"trial {str(self.trial_nr)}_{str(round(np.nanmean(self.session.RTs) * 1000))}_{str(np.round(np.nanmean(self.to_target_list) * 100))}")

            self.session.fixation.draw()
            self.session.win.flip()
            self.session.tracker.log("fix_display")  # Logging to EDF

            # logging vars to edf file
            # TODO: to ensure logging works - CHECK everything is being logged. Maybe can make one long string with returns
            # TODO: in it. Check with Elle if this would still work for Analysis
            # Maybe check docs to see how many log calls can be made sequentially
            initialization_log_string = f"start_trial" \
                                        f"\nTRIAL_VAR trial_nr {self.trial_nr}" \
                                        f"\nTRIAL_VAR target_pos {self.target_pos}" \
                                        f"\nTRIAL_VAR dist_pos {self.distractor_pos}" \
                                        f"\nTRIAL_VAR target_co_x {self.parameters['target_pos'][0]}" \
                                        f"\nTRIAL_VAR target_co_y {self.parameters['target_pos'][1]}" \
                                        f"\nTRIAL_VAR distractor_co_x {self.parameters['distractor_pos'][0]}" \
                                        f"\nTRIAL_VAR distractor_co_y {self.parameters['distractor_pos'][1]}" \
                                        f"\nTRIAL_VAR background_orientation {self.bg_orientation}" \
                                        f"\nTRIAL_VAR ISI {self.SOA}" \
                                        f"\nTRIAL_VAR practice {self.parameters['practice']}" \
                                        f"\nTRIAL_VAR target_salience {self.parameters['target_salience']}"

            self.session.tracker.log(initialization_log_string)

            # session.tracker.log("TRIAL_VAR target_salience %s" % target_salience)
            # exp.sleep(2)

        elif self.phase_names[int(self.phase)] == 'stimulus1':
            self.stimulus1.draw()
            self.session.eyetracker.log(self.parameters["stimulus1_log"])
            self.t0 = self.session.clock.getTime()  # TODO: Check if t0 needs to be accurate (target displayed rather than either)

        elif self.phase_names[int(self.phase)] == 'stimulus2':
            self.stimulus2.draw()
            self.session.eyetracker.log(self.parameters["stimulus2_log"])
            self.endtime, self.startpos, self.endpos = self.session.tracker.wait_for_saccade_end()

            self.response, self.RT = self.response_check(self.t0, self.practice)  # TODO: Practice stuff
            if self.response is not None:
                self.stop_phase()

            self.session.RTs.append(self.RT)

        elif self.phase_names[int(self.phase)] == 'ITI':
            # now check if the eye movement was made correctly
            self.success = self.check_saccade(self.endpos)
            self.session.tracker.log('TRIAL_VAR end_position ' + str(self.endpos))
            self.session.tracker.log('TRIAL_VAR succes ' + str(self.success))

            self.session.tracker.log('stop_trial')
            self.session.tracker.stop_recording()
            self.session.results[self.trial_nr, :] = [self.parameters["subject_number"],
                                                      self.trial_nr,
                                                      self.parameters["target_orientation"],
                                                      self.parameters["distractor_orientation"],
                                                      self.bg_orientation,
                                                      self.target_pos,
                                                      self.distractor_pos,
                                                      self.SOA,
                                                      self.success,
                                                      self.endpos,
                                                      self.parameters["target_salience"],
                                                      self.RT]

            self.save_results()

            # TODO: Check this
            # get current gaze
            curr_gaze = utils.getCurSamp(self.session.tracker, screen=self.session.screen)
            # calculate distance to center of screen
            dist2center = utils.distBetweenPoints(curr_gaze, (0, 0))

            # Check if sample is within boundary
            if dist2center < self.session.maxDist:
                self.session.gaze_sampleCount += 1

            # If enough samples within boundary
            if self.session.gaze_sampleCount >= self.session.settings['visual_search']['gaze_sampleCount']:
                # print('correctFixation')
                self.session.gaze_sampleCount = 0
                self.stop_phase()

        elif self.phase_names[int(self.phase)] == 'end_of_block':
            avg_RT = np.round(np.nanmean(self.session.RTs) * 1000)
            self.give_feedback(avg_RT)
            self.session.RTs = []

        # saving the results matrix

        self.save_results()

    def get_events(self):
        """ Logs responses/triggers """
        for ev, t in event.getKeys(
                timeStamped=self.session.clock):  # list of of (keyname, time) relative to Clock’s last reset
            if len(ev) > 0:
                if ev in ['q']:
                    print('trial canceled by user')
                    self.session.close()
                    self.session.quit()

    def response_check(self, t0, practice: bool):
        """
        Can just check if response was correct
        Parameters
        ----------
        practice
        t0

        Returns
        -------

        """
        warning_text = visual.TextStim(self.session.win, text='Wrong line!')

        curPos = self.session.tracker.sample()
        x_sample = curPos[0] - self.session.settings["self.session.windows_extra"]["size"][0] / 2
        y_sample = (curPos[1] - (self.session.settings["self.session.windows_extra"]["size"][1] / 2)) * -1
        curPos = (x_sample, y_sample)
        t1 = self.session.clock.getTime()
        RT = t1 - t0

        if self.target_circle.contains(curPos):
            response = 'target'
            self.to_target_list.append(1)
        elif self.dist_circle.contains(curPos):
            response = 'distractor'
            self.to_target_list.append(0)
            self.mySound.play(when=ptb.GetSecs())
            if practice:
                warning_text.draw()
                self.session.win.flip()
                # core.wait(1) TODO
        else:
            response = 'none'
        return response, RT

    def draw_fixation(self, pos, lineSize, draw_cross=True, color='white'):
        t = lineSize / 2.0
        vertical_line = visual.Line(self.session.win, start=(pos[0], pos[1] - t), end=(pos[0], pos[1] + t),
                                    lineColor=color)
        horizontal_line = visual.Line(self.session.win, start=(pos[0] - t, pos[1]), end=(pos[0] + t, pos[1]),
                                      lineColor=color)

        if draw_cross:
            vertical_line.draw()
            horizontal_line.draw()

    def save_results(self):
        # change this based on folder name!
        filename = f"C:/Users/RA-Eyelink/Data/behavioral_data_mieke{self.parameters['subject_number']}.pickle"
        with open(filename, 'wb') as handle:
            pickle.dump(self.session.results, handle, protocol=pickle.HIGHEST_PROTOCOL)
        filename = f"//labssrv/Labs/RA-Eyelink/Data/behavioral_data_mieke{self.parameters['subject_number']}.pickle"
        with open(filename, 'wb') as handle:
            pickle.dump(self.session.results, handle, protocol=pickle.HIGHEST_PROTOCOL)
        print('Saved the pickle file!')

    def give_feedback(self, avg_RT):
        text1 = visual.TextStim(self.session.win, text='This was block ' + str(self.block_num), pos=(0, 200))
        text2 = visual.TextStim(self.session.win, text='Your average response time was ' + str(avg_RT) + ' ms',
                                pos=(0, -0))
        text3 = visual.TextStim(self.session.win, text='Press c to continue', pos=(0, -200))

        text1.draw()
        text2.draw()
        text3.draw()
        self.session.win.flip()

        event.waitKeys(keyList=['c'])

    def check_saccade(self, endpos):
        """
        Check whether movement is one continuous saccade, not broken up
        Parameters
        ----------
        endpos

        Returns
        -------

        """
        x_sample = endpos[0] - self.settings['window_extra']['size'][0] / 2
        y_sample = (endpos[1] - (self.settings['window_extra']['size'][1] / 2)) * -1
        curPos = (x_sample, y_sample)

        # TODO: Use distance instead of circle
        if self.target_circle_big.contains(curPos) | self.dist_circle_big.contains(curPos):
            success = True
        else:
            success = False

        if not success:
            warning_text = visual.TextStim(self.session.win, text='Make one eye movement to the target!')
            warning_text.draw()
            self.session.win.flip()
        return success

    def get_events(self):
        """ Logs responses/triggers """
        for ev, t in event.getKeys(
                timeStamped=self.session.clock):  # list of of (keyname, time) relative to Clock’s last reset
            if len(ev) > 0:
                if ev in ['q']:
                    print('trial canceled by user')
                    self.session.close_all()
                    self.session.quit()

#
# class PracticeSingletonTrial(SingletonTrial):
#     def __init__(self,
#                  session: Session,
#                  trial_nr: int,
#                  phase_durations: Tuple[int],
#                  phase_names: Optional[Tuple[str]],
#                  parameters: Optional[Dict]):
#         super().__init__(session, trial_nr, phase_durations, phase_names, parameters)
#
#         """
#         if trial_nr == session.settings.nr_practice_trials:
#             end_of_practice()
#         """
#
#         def end_of_practise():
#             text1 = visual.TextStim(self.session.win, text='This is the end of the practise block', pos=(0, 200))
#             text2 = visual.TextStim(self.session.win,
#                                     text='From now on you will only hear a beep if you select the distractor (no longer the text Wrong Line!)',
#                                     pos=(0, -0))
#             text3 = visual.TextStim(self.session.win, text='Press c to continue', pos=(0, -200))
#
#             text1.draw()
#             text2.draw()
#             text3.draw()
#             self.session.win.flip()
#
#             event.waitKeys(keyList=['c'])
#             # core.wait(0.5)
