import numpy as np

from exptools2.core import Trial, Session
import pylink
from typing import Dict, List, Optional, Tuple

import pickle
import psychtoolbox as ptb
from psychopy.visual.elementarray import ElementArrayStim
from psychopy.visual.circle import Circle
from psychopy import event, visual
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
                 tone: sound.Sound,
                 behavioural_file: str,
                 debug: bool,
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

        # Delete tuples from parameters since they cannot be logged
        del parameters["target_pos"]
        del parameters["distractor_pos"]


        self.stimulus1 = stimulus1
        self.stimulus2 = stimulus2

        self.trial_nr = trial_nr
        self.block_num = block_num

        self.parameters = parameters
        self.practice = self.parameters["practice"]

        self.a_sound = tone

        self.behavioural_file = behavioural_file
        self.debug=debug

        # self.singletons = singletons
        # self.target_stim = target_stim
        # self.distractor_stim = distractor_stim

    def draw(self):
        # print(self.phase)
        # for i in range(len(self.phase_names)):
        #     print(self.phase_names[i], self.phase_durations[i])
        # print('phase', self.phase_names[int(self.phase)])
        # print(self.parameters)
        # TODO: CHECK + What happens at block end??
        # check to see if it's time for a break

        if self.practice and self.trial_nr == self.session.settings["practice_trial_number"]:
            self.end_of_practise()

        if self.phase_names[int(self.phase)] == 'initialization':
            # self.session.tracker deals with everything here
            driftCheck = False
            if not self.debug:
                while not driftCheck:
                    print("Drift check")
                    driftCheck = self.session.tracker.doDriftCorrect(
                        self.session.settings["window_extra"]["size"][0] // 2,
                        self.session.settings["window_extra"]["size"][1] // 2,
                        1,
                        1)
                    # self.session.win.flip()
                    print("drift done")

            # self.session.tracker.sendMessage(
            #     f"trial {str(self.trial_nr)}_{str(round(np.nanmean(self.session.RTs) * 1000))}_{str(np.round(np.nanmean(self.to_target_list) * 100))}")

            # self.session.win.flip()
            self.session.tracker.sendMessage("fix_display")  # Logging to EDF

            # logging vars to edf file
            # TODO: to ensure logging works - CHECK everything is being logged. Maybe can make one long string with returns
            # TODO: in it. Check with Elle if this would still work for Analysis
            # Maybe check docs to see how many log calls can be made sequentially
            initialization_log_string = f"start_trial" \
                                        f"\nTRIAL_VAR trial_nr {self.trial_nr}" \
                                        f"\nTRIAL_VAR target_pos {self.target_pos}" \
                                        f"\nTRIAL_VAR dist_pos {self.distractor_pos}" \
                                        f"\nTRIAL_VAR target_co_x {self.target_pos[0]}" \
                                        f"\nTRIAL_VAR target_co_y {self.target_pos[1]}" \
                                        f"\nTRIAL_VAR distractor_co_x {self.distractor_pos[0]}" \
                                        f"\nTRIAL_VAR distractor_co_y {self.distractor_pos[1]}" \
                                        f"\nTRIAL_VAR background_orientation {self.bg_orientation}" \
                                        f"\nTRIAL_VAR ISI {self.SOA}" \
                                        f"\nTRIAL_VAR practice {self.parameters['practice']}" \
                                        f"\nTRIAL_VAR target_salience {self.parameters['target_salience']}"

            self.session.tracker.sendMessage(initialization_log_string)

            # print("End initialization")
            self.session.fixation.draw()


            # session.tracker.sendMessage("TRIAL_VAR target_salience %s" % target_salience)
            # exp.sleep(2)

        elif self.phase_names[int(self.phase)] == 'stimulus1':
            self.stimulus1.draw()
            # print("Draw stim1")
            self.session.tracker.sendMessage(self.parameters["stimulus1_log"])
            self.t0 = self.session.clock.getTime()  # TODO: Check if t0 needs to be accurate (target displayed rather than either)

        elif self.phase_names[int(self.phase)] == 'stimulus2':
            self.stimulus1.draw()
            self.stimulus2.draw()
            self.session.tracker.sendMessage(self.parameters["stimulus2_log"])
            eye_ev = self.session.tracker.getNextData()
            got_sac = False
            if (eye_ev is not None) and (eye_ev == pylink.ENDSACC):
                eye_dat = self.session.tracker.getFloatData()
                self.endtime = eye_dat.getEndTime()  # offset time
                self.startpos = eye_dat.getStartGaze()  # start position
                self.endpos = eye_dat.getEndGaze()  # end position
                got_sac = True

            if got_sac:
                self.response, self.RT = self.response_check(self.t0, self.practice)  # TODO: Practice stuff             
                self.session.RTs.append(self.RT)
                self.stop_phase()

        elif self.phase_names[int(self.phase)] == 'ITI':
            # now check if the eye movement was made correctly
            

            self.session.fixation.draw()
            self.success = self.check_saccade(self.endpos)
            self.session.tracker.sendMessage('TRIAL_VAR end_position ' + str(self.endpos))
            self.session.tracker.sendMessage('TRIAL_VAR succes ' + str(self.success))

            self.session.tracker.sendMessage('stop_trial')
            # self.session.tracker.stop_recording()
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

        elif self.phase_names[int(self.phase)] == 'end_of_block':
            avg_RT = np.round(np.nanmean(self.session.RTs) * 1000)
            this_instruction_string = (f"This was block {self.block_num}"
                                   f"\n\n Your average response time was {avg_RT} ms."
                                   f"\n Well done! "
                                   f"\n\nPlease press the -spacebar- to begin.")


            text1 = visual.TextStim(self.session.win, text=this_instruction_string, pos=(0, 0), font="Helvetica Neue", height=30,
                                    anchorHoriz='center', anchorVert='center')

            text1.draw()
            self.session.RTs = []
            self.save_results()

        # print("Draw fixaxtion")
        # utils.draw_instructions(self.session.win, "TEST", keys='space')
        
        # print(self.parameters)

        # saving the results matrix
        # print(blue)

        


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
        with open(self.behavioural_file, 'wb') as handle:
            pickle.dump(self.session.results, handle, protocol=pickle.HIGHEST_PROTOCOL)



    def check_saccade(self, endpos):
        """
        Check whether movement is one continuous saccade, not broken up
        Parameters
        ----------
        endpos

        Returns
        -------

        """
        if endpos is None:
            return False
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
                    self.session.close()
                    self.session.quit()

                elif ev in['space'] and self.phase_names[int(self.phase)] == "end_of_block":
                    self.stop_phase()

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
