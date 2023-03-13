import numpy as np
import os, sys
import os.path as op
import math
import random
import pandas as pd
import yaml

from psychopy import visual, tools, colors, event
import psychopy.tools.colorspacetools as ct
from psychopy.tools.monitorunittools import deg2pix
import itertools
from typing import Optional

import time
import colorsys
import seaborn as sns

import pylink

# TODO: This should return two sizes, one for vert, one for hori
def dva_per_pix(height_cm=30, distance_cm=70, vert_res_pix=1080):
    """ calculate degrees of visual angle per pixel,
    to use for screen boundaries when plotting/masking
    Parameters
    ----------
    height_cm : int
        screen height
    distance_cm: float
        screen distance (same unit as height)
    vert_res_pix : int
        vertical resolution of screen
    
    Outputs
    -------
    deg_per_px : float
        degree (dva) per pixel
    
    """

    # screen size in degrees / vertical resolution
    deg_per_px = (2.0 * np.degrees(np.arctan(height_cm / (2.0 * distance_cm)))) / vert_res_pix

    return deg_per_px


def circle_points(radius, n_points):
    """ define positions in circle

    Parameters
    ----------
    radius : list/arr
        list of radius
    n_points: list/arr
        number of points per radius
    
    Outputs
    -------
    circles : list
        list of [x,y] positions per radius
    
    """
    circles = []

    for r, n in zip(radius, n_points):
        t = np.arange(0, 2 * np.pi, 2 * np.pi / float(n))  # np.linspace(0, 2 * np.pi, n)
        x = r * np.cos(t)
        y = r * np.sin(t)
        circles.append(np.c_[x, y])

    return circles


def grid_coordinates(x_count: int, y_count: int, pixel_spacing: float):
    """
    Given the number of x, and y points and the desired spacing in degrees, return a list of pixel coordiantes
    representing a rectangular grid
    Parameters
    ----------
    x_count
    y_count
    pixel_spacing

    Returns
    -------

    """
    half_x = x_count // 2
    half_y = y_count // 2

    x_range = np.arange(-pixel_spacing * half_x, pixel_spacing * (1 + half_x), pixel_spacing, dtype=np.int32)
    y_range = np.arange(-pixel_spacing * half_y, pixel_spacing * (1 + half_y), pixel_spacing, dtype=np.int32)
    X, Y = np.meshgrid(x_range, y_range)
    return list(zip(X.ravel(), Y.ravel()))


def get_grid_array(positions, ecc_range, convert2pix=True, screen=[1920, 1080],
                   height_cm=30, distance_cm=70,
                   constraint_type='ellipse', constraint_bounds_pix=[500, 700]):
    """ get position array
    needs postion list with positions per ecc
    and ecc range

    Parameters
    ----------
    positions : list/arr
        list of [x,y] positions per ecc
    ecc_range: list/arr
        list with eccs in position
    convert2pix: bool
        if outputted list in pixels or not
    constrain_type: str
        type of position contraint to use eg: 'ellipse', 'square', 'rectangle'
    constraint_bounds_pix: list/arr
        bounds to constraint positions to
    
    Outputs
    -------
    pos_list : arr
        list of [x,y] positions (pix if convert2pix == True)
    ecc_list: arr
        list of ecc per position pair (dva)
    
    """

    pos_list = []
    ecc_list = []

    # if converting to pix, then need to convert the bounds to deg
    if convert2pix:
        constraint_bounds = constraint_bounds_pix * dva_per_pix(height_cm=height_cm,
                                                                distance_cm=distance_cm,
                                                                vert_res_pix=screen[-1])
    else:
        constraint_bounds = constraint_bounds_pix

    for ind, e in enumerate(positions):

        # append list of positions
        for pos in e:
            # check if within bounds
            if constraint_type == 'ellipse' and \
                    (((pos[0] ** 2) / (max(constraint_bounds) ** 2) + (pos[1] ** 2) / (
                            min(constraint_bounds) ** 2)) <= 1):
                pos_list.append(list(pos))

                # append eccentricity of these positions
                ecc_list.append(ecc_range[ind])

    if convert2pix:
        pos_list = pos_list / dva_per_pix(height_cm=height_cm,
                                          distance_cm=distance_cm,
                                          vert_res_pix=screen[-1])
    else:
        pos_list = np.array(pos_list)

    return pos_list, np.array(ecc_list)


def draw_instructions(win, instructions, keys=['b'], visual_obj=[], image_path=[],
                      color=(1, 1, 1), font='Helvetica Neue', pos=(0, 0), height=30,  # .65,
                      italic=True, anchorHoriz='center', anchorVert='center'):
    """ draw instructions on screen
    
    Parameters
    ----------
    win : object
        window object to draw on
    instructions : str
        instruction string to draw 
    key: list
        list of keys to skip instructions
    visual_obj: list
        if not empty, should have psychopy visual objects (to add to the display ex: side rectangles to limit display)
        
    """

    text = visual.TextStim(win=win,
                           text=instructions,
                           color=color,
                           font=font,
                           pos=pos,
                           height=height,
                           italic=italic,
                           anchorHoriz=anchorHoriz,
                           anchorVert=anchorVert
                           )

    # draw text 
    text.draw()

    if len(visual_obj) > 0:
        for w in range(len(visual_obj)):
            visual_obj[w].draw()

    if len(image_path) > 0:
        for _, img in enumerate(image_path):
            img_stim = visual.ImageStim(win=win,
                                        image=img,
                                        pos=(0, 100))
            img_stim.draw()

    win.flip()

    key_pressed = event.waitKeys(keyList=keys)

    return key_pressed


def rgb255_2_hsv(arr):
    """ convert RGB 255 to HSV
    
    Parameters
    ----------
    arr: list/array
        1D list of rgb values
        
    """

    rgb_norm = np.array(arr) / 255

    hsv_color = np.array(colorsys.rgb_to_hsv(rgb_norm[0], rgb_norm[1], rgb_norm[2]))
    hsv_color[0] = hsv_color[0] * 360

    return hsv_color


def near_power_of_2(x, near='previous'):
    """ Get nearest power of 2
    
    Parameters
    ----------
    x : int/float
        value for which we want to find the nearest power of 2
    near : str
        'previous' or 'next' to indicate if floor or ceiling power of 2        
    """
    if x == 0:
        val = 1
    else:
        if near == 'previous':
            val = 2 ** math.floor(math.log2(x))
        elif near == 'next':
            val = 2 ** math.ceil(math.log2(x))

    return val


def update_elements(ElementArrayStim, elem_positions=[], grid_pos=[],
                    elem_color=[204, 204, 204], elem_ori=[353, 7],
                    elem_sf=4, elem_names=['bR', 'bL', 'pL', 'pR'],
                    key_name=['bR', 'bL']):
    """ update element array settings
    
    Parameters
    ----------
    ElementArrayStim: Psychopy object
    	ElementArrayStim to be updated 
    condition_settings: dict
        dictionary with all condition settings
    this_phase: str
        string with name of condition to be displayed
    elem_positions: arr
         numpy array with element positions to be updated and shown (N,2) -> (number of positions, [x,y])
         to be used for opacity update
    grid_pos: arr
        numpy array with element positions (N,2) of whole grid -> (number of positions, [x,y])
    monitor: object
        monitor object (to get monitor references for deg2pix transformation)
    screen: arr
        array with display resolution
    luminance: float or None
        luminance increment to alter color (used for flicker task)
    update_settings: bool
        choose if we want to update settings or not (mainly for color changes)
    new_color: array
        if we are changing color to be one not represented in settings (ca also be False if no new color used)
        
    """

    # set number of elements
    if len(grid_pos) == 0:
        nElements = 1
    else:
        nElements = grid_pos.shape[0]

    ## to make colored gabor, need to do it a bit differently (psychopy forces colors to be opposite)
    # get rgb color and convert to hsv
    hsv_color = rgb255_2_hsv(elem_color)
    grat_res = near_power_of_2(ElementArrayStim.sizes[0][0],
                               near='previous')  # use power of 2 as grating res, to avoid error

    # initialise grating
    grating = visual.filters.makeGrating(res=grat_res)
    grating_norm = (grating - np.min(grating)) / (np.max(grating) - np.min(grating))  # normalize between 0 and 1

    # initialise a base texture 
    colored_grating = np.ones((grat_res, grat_res, 3))

    # replace the base texture red/green channel with the element color value, and the value channel with the grating
    colored_grating[..., 0] = hsv_color[0]
    colored_grating[..., 1] = hsv_color[1]
    colored_grating[..., 2] = grating_norm * hsv_color[2]

    elementTex = ct.hsv2rgb(colored_grating)  # convert back to rgb

    # update element colors to color of the patch 
    element_color = np.ones((int(np.round(nElements)), 3))

    # update element spatial frequency
    element_sfs = np.ones((nElements)) * elem_sf  # in cycles/gabor width

    # update element orientation
    element_ori = np.ones((nElements))

    if nElements > 1:
        # get left and right indices from keys names
        L_indices = [ind for ind, k in enumerate(elem_names) if k in key_name and 'L' in k]
        R_indices = [ind for ind, k in enumerate(elem_names) if k in key_name and 'R' in k]

        # make grid and element position lists of lists
        list_grid_pos = [list(val) for _, val in enumerate(grid_pos)]

        if len(L_indices) > 0:
            list_Lelem_pos = [list(val) for _, val in enumerate(elem_positions[L_indices])]
            # get left and right global indices (global, because indices given grid pos)
            L_glob_indices = [list_grid_pos.index(list_Lelem_pos[i]) for i in range(len(list_Lelem_pos))]

            element_ori[L_glob_indices] = np.array(elem_ori)[L_indices][0]
        else:
            L_glob_indices = []

        if len(R_indices) > 0:
            list_Relem_pos = [list(val) for _, val in enumerate(elem_positions[R_indices])]
            # get left and right global indices (global, because indices given grid pos)
            R_glob_indices = [list_grid_pos.index(list_Relem_pos[i]) for i in range(len(list_Relem_pos))]

            element_ori[R_glob_indices] = np.array(elem_ori)[R_indices][0]
        else:
            R_glob_indices = []

        # combine left and right global indices
        glob_indices = L_glob_indices + R_glob_indices

    else:
        glob_indices = 0
        element_ori[glob_indices] = elem_ori[glob_indices]
        element_positions = np.array([elem_positions])

    # set element contrasts
    element_contrast = np.zeros(nElements)
    element_contrast[glob_indices] = 1

    # set element opacities
    element_opacities = np.zeros(nElements)
    element_opacities[glob_indices] = 1

    # set all of the above settings
    ElementArrayStim.setTex(elementTex)
    ElementArrayStim.setSfs(element_sfs)
    ElementArrayStim.setOpacities(element_opacities)
    ElementArrayStim.setOris(element_ori)
    ElementArrayStim.setColors(element_color, 'rgb')
    ElementArrayStim.setContrs(element_contrast)
    if nElements == 1:
        ElementArrayStim.setXYs(element_positions)

    return (ElementArrayStim)


def update_grating(GratingStim,
                   elem_color=[204, 204, 204], elem_ori=7,
                   elem_sf=4,
                   elem_pos=(0, 0)):
    """ update grating stim settings
    
    Parameters
    ----------
    GratingStim: Psychopy object
    	GratingStim to be updated 
      
    """

    ## to make colored gabor, need to do it a bit differently (psychopy forces colors to be opposite)
    # get rgb color and convert to hsv
    hsv_color = rgb255_2_hsv(elem_color)
    grat_res = near_power_of_2(GratingStim.size[0], near='previous')  # use power of 2 as grating res, to avoid error

    # initialise grating
    grating = visual.filters.makeGrating(res=grat_res)
    grating_norm = (grating - np.min(grating)) / (np.max(grating) - np.min(grating))  # normalize between 0 and 1

    # initialise a base texture 
    colored_grating = np.ones((grat_res, grat_res, 3))

    # replace the base texture red/green channel with the element color value, and the value channel with the grating
    colored_grating[..., 0] = hsv_color[0]
    colored_grating[..., 1] = hsv_color[1]
    colored_grating[..., 2] = grating_norm * hsv_color[2]

    elementTex = ct.hsv2rgb(colored_grating)  # convert back to rgb

    # update element colors to color of the patch 
    element_color = np.ones((1, 3))

    # set all of the above settings
    GratingStim.tex = elementTex
    GratingStim.pos = elem_pos
    GratingStim.sf = elem_sf
    GratingStim.ori = elem_ori
    GratingStim.setColor(element_color, 'rgb')
    GratingStim.mask = 'gauss'

    return (GratingStim)


def getCurSamp(tracker, screen=[1920, 1080]):
    """
    Gets the most recent gaze position sample from the eyelink. This
    sample might be a couple of ms delayed, depending on the eyelink
    settings used.
    The eyetracker needs to be in recording mode for this to work.

   Parameters
    ----------
    tracker: Eyelink
    	Pylink 'Eyelink' object (tracker = pylink.Eyelink)

    Returns
    -------
    curSamp : tuple
        The (x,y) gaze position on the screen. In center-based coordinates.
    Examples
    --------
    >>> curSamp = tracker.getCurSamp()
    >>> curSamp
    (100,250)
   
    """
    curSamp = tracker.getNewestSample()

    if curSamp is not None:

        if curSamp.isRightSample():
            gazePos = curSamp.getRightEye().getGaze()

        if curSamp.isLeftSample():
            gazePos = curSamp.getLeftEye().getGaze()

        newGazePos = [0, 0]
        newGazePos[0] = gazePos[0] - screen[0] / 2
        newGazePos[1] = -(gazePos[1] - screen[1] / 2)
        curSamp = newGazePos

    return curSamp


def distBetweenPoints(p1, p2):
    """
    Calculates the distance between two points in a grid
    Parameters
    ----------
    p1 : tuple
        A tuple containing the (x,y) coordinates of the first point
    p2 : tuple
        A tuple containing the (x,y) coordinates of the second point
    Returns
    -------
    dist : float
        The Euclidian distance between the two points, the function assumes
        that the y-scaling and x-scaling are the same
    Examples
    --------
    >>> dist = distBetweenPoints((0,0), (10,10))
    >>> dist
    14.142135623730951
    """
    dist = np.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)
    return dist


def get_flanker_name(target_name='bL', list_cond=['bL', 'pL', 'bR', 'pR'],
                     num_fl=4,
                     same_ori=True, same_color=True):
    """
    Get flanker name given target name
    """

    if same_ori and same_color:  # everything the same (should not be the case)
        flank_name = list(np.repeat(target_name, num_fl))

    elif same_ori and not same_color:  # same orientation
        flank_name = [val for val in list_cond if val.endswith(target_name[-1])]
        flank_name = list(np.repeat(flank_name, num_fl / len(flank_name)))

    elif not same_ori and same_color:  # same color
        flank_name = [val for val in list_cond if val.startswith(target_name[0])]
        flank_name = list(np.repeat(flank_name, num_fl / len(flank_name)))

    else:  # all different
        flank_name = list_cond

    np.random.shuffle(flank_name)  # randomize

    if num_fl < len(
            list_cond):  # if less flankers, subselect --> should change this, and add condition where more flankers available
        flank_name = flank_name[:num_fl]

    return flank_name


def get_response4staircase(event_key=[], target_key=[]):
    """ helper function to get responses for crowding task """

    if event_key in target_key:
        response = 1
        print('correct answer')
    else:
        response = 0
        print('wrong answer')

    return response


class StaircaseCostum():
    """
    Costum staircase class - X Up Y Down
    """

    def __init__(self,
                 startVal,
                 stepSize=.1,  # stepsize
                 nUp=1,
                 nDown=3,  # correct responses before stim goes down
                 minVal=0,
                 maxVal=1):

        """ Initializes object and set relevant variables """

        # input variables
        self.startVal = startVal

        self.nUp = nUp
        self.nDown = nDown

        self.stepSize = stepSize

        self.minVal = minVal
        self.maxVal = maxVal

        self.data = []
        self.intensities = [startVal]
        self.reversalIntensities = []

        # correct since last stim change (minus are incorrect):
        self.correctCounter = 0
        self.incorrectCounter = 0
        self._nextIntensity = startVal

        self.increase = False
        self.decrease = False

    def addResponse(self, result):

        """ add pp response to staircase """

        # add response to data
        self.data.append(result)

        # increment the counter of correct scores
        if result == 1:

            self.correctCounter += 1

            if self.correctCounter >= self.nDown:
                self.decrease = True
                # reset counter
                self.correctCounter = 0

        elif result == 0:

            self.incorrectCounter += 1

            if self.incorrectCounter >= self.nUp:
                self.increase = True
                # reset both counters
                self.correctCounter = 0
                self.incorrectCounter = 0

        # calculate next intensity
        self.calculateNextIntensity()

    def calculateNextIntensity(self):

        """ calculate next value to use """

        # add reversal info
        if self.increase or self.decrease:
            self.reversalIntensities.append(self.intensities[-1])

        if self.increase:

            self._nextIntensity += self.stepSize

            # check we haven't gone out of the legal range
            if (self.maxVal is not None) and (self._nextIntensity > self.maxVal):
                self._nextIntensity = self.maxVal

            self.increase = False

        elif self.decrease:

            self._nextIntensity -= self.stepSize

            # check we haven't gone out of the legal range
            if (self.minVal is not None) and (self._nextIntensity < self.minVal):
                self._nextIntensity = self.minVal

            self.decrease = False

        # append intensities
        self.intensities.append(self._nextIntensity)

    def mean(self):
        return np.array(self.intensities).mean()

    def sd(self):
        return np.array(self.intensities).std()


def get_flanker_pos(num_fl=4, offset_ang=45, distance_r=.8, hemi='right',
                    ecc=8):
    """ define distractor positions

    Parameters
    ----------
    num_fl : int
        number of flankers
    offset_ang: float
        angle in degrees to offset from 0
    distance_r: float
        ratio of ecc (to calculate radial distance between target and flank)
    hemi: str
        visual hemifield we're plotting in
    ecc: int/float
        eccentricity in dva
    
    Outputs
    -------
    fl_pos : list
        list of [x,y] positions per distractor
    
    """

    hypotenuse = distance_r * ecc
    fl_angles = [offset_ang + (360 / num_fl) * i for i in np.arange(num_fl)]

    fl_pos = []

    for num in range(num_fl):
        fl_pos.append(list([hypotenuse * np.cos(np.deg2rad(fl_angles[num])),
                            hypotenuse * np.sin(np.deg2rad(fl_angles[num]))]))

        # update x position given hemifield of stim
        fl_pos[num][0] += ecc if hemi == 'right' else -ecc

    return fl_pos


def update_dots(ElementArrayStim, elem_positions=[], grid_pos=[], contrast=.2, opac=.3):  # 4):

    """ quick fix func to update dot element array settings
    should refurbish
    
    Parameters
    ----------
    ElementArrayStim: Psychopy object
    	ElementArrayStim to be updated 
    elem_positions: arr
         numpy array with element positions to be updated and shown (N,2) -> (number of positions, [x,y])
         to be used for opacity update
    grid_pos: arr
        numpy array with element positions (N,2) of whole grid -> (number of positions, [x,y])
        
    """

    # set number of elements
    nElements = grid_pos.shape[0]

    element_positions = grid_pos.copy()

    # update dot positions
    element_positions[:elem_positions.shape[0]] = np.array([elem_positions])

    # set element contrasts
    element_contrast = np.zeros(nElements)
    element_contrast[:elem_positions.shape[0]] = contrast

    # set element opacities
    element_opacities = np.zeros(nElements)
    element_opacities[:elem_positions.shape[0]] = opac

    # set all of the above settings
    ElementArrayStim.setOpacities(element_opacities)
    ElementArrayStim.setContrs(element_contrast)
    ElementArrayStim.setXYs(element_positions)

    return (ElementArrayStim)
