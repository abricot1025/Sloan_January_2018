"""

MOONS OPS to FPS Shared Library prototype

A prototype module containing functions shared between the MOONS
Observation Preparation Software, Fibre Positioner Control Software
and Fibre Positioner Collision Avoidance Software

This prototype is written in Python. It is expected that the final version
will be translated into C or C++ for greater efficiency and compatibility
with ESO/VLT software.

08 Jul 2013: First version as fibre positioner simulator.
24 Jul 2014: Converted into fibre positioner control prototype.
26 Jan 2015: Second version as OPS to FPS shared library prototype.
12 May 2015: Beta arm travel limited to -180 to +150 degrees.
11 Jan 2016: Third version with new avoidance zones and new fibre positioner
             characteristics.
15 Jan 2016: Tested and debugged. Log messages improved. Graphical debugging
             option added. Positioner plots now show the ID.
18 Jan 2016: Added a collision zone representing the "foot" of the beta arm,
             which actuates the beta datum switch.
19 Jan 2016: Added path analysis ellipses and a safety padding around each
             avoidance zone. Added conflict type counters.
02 Feb 2016: Corrected allocation of near and far neighbours within a grid
             and corrected limiting reach for checking far neighbours.
             Changed "merit" to "degree of difficulty". Added the option to
             make the fibre snag conflict check symmetrical (i.e. regard both
             positioners as in conflict). check_pair function added.
03 Feb 2016: Positioners looked up by ID rather than row and column.
             check_pair function modified to use target tuples.
17 Feb 2016: check_targets function reworked to choose the best parity for
             each positioner. Low-level functions are now check_pair and
             test_pair.
24 Mar 2016: Fail gracefully if the plotting module cannot be imported.
10 May 2016: Improved logging. Make sure angles are reported consistently
             in degrees.
28 Jun 2016: Temporarily added Motor code from path analysis, until the
             path analysis code is successfully merged.
01 Jul 2016: Corrected a typo in Positioner.target_angles().
19 Jul 2016: Ensure that target_angles() keeps the alpha and beta angles
             within their defined limits. Calculate the alpha and beta
             movements as differences from their default positions.
             Some modifications to equalise the outputs from the Python
             and C++ versions of the code.
26 Jul 2016: Added set_arm_angles() function, so that the starting angles
             can be defined as parameters. Added weights parameter to
             difficulty_values function. Adjusted avoidance zone sizes and
             removed the circular zone at the end of the datum foot.
27 Jul 2016: target_angles() renamed to get_arm_angles(). Implemented the
             alpha zero-point angle, self.orient, for each positioner.
             Added arm motor tests.
27 Jan 2018: set_target() updated to correct the target position of the
             positioner to depends upon the fiber used. Furthermore, the
             set_target() have been updated to contain the fiber attribution
             necessary for the Sloan project.

             TODO: Update the correction added to the target positioner to depend on the parity of the target.
             
THIS CODE HAS BEEN CONVERTED TO C++.

@author: Steven Beard (UKATC)

"""
# Python 3 emulation under Python 2.7
from __future__ import absolute_import, unicode_literals, division, print_function

# Select a logging level.
import logging
import circles
import path_generator as pg
# logging.basicConfig(level=logging.ERROR)  # Errors only
# logging.basicConfig(level=logging.WARN)   # Warnings and errors only
logging.basicConfig(level=logging.INFO)   # Informational output 
# logging.basicConfig(level=logging.DEBUG)  # Some debugging output
logger = logging.getLogger("fps_classes")

# Use this flag to turn graphical debugging on and off.
# When this flag is True, plots of conflicting positioners will be
# shown in more detail. This can generate a lot of graphics.
GRAPHICAL_DEBUGGING = False

# The beta and datum avoidance zones overlap considerably. They clutter up
# the plots, and and it may be a waste of time including both of them.
# These flags can be used to see the result of turning the zones on and
# off. The flags can also be used to turn plots of the ellipses and padding
# on and off.
INCLUDE_BETA_ZONE = True
INCLUDE_DATUM_FOOT = True
PLOT_ELLIPSES = False
PLOT_PADDING = False

# Import standard Python libraries
import math, time

# Import common utility and plotting libraries.
# If any package is not available, try the current directory.
try:
    import mocpath.util as util
except ImportError:
    import util as util
try:
    import mocpath.plotting as plotting
except ImportError:
    try:
        import plotting as plotting
    except ImportError:
        plotting = None
        logger.warning("Could not import plotting utility. Plotting disabled.")

import numpy as np
import parameters as param

# Import constants from the utility module.
EPS = util.EPS            # Smallest possible floating point increment
MAXFLOAT = util.MAXFLOAT  # Largest possible floating point value.

ROOT3 = util.ROOT3
ROOT3BY2 = util.ROOT3BY2
PI2 = util.PI2
PIBY2 = util.PIBY2

PARITY_RIGHT = util.PARITY_RIGHT
PARITY_LEFT = util.PARITY_LEFT

# Global variables that define the default lengths and widths
# of the alpha and beta arms in millimetres.
#
# NOTE: In the final software, these parameters will be imported
# from a configuration file, not defined as constants.
#
# TODO: READ FROM CONFIGURATION FILE
INS_POS_LENGTH1 = param.INS_POS_LENGTH1   # Length of alpha arm (mm).
INS_POS_WIDTH1 = param.INS_POS_WIDTH1      # Width of alpha arm (mm).
ALPHA_LIMITS = [param.ALPHA_TRAVEL_MIN, param.ALPHA_TRAVEL_MAX]  # Min and max travel limits of alpha arm (deg)
ALPHA_DEFAULT = param.ALPHA_DEFAULT       # Default (starting) location for alpha arm (deg).
ALPHA_RPM_LIMITS = [param.ALPHA_RPM_MIN, param.ALPHA_RPM_MAX] # Min and max motor speeds for alpha arm

INS_POS_LENGTH2 = param.INS_POS_LENGTH2    # Length of beta arm (mm).
INS_POS_WIDTH2 = param.INS_POS_WIDTH2      # Width of beta arm (mm).
BETA_LIMITS  = [param.BETA_TRAVEL_MIN, param.BETA_TRAVEL_MAX]          # Min and max travel limits of beta arm (deg)
BETA_DEFAULT = param.BETA_DEFAULT     # Default (starting) location for beta arm (deg).
BETA_RPM_LIMITS = [param.BETA_RPM_MIN, param.BETA_RPM_MAX] # Min and max motor speeds for beta arm

MAX_WAVEFORM_STEPS = param.MAX_WAVEFORM_STEPS  # Maximum number of time steps in a motor waveform

# Avoidance zone parameters
INS_POS_B1  = param.INS_POS_B1  # Length of new metrology zone (mm)
INS_POS_B2  = param.INS_POS_B2  # Rear protruding length of metrology zone (mm)
INS_POS_B3  = param.INS_POS_B3  # Additional length of new triangular avoidance zone.
INS_POS_B4  = INS_POS_LENGTH2 - (INS_POS_B1 + INS_POS_B2 + INS_POS_B3)
INS_POS_TB2 = param.INS_POS_TB2  # Length of small triangular avoidance zone (mm)
INS_POS_TB3 = param.INS_POS_TB3  # Length of large triangular avoidance zone.
INS_POS_TW1 = INS_POS_WIDTH2  # Width of small triangular avoidance zone.
INS_POS_TW2 = param.INS_POS_TW2   # Width of large triangular avoidance zone.
INS_POS_DL  = param.INS_POS_DL    # Datum actuator length
INS_POS_DW  = param.INS_POS_DW   # Datum actuator width
INS_POS_SAFETY = param.INS_POS_SAFETY  # Safety tolerance added to all avoidance zones.

# Fibre holder parameters
INS_POS_MINDIST = param.INS_POS_MINDIST    # Closest fibre approach distance (mm).
INS_POS_TOLER   = param.INS_POS_TOLER   # Fibre positioning tolerance (microns)
FIBRE_RADIUS    = INS_POS_MINDIST/2.0  # Radius of fibre holder in mm

# Focal plane radius of curvature (mm)
FOCAL_PLANE_CURVATURE = param.FOCAL_PLANE_CURVATURE
# Focal plane diameter (mm)
FOCAL_PLANE_DIAMETER = param.FOCAL_PLANE_DIAMETER


# ***
# ORIGINAL (EXTREMELY INEFFICIENT) CONFLICT FUNCTION
# SUPERCEDED BY THE NEW PositionerGrid.test_pair function.
def check_for_conflict(rcenpos1, thcenpos1, orient1,
                       rtarget1, thtarget1, parity1,
                       rcenpos2, thcenpos2, orient2,
                       rtarget2, thtarget2, parity2):
    """
    
    A function which implements the conflict check with a non-object-oriented 
    interface. This version of the API would be used if the shared library is
    implemented in C.
    
    :Parameters:
        
    rcenpos1: float
        The radial distance of the centre of the first FPU,
        in polar coordinates, with respect to the centre of
        the focal plane (in mm).
    thcenpos1: float
        The theta angle of the centre of the first FPU,
        in polar coordinates, with respect to the centre of
        the focal plane (in radians)
    orient1: float.
        The rotation angle between the focal plane coordinate system
        and the local coordinate system of the first positioner
        (in radians). 0.0 means no local rotation.
    rtarget1: float
        The radial distance of the target assigned to the first FPU,
        in polar coordinates, with respect to the centre of
        the focal plane (in mm).
    thtarget1: float
        The theta angle of the target assigned to the first FPU,
        in polar coordinates,  with respect to the centre of
        the focal plane (in radians).
    parity1: int
        The arm elbow orientation to be adopted by the target
        assigned to the first FPU:
            
        * 1 means elbow right armed
        * -1 means elbow left armed
        
     rcenpos2: float
        The radial distance of the centre of the second FPU,
        in polar coordinates, with respect to the centre of
        the focal plane (in mm).
    thcenpos2: float
        The theta angle of the centre of the second FPU,
        in polar coordinates, with respect to the centre of
        the focal plane (in radians)
    orient2: float.
        The rotation angle between the focal plane coordinate system
        and the local coordinate system of the second positioner
        (in radians). 0.0 means no local rotation.
    rtarget2: float
        The radial distance of the target assigned to the second FPU,
        in polar coordinates, with respect to the centre of
        the focal plane (in mm).
    thtarget2: float
        The theta angle of the target assigned to the second FPU,
        in polar coordinates,  with respect to the centre of
        the focal plane (in radians).
    parity2: int
        The arm elbow orientation to be adopted by the target
        assigned to the second FPU:
            
        * 1 means elbow right armed
        * -1 means elbow left armed
   
    :Returns:

    in_conflict: boolean
        True if the two positioners are in conflict, or if any of the
        targets cannot be reached.
        False if both targets are ok and the positioners are not in conflict.
    
    """
    # NOTE: very inefficient because it creates, configures an destroys
    # the Positioner objects each time.
    
    # Create two positioner objects
    positioner1 = Positioner(1, r_centre_focal=rcenpos1,
                             theta_centre_focal=thcenpos1,
                             orient=orient1, column=0, row=0)
    positioner2 = Positioner(2, r_centre_focal=rcenpos2,
                             theta_centre_focal=thcenpos2,
                             orient=orient2, column=0, row=1)
    
    # Assign the targets
    result1 = positioner1.set_target( rtarget1, thtarget1, parity1 )
    result2 = positioner2.set_target( rtarget2, thtarget2, parity2 )
    
    if result1 and result2:
        # Check for conflict
        in_conflict = positioner1.in_conflict_with( positioner2 )
        return in_conflict
    else:
        # Targets cannot be reached.
        return True



# FROM PATH ANALYSIS RELEASE - 13 JUN 2016
# FIXME: NEEDS REWORK. USE ALPHA/BETA CONSTANTS. CAN BE REPLACED BY ONE CLASS.
class Motor(object):
    # STEVEN: THIS WORKS FOR BOTH MOTORS
    def __init__(self,sim_length, position_target, max_speed, min_speed):
        self.position_array = [0.0] * sim_length
        self.speed_array = [0.0] * (sim_length - 1)
        self.position_array[sim_length - 1] = position_target
        self.maximum_speed = max_speed
        self.minimum_speed = min_speed


class Positioner(object):
    """
    
    Class representing a MOONS fibre positioner unit (FPU).
    
    Each positioner is identified by a unique identifier and a location
    in MOONS focal plane coordinates. Coordinates within a MOONS fibre
    positioner can be described in 4 different ways:
    
    1) Focal plane polar coordinates: r,theta relative to the centre
       of the focal plane.
       
    2) Focal plane Cartesian coordinates: x,y relative to the centre
       of the focal plane.
       
    3) Local polar coordinates: r,theta relative to the centre of
       the positioner.
       
    4) Local Cartesian coordinates: x,y relative to the centre of
       the positioner.

    External interfaces use focal plane coordinates (1) or (2) because
    those are universal over the whole focal plane, and can detect
    conflicts between positioners.
    
    Internal interfaces use local coordinates (3) or (4), which are
    valid within the patrol zone of one positioner only.
    
    """
    # Conflict types
    CONFLICT_OK = 0
    CONFLICT_UNREACHABLE = 1
    CONFLICT_LOCKED = 2
    CONFLICT_TOOCLOSE = 3
    CONFLICT_METROLOGY = 4
    CONFLICT_BETA = 5
    CONFLICT_DATUM = 6
    CONFLICT_FIBRE_SNAG = 7
    
    # To start with, common parameter are not defined.
    _params_defined = False
    
    # ***
    @classmethod
    def define_common_params(cls, length1, length2):
        """
        
        Class method which defines the parameters which are common to all
        fibre positioners. These parameters only need to be calculated once.
        
        :Parameters:
        
        """
        # Define the parameters valid for all positioner arms (saving the
        # square of those parameters for optimisation purposes).
        cls.length1 = float(length1)
        cls.length2 = float(length2)
        cls._length1sq = cls.length1 * cls.length1
        cls._length2sq = cls.length2 * cls.length2
        # Inner and outer limits of the positioner patrol zone.
        cls.outer = cls.length1 + cls.length2
        cls.inner = abs(cls.length2 - cls.length1)
        cls._outersq = cls.outer * cls.outer
        cls._innersq = cls.inner * cls.inner
        # A safe distance beyond which positioners cannot conflict
        cls.safedist = cls.outer + INS_POS_MINDIST
        cls._safedistsq = cls.safedist * cls.safedist
        logger.debug("Target distance beyond which positioners cannot " + \
                     "conflict is %f" % cls.safedist)
        # A limiting distance inside which far neighbours cannot conflict
        cls.limitdist = (cls.outer - INS_POS_MINDIST) * (ROOT3 - 1.0)
#         cls._limitdistsq = cls.limitdist * cls.limitdist
        logger.debug("Reach distance inside which far neighbours cannot " + \
                     "conflict is %f" % cls.limitdist)

        # The first (length2-length1)/2.0 of the beta arm is inside the safe
        # zone. Since the boundary of the safe zone is a circular arc, the
        # following calculation takes into account the finite width of the
        # beta arm and calculates the length at the edges of the arm.
        cls._bsafesq = (((cls.length2-cls.length1)/2.0) * \
                            ((cls.length2-cls.length1)/2.0)) - \
                        (INS_POS_WIDTH2/2.0 * INS_POS_WIDTH2/2.0)
        cls.bsafe = math.sqrt(cls._bsafesq)
        logger.debug("Portion of beta arm safe from conflict is %f" % cls.bsafe)
        
        # Determine the critical extension length beyond which a
        # limited range of travel for the beta arm makes a difference
        # to the choice of parity.
        if BETA_LIMITS and BETA_LIMITS[1] < 180.0:
            upperlimit = math.radians(BETA_LIMITS[1])
            lcritsq = (cls.length1 * cls.length1) + \
                      (cls.length2 * cls.length2) - \
                      (2.0 * cls.length1 * cls.length2 * math.cos(upperlimit))
            cls.lengthcrit = math.sqrt(lcritsq)
            cls.paritycrit = PARITY_RIGHT
        elif BETA_LIMITS and BETA_LIMITS[0] > -180.0:
            upperlimit = math.radians(BETA_LIMITS[0])
            lcritsq = (cls.length1 * cls.length1) + \
                      (cls.length2 * cls.length2) - \
                      (2.0 * cls.length1 * cls.length2 * math.cos(upperlimit))
            cls.lengthcrit = math.sqrt(lcritsq)
            cls.paritycrit = PARITY_LEFT
        else:
            cls.lengthcrit = cls.length2
            cls.paritycrit = PARITY_RIGHT
        logger.debug("Only %s parity is possible beyond %f" % \
                    (util.elbow_parity_str(cls.paritycrit), cls.lengthcrit))
        # The common parameters are now defined.
        cls._params_defined = True
    
    def __init__(self, ident, r_centre_focal, theta_centre_focal, orient,
                 column=None, row=None, simulated=True, locked=False,
                 rfibre=None, thfibre=None, pfibre=PARITY_RIGHT,
                 length1=INS_POS_LENGTH1, length2=INS_POS_LENGTH2):
        """
        
        Constructor. The location of the positioner can be provided in
        either polar coordinates or Cartesian coordinates. If both are
        provided, the Cartesian coordinates are ignored.
        
        :Parameters:
        
        ident: int
            Unique identifier for this fibre positioner
        r_centre_focal: float
            The radial distance of the FPU centre, in polar coordinates,
            with respect to the centre of the focal plane (in mm)
        theta_centre_focal: float
            The theta angle of the FPU centre, in polar coordinates,
            with respect to the centre of the focal plane (in radians)
        orient: float
            The rotation angle between the focal plane coordinate system
            and the FPU local coordinate system (in radians).
            Defaults to 0.0 (no local rotation)
        column: int (optional)
            Column number of this positioner (if part of a grid)
        row: int (optional)
            Row number of this positioner (if part of a grid)
        simulated: bool (optional)
            True if the positioner is simulated. The default is True
            (since this code is a simulation).
        locked: bool (optional)
            True if the positioner is locked (broken). The default is False.
        rfibre: float (optional)
            Starting location for the fibre.
            Defaults to the initialised position, near datum.
        thfibre: float (optional)
            Starting location for the fibre
            Defaults to the initialised position, near datum.
        pfibre: int (optional)
            Starting parity for the fibre.
            Defaults to right-armed parity.
        length1, float (optional)
            Length of alpha arm in mm.
            Defaults to the designed length, INS_POS_LENGTH1
        length2, float (optional)
            Length of beta arm in mm.
            Defaults to the designed length, INS_POS_LENGTH2
   
        :Returns:

        New Positioner object
                    
        """
        # Define a unique name for this positioner.
        self.ident = int(ident)
        self.name = "POS%d" % self.ident
        logger.debug("Creating positioner %s" % self.name)
        
        # Define the location of the centre and the orientation of the
        # positioner on the MOONS focal plane, in polar and Cartesian
        # coordinates.
        self.r_centre_focal = float(r_centre_focal)
        self.theta_centre_focal = float(theta_centre_focal)
        (self.x_centre_focal, self.y_centre_focal) = \
            util.polar_to_cartesian(self.r_centre_focal, self.theta_centre_focal)
        self.orient = float(orient)
            
        # If the positioner is part of a grid, record its  location
        # in hexagonal grid coordinates (if defined).
        if column is not None and row is not None:
            self.column = int(column)
            self.row = int(row)
        else:
            self.column = None
            self.row = None
        
        # Is the positioner simulated or broken?
        self.simulated = bool(simulated)
        self.locked = bool(locked)
     
        # If not already defined, define the parameters common to all
        # fibre positioners.
        if not Positioner._params_defined:
            Positioner.define_common_params(length1, length2)

        # Start with no target assigned.
        # Start at the initialised location, near datum, unless a specific
        # starting location has been defined.
        self.starting_rfibre = rfibre
        self.starting_thfibre = thfibre
        self.starting_parity = pfibre
        self.initialise()
        
        # Initialise the list of pointers to neighbouring positioners
        self.near_neighbours = []
        self.far_neighbours = []

        # priority:
        #  by default the priority of each positioner is defined at 1 (the lowest priority),
        #  the highest priority is 4 and the lowest is 1
        # [0]: is the priority that we will used for the algorithm (the one from OPS file), ITS REMAINS CONSTANT !!!
        # [1]: is the pseudo_algorithm used in case two positioners have the same priority !
        self.OPS_priority=1
        self.pseudo_priority=1

        # STATUS (+ some explanation...)
        # each positioner has a status at its creation
        # - the first indicates the positioner status(ON=1,OFF=0):
        #   ON: has its attractive force to target activated, OFF: have its attractive force deactivated
        # - the second term is the ident of the neighbour positioner that put this agent positioner statusForceAttr=OFF
        #   0: if NO neigh positioners changed the agent positioners statusForceAttr to OFF
        # - The third element is the copy of the priority of the neighbour positioner that put
        #   this agent positioner statusForceAttr to OFF
        # - the fourth element:
        #  when the positioner arrives at destination, we remove its capacity to put the other positioners status[0,1,2]
        #  basically tell if a positioners arrives to its target
        #  1 ==INFLUENCABLE ==>  can influence the other positioners's status[0,1,2]
        #  0 ==NOT_INFLUENCABLE ==> cannot influence the other positioners's status[0,1,2]
        # - The fifth element:
        # when there is a deadlock cause by positioners with same priority[0], then statusDeadlock_sp becomes 1,
        # meaning that it goes a layers of priority above where additional action are used to solve the deadlock.
        # Values of statusDeadlock_sp, which is the same as its neigh.statusDeadlock_sp positioners, is either the agent ident or the
        # neigh ident value

        self.statusForceAttr = 1
        self.statusMemID = 0
        self.statusMemPr = 0
        self.statusArr = 1
        self.statusDeadlock_sp=0 # deadlock_samepriority




        #Average velocity of the motor (motor 1 and motor 2)
        # used in pd_dnf to determine when
        self.MeanVelMotor=[0,0]

        self.repulsiveTarget = np.array((0,0))

    def __str__(self):
        """
        
        Return a readable string describing the positioner.
   
        :Returns:

        strg: str
            Readable string
        
        """
        strg = "Positioner object \'%s\' " % self.name
        if self.column is not None and self.row is not None:
            strg += "[%d,%d] " % (self.column, self.row)
        theta_degrees = math.degrees(self.theta_centre_focal)
        strg += "located on the focal plane at (r=%.3f, theta=%.3f deg) " % \
            (self.r_centre_focal, theta_degrees)
        strg += "or (x=%.3f, y=%.3f).\n" % \
            (self.x_centre_focal, self.y_centre_focal)
        if abs(self.orient) > 0.0:
            orient_degrees = math.degrees(self.orient)
            strg += "  Positioner oriented by a0=%.3f (deg)." % orient_degrees
        strg += "  There are %d near and %d far neighbouring positioners.\n" % \
            (len(self.near_neighbours), len(self.far_neighbours))
            
        if self.locked:
            strg += "  Positioner is LOCKED. No target can be assigned."
        else:
            if self.target_assigned:
                strg += "  Target assigned. "
            else:
                strg += "  Target NOT assigned. "
        theta_degrees = math.degrees(self.theta_fibre_focal)
        strg += "Fibre located on the focal plane at (r=%.3f, theta=%.3f deg) " % \
            (self.r_fibre_focal, theta_degrees)
        strg += "or (x=%.3f, y=%.3f),\n" % \
            (self.x_fibre_focal, self.y_fibre_focal)
        theta_local = math.degrees(self.theta_fibre_local)
        strg += "    equivalent to local coordinates (r=%.3f, theta=%.3f deg)" % \
            (self.r_fibre_local, theta_local)
        strg += " or (x=%.3f, y=%.3f).\n" % \
            (self.x_fibre_local, self.y_fibre_local)

        strg += "  Elbow %s parity leads to elbow location" % \
            util.elbow_parity_str(self.target_parity)
        strg += " of (x=%.3f, y=%.3f) in local coordinates." % \
            (self.x_elbow_local, self.y_elbow_local)
            
        if self.in_conflict:
            strg += "\n  Positioner IN CONFLICT: " + self.conflict_reason
        return strg

# FROM PATH ANALYSIS RELEASE - 13 JUN 2016
#     def is_neighbour(self, other):
#         """
#         
#         Returns true if the two positioners are neighbours
# 
#         :Parameters:
# 
#         other: Positioner object
#             other positioner (ignored if None)
# 
#         """
#         if other is not None:
#             assert isinstance(other, Positioner)
#             if math.sqrt((self.x_centre_focal - other.x_centre_focal)**2 + \
#                (self.y_centre_focal - other.y_centre_focal)**2) < 2*PITCH + 1:
#                 return True
#             else:
#                 return False

    def add_neighbour(self, neighbour, near=True):
        """
        
        Add a link to a neighbouring positioner whose patrol zone
        overlaps with this one.
        
        :Parameters:

        neighbour: Positioner object
            Neighbouring positioner (ignored if None)
        near: bool (optional)
            Set True if the other positioner is a near neighbour
            or False if the other positioner is a far neighbour.
            Defaults to True.
        
        """
#         print("%s: Adding neighbour %s" % (self.name, neighbour.name))
        if neighbour is not None:
            assert isinstance(neighbour, Positioner)
            if near:
                self.near_neighbours.append(neighbour)
            else:
                self.far_neighbours.append(neighbour)

    def get_neighbours(self):
        """
        
        Return a list of neighbours of the current positioner.
        
        """
        # Append the lists of near and far neighbours.
        return self.near_neighbours + self.far_neighbours

    def initialise(self):
        """
        
        Move the positioner to its initial location.
        
        :Parameters:

        None
        
        """
        # Start with no conflict and a blank conflict message string
        self.in_conflict = False
        self.conflict_type = self.CONFLICT_OK
        self.conflict_reason = ''

        # At first there is no target assigned
        self.target_assigned = False
        # The positioner starts at its default starting position,
        # unless a different starting location has been specified.
        if self.starting_rfibre is not None and self.starting_thfibre is not None:
            # An alternative start location has been defined.
            result = self.set_target(self.starting_rfibre, self.starting_thfibre,
                            self.starting_parity, starting_location=True)
            if not result:
                strg = "Cannot move positioner %d to its initial location (%f,%f)" \
                    % (self.ident, self.starting_rfibre, self.starting_thfibre)
                raise ValueError(strg)
        else:
            # Move to default location
            self.set_arm_angles( ALPHA_DEFAULT, BETA_DEFAULT )
            # Move to the starting location, doubled back to minimum extension.
#             self.x_fibre_local = -self.inner
#             self.y_fibre_local = 0.0
#             self.target_parity = PARITY_RIGHT
#             (self.r_fibre_local, self.theta_fibre_local) = \
#                 util.cartesian_to_polar(self.x_fibre_local, self.y_fibre_local)
            strg = "Target initialised at local x=%f,y=%f --> r=%f,theta=%f at %s parity" % \
                (self.x_fibre_local, self.y_fibre_local,
                 self.r_fibre_local, math.degrees(self.theta_fibre_local),
                 util.elbow_parity_str(self.target_parity))
            if abs(self.orient) > 0.0:
                strg += " (orient=%f (deg))" % math.degrees(self.orient)
            logger.debug(strg)
        
#             self.x_fibre_focal = self.x_centre_focal + self.x_fibre_local
#             self.y_fibre_focal = self.y_centre_focal + self.y_fibre_local
#             (self.r_fibre_focal, self.theta_fibre_focal) = \
#                 util.cartesian_to_polar(self.x_fibre_focal, self.y_fibre_focal)
#             
#             # Solve the elbow location for this target
#             (self.x_elbow_local, self.y_elbow_local) = \
#                 util.solve_elbow_xy(self.x_fibre_local, self.y_fibre_local,
#                             self.target_parity, self.length1, self._length1sq,
#                             self.length2, self._length2sq)

# ***
    def can_move_to(self, r_target_focal, theta_target_focal, parity):
        """
        
        Determine whether this positioner is capable of moving to
        the specified target.
        
        NOTE: No collisions are detected by this function. It merely tests
        whether a target is within the patrol zone of this positioner.
        
        :Parameters:
        
        r_target_focal: float
            The radial distance of the target, in polar coordinates,
            with respect to the centre of the focal plane (in mm)
        theta_target_focal: float
            The theta angle of the target, in polar coordinates,
            clockwise from Y axis,
            with respect to the centre of the focal plane (in radians)
        parity: int
            The arm elbow orientation to be adopted:
            
            * 1 means elbow right armed
            * -1 means elbow left armed
   
        :Returns:

        target_reachable: boolean
            True if the target can be reached
            False if the target cannot be reached
        
        """
        logger.debug("Can positioner %s move to (r=%f, theta=%f) with %s parity?" %
                      (self.name, r_target_focal,
                       math.degrees(theta_target_focal),
                       util.elbow_parity_str(parity)))
        # Solve the triangle cosine rule to determine the local R
        # for the target
        rls = r_target_focal * r_target_focal + \
              self.r_centre_focal * self.r_centre_focal - \
              2.0 * r_target_focal * self.r_centre_focal * \
              math.cos(theta_target_focal-self.theta_centre_focal)
        r_target_local = math.sqrt(rls)
        
        # Can the positioner reach this target?
        # TODO: Should INS_POS_TOLER be taken into account?
        if r_target_local >= self.inner and r_target_local <= self.outer:
            # If the range of travel of the beta arm is limited, targets beyond
            # a certain critical distance can only be reached at one parity.
            if r_target_local > self.lengthcrit and parity != self.paritycrit:
                self.conflict_type = self.CONFLICT_UNREACHABLE
                self.conflict_reason = "%s: Wrong target parity at this reach " % \
                    self.name
                self.conflict_reason += "%.3f > %.3f" % (r_target_local, self.lengthcrit)
                logger.info(self.conflict_reason)
                return False
            else:
                logger.debug("YES.")
                return True
        else:
            self.conflict_type = self.CONFLICT_UNREACHABLE
            self.conflict_reason = "%s: Target outside patrol zone." % self.name
            #logger.info(self.conflict_reason)
            return False

    # ***
    def set_target(self, r_target_focal, theta_target_focal, parity,
                   priority, fiber, positioner, starting_location=False):
        """
        
        Propose a target for this positioner. If the target is reachable,
        the simulated positioner is moved immediately to the target.
        
        :Parameters:
        
        r_target_focal: float
            The radial distance of the target, in polar coordinates,
            with respect to the centre of the focal plane (in mm)
        theta_target_focal: float
            The theta angle of the target, in polar coordinates,
            clockwise from Y axis,
            with respect to the centre of the focal plane (in radians)
        parity: int
            The arm elbow orientation to be adopted:
            
            * 1 means elbow right armed
            * -1 means elbow left armed
        priority: int
            The priority of the target to be reached by the positioner
            4 is the highest priority and 1 is the lowest
   
        starting_location: bool (optional)
            If True, defines starting location, not a new target.
            The default is False.
   
        :Returns:

        target_reachable: boolean
            True if the target can be reached (target location stored)
            False if the target cannot be reached (target ignored)
        """

        # Laleh: It would be easier for PAS to handle targets already presented in local frame of each FPU
        # Can the positioner reach this target?
        if self.can_move_to(r_target_focal, theta_target_focal, parity):
            # Yes, target can be reached. If this is not the starting position,
            # and the positioner is not locked, the previous target is
            # overwritten.
            if not starting_location:
                if not self.locked:
                    self.target_assigned = True
                else:
                    # The positioner cannot be moved from its starting location.
                    self.in_conflict = True
                    self.conflict_type = self.CONFLICT_LOCKED
                    self.conflict_reason = "%s: Positioner locked." % self.name
                    return False

            # Update the fibre location to match the target, in focal plane
            # polar and Cartesian coordinates.

            fiber = int(fiber)
            self.fiber=fiber

            # The following lines allow the program to continue only if fiber conditions are filled
            # By default the program continues to work as it was working for MOONS
            # If conditions are not filled, the function returns False and the target is not set to this positioner

            positioner_IR = False
            positioner_fiducial = False
            positioner_empty = False

            for item in circles.IR:                         # Check in circles.py if the positioner is IR
                if positioner == item:
                    positioner_IR = True

            for item in circles.fiducial:                   # Check in circles.py if the positioner is a fiducial
                if positioner == item:
                    positioner_fiducial = True

            for item in circles.empty:                      # Check in circles.py if the positioner is empty
                if positioner == item:
                    positioner_empty = True

            if fiber == 1 and positioner_IR == False:       # If target must me surveyed with an IR fiber
                return False                                # but the positioner does not have one -> Stop targeting

            if positioner_fiducial == True:                 # If the positioner is a fiducial -> Stop targeting
                return False

            if positioner_empty == True:                    # If the positioner is a empty -> Stop targeting
                return False

            self.r_fibre_focal = r_target_focal
            self.theta_fibre_focal = theta_target_focal

            (self.x_fibre_focal, self.y_fibre_focal) = \
            util.polar_to_cartesian(self.r_fibre_focal,
                                        self.theta_fibre_focal)
            strg = "Target set to focal coordinates " + \
            "(r=%f,theta=%f) or (x=%f,y=%f)" % \
            (self.r_fibre_focal, math.degrees(self.theta_fibre_focal),
                 self.x_fibre_focal, self.y_fibre_focal)
            logger.debug(strg)
            
            # Determine the local (wrt positioner centre) coordinates for the
            # target.
            self.x_fibre_local = self.x_fibre_focal - self.x_centre_focal
            self.y_fibre_local = self.y_fibre_focal - self.y_centre_focal
            (self.r_fibre_local, self.theta_fibre_local) = \
                util.cartesian_to_polar(self.x_fibre_local, self.y_fibre_local)

            # Correct the position depending on the fiber used
            # Has to be updated to have the correction depending on the parity
            if fiber == 1:
                (self.r_fibre_local, self.theta_fibre_local) =\
                    (self.r_fibre_local + param.R_CORRECTION_IR, self.theta_fibre_local + param.THETA_CORRECTION_IR)
                self.fiber = 1
            elif fiber == 0:
                (self.r_fibre_local, self.theta_fibre_local) =\
                    (self.r_fibre_local + param.R_CORRECTION_VISIBLE, self.theta_fibre_local + param.THETA_CORRECTION_VISIBLE)
                self.fiber = 0
            elif fiber == -1:
                (self.r_fibre_local, self.theta_fibre_local) = \
                    (self.r_fibre_local + param.R_CORRECTION_CALIBRATION, self.theta_fibre_local + param.THETA_CORRECTION_CALIBRATION)
                self.fiber = -1

            self.target_parity = parity
            strg = "Target is at local coordinates " + \
                    "(r=%f,theta=%f) or (x=%f,y=%f) at %s parity" % \
                    (self.r_fibre_local, math.degrees(self.theta_fibre_local),
                    self.x_fibre_local, self.y_fibre_local,
                    util.elbow_parity_str(self.target_parity))
            logger.debug(strg)
         
            # Solve the elbow location for this target and save the result.
            (self.x_elbow_local, self.y_elbow_local) = \
                util.solve_elbow_xy(self.x_fibre_local, self.y_fibre_local,
                            self.target_parity, self.length1, self._length1sq,
                            self.length2, self._length2sq)
            self.in_conflict = False

            # assign the priority
            self.OPS_priority=priority
            self.pseudo_priority=priority

            pg.positioners.append(positioner)

            return True

        else:
            # No, target cannot be reached. The previous target is not overwritten.
            return False

    # ***
    def set_arm_angles(self, angle1, angle2):
        """
        
        Set the alpha and beta motor angles for this positioner.
        The simulated positioner is moved immediately to these angles.
        
        :Parameters:
        
        angle1: float
            Alpha motor angle, with respect to its datum location (in degrees)
        angle2: float
            Beta motor angle, with respect to its datum location (in degrees)
   
        :Returns:

        None
                
        """
        # Make sure alpha and beta angles are within the defined range
        # Ensure that the requested angles are within their allowed limits.
        while angle1 < ALPHA_LIMITS[0] :
            angle1 += 360.0
        while angle1 > ALPHA_LIMITS[1]:
            angle1 -= 360.0
        while angle2 < BETA_LIMITS[0] :
            angle2 += 360.0
        while angle2 > BETA_LIMITS[1]:
            angle2 -= 360.0
        angle1_rad = math.radians(angle1)
        angle2_rad = math.radians(angle2)
        
        if angle2 > 0.0:
            self.target_parity = PARITY_LEFT
        else:
            self.target_parity = PARITY_RIGHT
        
        alpha_angle = angle1_rad + self.orient
        self.x_elbow_local = self.length1 * math.cos(alpha_angle)
        self.y_elbow_local = self.length1 * math.sin(alpha_angle)
        
        elbow_angle = math.pi + alpha_angle + angle2_rad
        self.x_fibre_local = self.x_elbow_local + self.length2 * math.cos(elbow_angle)
        self.y_fibre_local = self.y_elbow_local + self.length2 * math.sin(elbow_angle)

        (self.r_fibre_local, self.theta_fibre_local) = \
            util.cartesian_to_polar(self.x_fibre_local, self.y_fibre_local)
        strg = "Arms moved to (alpha=%f,beta=%f) at orient=%f\n" % \
            (angle1,angle2, math.degrees(self.orient))
        strg += "  which moves fibre to local coordinates " + \
            "(r=%f,theta=%f) or (x=%f,y=%f) at %s parity" % \
            (self.r_fibre_local, math.degrees(self.theta_fibre_local),
                self.x_fibre_local, self.y_fibre_local,
                util.elbow_parity_str(self.target_parity))
        logger.debug(strg)
        
        self.x_fibre_focal = self.x_centre_focal + self.x_fibre_local
        self.y_fibre_focal = self.y_centre_focal + self.y_fibre_local
        (self.r_fibre_focal, self.theta_fibre_focal) = \
            util.cartesian_to_polar(self.x_fibre_focal, self.y_fibre_focal)
        strg = "Fibre is at focal coordinates " + \
                "(r=%f,theta=%f) or (x=%f,y=%f)" % \
                (self.r_fibre_focal, math.degrees(self.theta_fibre_focal),
                 self.x_fibre_focal, self.y_fibre_focal)
        logger.debug(strg)

    def get_arm_angles(self):
        """
        
        Calculate the alpha and beta motor angles that would be adopted for
        the current positioner target and parity.
        
        :Parameters:
        
        None
                
        :Returns:
        
        (angle1,angle2)
            Alpha and beta motor angles in degrees.
        
        """
        reachsq = self.r_fibre_local * self.r_fibre_local
        # Determine the third length of the triangle that needs
        # to be made by the arm's shoulder and elbow to reach
        # the fibre position.
        # The EPS in the following tests accounts for floating point rounding.
        if self.r_fibre_local > (self.outer - EPS):
            # Special case of arm at full stretch.
            angle1 = math.radians(90.0) - self.theta_fibre_local
            angle2 = -math.radians(180.0)
        elif self.r_fibre_local < (self.inner + EPS):
            # Special case of arm completely doubled back.
            angle1 =  (math.radians(90.0) - self.theta_fibre_local) + \
                        math.radians(180.0)
            angle2 = 0.0
        else:
            # Solve the triangle cosine rule to determine the angle between
            # the shoulder and the fibre target.
            shoulder_fibre = util.solve_triangle(self.r_fibre_local, reachsq,
                                self.length1, self._length1sq, self._length2sq )
            if shoulder_fibre is None:
                # Equation cannot be solved
                logger.error("+++ Shoulder to fibre angle equation cannot be solved")
                return (None, None)
#             print("shoulder_fibre=", math.degrees(shoulder_fibre))

            # Convert shoulder to fibre angle into alpha angle
            angle1 = util.solve_shoulder(self.theta_fibre_local, shoulder_fibre,
                                    self.target_parity)        
#             print("angle1=", math.degrees(angle1))

            # Solve the cosine rule again to determine the angle between
            # the shoulder and elbow.
            shoulder_elbow = util.solve_triangle(self.length1, self._length1sq,
                                            self.length2, self._length2sq,
                                            reachsq)
            if shoulder_elbow is None:
                # Equation cannot be solved
                logger.error("+++ Shoulder to elbow angle equation cannot be solved")
                return (None, None)
        
            # The beta angle is the negative of this angle. The order of
            # calculation depends on the parity adopted.
            if self.target_parity == PARITY_RIGHT:
                angle2 = -shoulder_elbow
            else:
                angle2 = shoulder_elbow
                
        # Subtract the alpha zero-point angle
        angle1 -= self.orient
                
        # Ensure that the calculated angles are within their allowed limits.
        angle1_deg = math.degrees(angle1)
        angle2_deg = math.degrees(angle2)
        while angle1_deg < ALPHA_LIMITS[0] :
            angle1_deg += 360.0
        while angle1_deg > ALPHA_LIMITS[1]:
            angle1_deg -= 360.0
        while angle2_deg < BETA_LIMITS[0] :
            angle2_deg += 360.0
        while angle2_deg > BETA_LIMITS[1]:
            angle2_deg -= 360.0

        return (angle1_deg, angle2_deg)              

    def difficulty_value(self, others=[],
                         weights=[0.4, 0.0, 0.1, 0.1, 0.2, 0.2]):
        """
        
        Calculate the difficulty value associated with the current positioner
        target and parity, given the targets assigned to the positioner's
        neighbours.
        
        :Parameters:
        
        others: list of Positioner (optional)
            A list of other positioners to be checked. If not specified,
            the positioner's neighbours will be checked.
        weights: list of 6 floats (optional)
            The weights to be applied to the difficulty calculation.
            
            * Alpha movement weight
            * Beta movement weight
            * Positioner extension/reach weight
            * Closest target-target approach weight
            * Closest metrology-metrology approach weight
            * Closest beta arm-beta arm approach weight
            
            The weights must add up to 1.
                
        :Returns:
        
        difficulty: float
            A number between 0.0 and 1.0, estimating the difficulty in
            achieving the given configuration with real hardware.
        
        """
        # TODO: INVESTIGATE THE BEST DIFFICULTY DEFINITION.
       
        # (1) Calculate the amount of alpha and beta axis movement required
        # to reach the target.
        (angle1,angle2) = self.get_arm_angles()
        difficulty_alpha = abs(angle1 - ALPHA_DEFAULT) / 360.0
        difficulty_beta = abs(angle2 - BETA_DEFAULT) / 360.0
        
        # (2) Calculate the degree to which this positioner's target intrudes
        # into the patrol field of a neighbouring positioner.
        # NOTE: Overlaps somewhat with the beta difficulty.
        reachsq = self.r_fibre_local * self.r_fibre_local
        middlesq = self._outersq / 4.0
        difficulty_reach = max(0.0, (reachsq - middlesq) / self._outersq)

        # (3) Calculate the closest distance between any two fibre
        # holders, belonging to this positioner and its neighbours, taking
        # into account the minimum target approach distance, INS_POS_MINDIST.
        closestsq = self._outersq
        if others:
            check_list = others
        else:
            check_list = self.near_neighbours + self.far_neighbours
        for neighbour in check_list:
            distsq = util.distance_squared(self.x_fibre_focal,
                                           self.y_fibre_focal,
                                           neighbour.x_fibre_focal,
                                           neighbour.y_fibre_focal)
            if distsq < closestsq:
                closestsq = distsq
        # Base the difficulty on the square of the distance ratio, so the
        # value rises more steeply at close distances.
        ratio = max(0.0, (math.sqrt(closestsq) - INS_POS_MINDIST)) / self.outer
        difficulty_target_target = max(0.0, 1.0 - math.sqrt(ratio))
        
        # (4) Calculate the closest approach distance between any two
        # metrology avoidance ellipses or beta avoidance ellipses
        closestsq_met = self._outersq
        closestsq_beta = self._outersq
        if others:
            check_list = others
        else:
            check_list = self.near_neighbours + self.far_neighbours
        my_met_ellipse = self.avoidance_ellipse_metrology()
        my_met_x = my_met_ellipse[0]
        my_met_y = my_met_ellipse[1]
        my_beta_ellipse = self.avoidance_ellipse_beta()
        my_beta_x = my_beta_ellipse[0]
        my_beta_y = my_beta_ellipse[1]
        for neighbour in check_list:
            other_met_ellipse = neighbour.avoidance_ellipse_metrology()
            other_met_x = other_met_ellipse[0]
            other_met_y = other_met_ellipse[1]
            other_beta_ellipse = neighbour.avoidance_ellipse_beta()
            other_beta_x = other_beta_ellipse[0]
            other_beta_y = other_beta_ellipse[1]
            distsq_met = util.distance_squared(my_met_x, my_met_y,
                                               other_met_x, other_met_y)
            distsq_beta = util.distance_squared(my_beta_x, my_beta_y,
                                               other_beta_x, other_beta_y)
            if distsq_met < closestsq_met:
                closestsq_met = distsq_met
            if distsq_beta < closestsq_beta:
                closestsq_beta = distsq_beta
        # Base the difficulties on the square root of the distance ratio, so the
        # value rises more steeply at close distances.
        ratio = max(0.0, (math.sqrt(closestsq_met) - INS_POS_WIDTH2)) / self.outer
        difficulty_ellipse_met = max(0.0, 1.0 - math.sqrt(ratio))
        ratio = max(0.0, (math.sqrt(closestsq_beta) - INS_POS_MINDIST)) / self.outer
        difficulty_ellipse_beta = max(0.0, 1.0 - math.sqrt(ratio))

        # Try varying these weights. The total must add up to 1.0
        if not isinstance(weights, (tuple,list)) or len(weights) < 6:
            weights=[0.4, 0.0, 0.1, 0.1, 0.2, 0.2]
        difficulty_alpha_weight = weights[0]
        difficulty_beta_weight = weights[1]
        difficulty_reach_weight = weights[2]
        difficulty_target_target_weight = weights[3]
        difficulty_ellipse_met_weight = weights[4]
        difficulty_ellipse_beta_weight = weights[5]

        # Return a combined difficulty estimate.
        strg = "Target/parity difficulties for %d: " % self.ident
        strg += "alpha=%f, beta=%f, reach=%f, target=%f, ellipse_met=%f, ellipse_bet=%f" % \
            (difficulty_alpha, difficulty_beta, difficulty_reach,
             difficulty_target_target, difficulty_ellipse_met,
             difficulty_ellipse_beta)
        logger.debug( strg )
        difficulty = difficulty_alpha_weight * difficulty_alpha + \
                difficulty_beta_weight * difficulty_beta + \
                difficulty_reach_weight * difficulty_reach + \
                difficulty_target_target_weight * difficulty_target_target + \
                difficulty_ellipse_met_weight * difficulty_ellipse_met + \
                difficulty_ellipse_beta_weight * difficulty_ellipse_beta
        return difficulty
        
    def avoidance_polygon_metrology(self, padding=0.0):
        """
        
        Return, in focal plane coordinates, the vertices of the
        polygon bounding the main metrology target area,
        defined as ((x1, y2), (x2, y2), (x3, y3), (x4, y4), (x5, y5)).
        
        :Parameters:
        
        padding: float (optional)
            Safety padding added around the edge of the avoidance zone.
            Defaults to 0.0
                
        :Returns:
        
        ((x1, y2), (x2, y2), (x3, y3), (x4, y4), (x5, y5))
        
        """
        # TODO: These avoidance zone functions need to be optimised.
        # Calculate the vertices of a rectangle which encloses the fibre
        # holder and the metrology targets.
        xdiff = self.x_fibre_local - self.x_elbow_local
        ydiff = self.y_fibre_local - self.y_elbow_local
        
        half_width = padding + INS_POS_WIDTH2/2.0
        xdelta = half_width * ydiff / self.length2
        ydelta = half_width * xdiff / self.length2

        metrology_fraction = INS_POS_B1 / self.length2
        padding_fraction = padding / self.length2
        xmetrology_end = self.x_fibre_focal - metrology_fraction * xdiff
        ymetrology_end = self.y_fibre_focal - metrology_fraction * ydiff

        triangle_fraction = INS_POS_TB2 / self.length2
        xtriangle_end = xmetrology_end - \
            (triangle_fraction + padding_fraction) * xdiff
        ytriangle_end = ymetrology_end - \
            (triangle_fraction + padding_fraction) * ydiff
        xmetrology_start = self.x_fibre_focal + padding_fraction * xdiff
        ymetrology_start = self.y_fibre_focal + padding_fraction * ydiff
            
        x1 = xmetrology_end - xdelta
        y1 = ymetrology_end + ydelta
        x2 = xtriangle_end
        y2 = ytriangle_end
        x3 = xmetrology_end + xdelta
        y3 = ymetrology_end - ydelta
        x4 = xmetrology_start + xdelta
        y4 = ymetrology_start - ydelta
        x5 = xmetrology_start - xdelta
        y5 = ymetrology_start + ydelta
        return ((x1,y1),(x2,y2),(x3,y3),(x4,y4),(x5,y5))

    def avoidance_rectangle_metrology(self, padding=0.0):
        """
        
        Return, in focal plane coordinates, the vertices of the
        rectangle bounding the main metrology target area,
        defined as ((x1, y2), (x2, y2), (x3, y3), (x4, y4)).
        
        :Parameters:
        
        padding: float (optional)
            Safety padding added around the edge of the avoidance zone.
            Defaults to 0.0
                                
        :Returns:
        
        ((x1, y2), (x2, y2), (x3, y3), (x4, y4))
        
        """
        # TODO: These avoidance zone functions need to be optimised.
        # TODO: REPLACE WITH AVOIDANCE_POLYGON_METROLOGY
        # Calculate the vertices of a rectangle which encloses the fibre
        # holder and the metrology targets.
        xdiff = self.x_fibre_local - self.x_elbow_local
        ydiff = self.y_fibre_local - self.y_elbow_local
        
        half_width = padding + INS_POS_WIDTH2/2.0
        xdelta = half_width * ydiff / self.length2
        ydelta = half_width * xdiff / self.length2

        metrology_fraction = INS_POS_B1 / self.length2
        padding_fraction = padding / self.length2
        xmetrology_end = self.x_fibre_focal - \
            (metrology_fraction + padding_fraction) * xdiff
        ymetrology_end = self.y_fibre_focal - \
            (metrology_fraction + padding_fraction) * ydiff
        xmetrology_start = self.x_fibre_focal + padding_fraction * xdiff
        ymetrology_start = self.y_fibre_focal + padding_fraction * ydiff
            
        x1 = xmetrology_end - xdelta
        y1 = ymetrology_end + ydelta
        x2 = xmetrology_end + xdelta
        y2 = ymetrology_end - ydelta
        x3 = xmetrology_start + xdelta
        y3 = ymetrology_start - ydelta
        x4 = xmetrology_start - xdelta
        y4 = ymetrology_start + ydelta
        return ((x1,y1),(x2,y2),(x3,y3),(x4,y4))

    def avoidance_triangle_metrology(self, padding=0.0):
        """
        
        Return, in focal plane coordinates, the vertices of the
        triangle bounding the extension of the metrology target area,
        defined as ((x1, y2), (x2, y2), (x3, y3)).
        
        :Parameters:
        
        padding: float (optional)
            Safety padding added around the edge of the avoidance zone.
            Defaults to 0.0
                                
        :Returns:
        
        ((x1, y2), (x2, y2), (x3, y3))
        
        """
        # TODO: These avoidance zone functions need to be optimised.
        # TODO: REPLACE WITH AVOIDANCE_POLYGON_METROLOGY
        # Calculate the vertices of a triangle bounding the extension
        # of the metrology target area.
        xdiff = self.x_fibre_local - self.x_elbow_local
        ydiff = self.y_fibre_local - self.y_elbow_local
        
        half_width = padding + INS_POS_TW1/2.0
        xdelta = half_width * ydiff / self.length2
        ydelta = half_width * xdiff / self.length2

        metrology_fraction = INS_POS_B1 / self.length2
        padding_fraction = padding / self.length2
        xmetrology_end = self.x_fibre_focal - metrology_fraction * xdiff
        ymetrology_end = self.y_fibre_focal - metrology_fraction * ydiff

        triangle_fraction = (INS_POS_TB2+padding) / self.length2
        xtriangle_end = xmetrology_end - \
            (triangle_fraction+padding_fraction) * xdiff
        ytriangle_end = ymetrology_end - \
            (triangle_fraction+padding_fraction) * ydiff
            
        x1 = xmetrology_end - xdelta
        y1 = ymetrology_end + ydelta
        x2 = xmetrology_end + xdelta
        y2 = ymetrology_end - ydelta
        x3 = xtriangle_end
        y3 = ytriangle_end
        return ((x1,y1),(x2,y2),(x3,y3))

    def avoidance_ellipse_metrology(self, padding=0.0):
        """
        
        Return, in focal plane coordinates, the parameters of
        an ellipse covering the metrology avoidance zone.
        
        :Parameters:
        
        padding: float (optional)
            Safety padding added around the edge of the avoidance zone.
            Defaults to 0.0
                                
        :Returns:
        
        (xcen, ycen, major, minor, tilt)
            The centre, major axis, minor axis and tilt angle of the ellipse
        
        """
        # TODO: These avoidance zone functions need to be optimised.
        # Determine the extent of the metrology target zone.
        # Expand the zone slightly in the direction of the fibre holder.
        xdiff = self.x_fibre_local - self.x_elbow_local
        ydiff = self.y_fibre_local - self.y_elbow_local
        metrology_fraction = (INS_POS_B1+INS_POS_TB2+padding) / self.length2
        fibre_holder_fraction = 1.15 * (FIBRE_RADIUS+padding) / self.length2
        xmetrology_end = self.x_fibre_focal - metrology_fraction * xdiff
        ymetrology_end = self.y_fibre_focal - metrology_fraction * ydiff
        xmetrology_start = self.x_fibre_focal + fibre_holder_fraction * xdiff
        ymetrology_start = self.y_fibre_focal + fibre_holder_fraction * ydiff
        
        # The ellipse is placed in the middle of this zone
        xcen = (xmetrology_start + xmetrology_end) / 2.0
        ycen = (ymetrology_start + ymetrology_end) / 2.0
        
        # The major and minor axes are determined by the length and width
        # of the zone. Over-size the width a little to ensure the zone
        # is mostly covered by the ellipse.
        major = (FIBRE_RADIUS+INS_POS_B1+INS_POS_TB2)/2.0
        minor = 1.25 * INS_POS_WIDTH2 /  2.0
        tilt = math.atan2(ydiff, xdiff)
        
        return (xcen, ycen, major, minor, tilt)

    def avoidance_triangle_fibre(self, padding=0.0):
        """
        
        Return, in focal plane coordinates, the vertices of the
        triangle bounding the fibre avoidance area,
        defined as ((x1, y2), (x2, y2), (x3, y3)).
        
        :Parameters:
        
        padding: float (optional)
            Safety padding added around the edge of the avoidance zone.
            Defaults to 0.0
                                
        :Returns:
        
        ((x1, y2), (x2, y2), (x3, y3))
        
        """
        # TODO: These avoidance zone functions need to be optimised.
        # Calculate the vertices of a triangle bounding the the fibre
        # avoidance area.
        xdiff = self.x_fibre_local - self.x_elbow_local
        ydiff = self.y_fibre_local - self.y_elbow_local
        
        half_width = padding + INS_POS_TW2/2.0
        xdelta = half_width * ydiff / self.length2
        ydelta = half_width * xdiff / self.length2

        metrology_fraction = INS_POS_B1 / self.length2
        padding_fraction = padding / self.length2
        # The padding is extended at the pointed end of the triangle
        # and reduced at the straight end by this ratio.
        triang_length_ratio = INS_POS_B1 / INS_POS_TW2
        xmetrology_end = self.x_fibre_focal - metrology_fraction * xdiff
        ymetrology_end = self.y_fibre_focal - metrology_fraction * ydiff
        xtriangle_start = xmetrology_end + \
            padding_fraction * xdiff/triang_length_ratio
        ytriangle_start = ymetrology_end + \
            padding_fraction * ydiff/triang_length_ratio

        triangle_fraction = INS_POS_TB3 / self.length2
        xtriangle_end = xmetrology_end - \
            (triangle_fraction+padding_fraction*triang_length_ratio) * xdiff
        ytriangle_end = ymetrology_end - \
            (triangle_fraction+padding_fraction*triang_length_ratio) * ydiff
            
        x1 = xtriangle_start - xdelta
        y1 = ytriangle_start + ydelta
        x2 = xtriangle_start + xdelta
        y2 = ytriangle_start - ydelta
        x3 = xtriangle_end
        y3 = ytriangle_end
        return ((x1,y1),(x2,y2),(x3,y3))

    def avoidance_rectangles_beta(self, nzones=3, padding=0.0):
        """
         
        Return, in focal plane cooordinates, the vertices of the rectangles
        representing successive zones along the beta arm. These zones
        approximate the curved, 3-D design of the beta arm as a "stair case"
        of successive cuboids at different heights.
        
        Two beta arms will clash if two corresponding zones (at the same
        height) overlap.
         
        :Parameters:
         
        nzones: int (optional)
            The number of zones in which to divide the vulnerable
            portion of the beta arm. Default 3.
        padding: float (optional)
            Safety padding added around the edge of the avoidance zone.
            Defaults to 0.0
                    
        :Returns:
 
        list of rectangles: list of ((x1,y1),(x2,y2),(x3,y3),(x4,y4)))
            The beta arm avoidance rectangles, in focal plane coordinates.
         
        """
        # TODO: These avoidance zone functions need to be optimised.
        rect_list = []

        # Calculate X and Y increments for the corners of each rectangle
        # from the length and orientation of the beta arm.
        xdiff = self.x_fibre_local - self.x_elbow_local
        ydiff = self.y_fibre_local - self.y_elbow_local
        half_width = padding + INS_POS_WIDTH2/2.0
        xdelta = half_width * ydiff / self.length2
        ydelta = half_width * xdiff / self.length2
        
        # Find the start of the vulnerable zone on the beta arm.
        padding_fraction = padding / self.length2
        safe_fraction = self.bsafe / self.length2
        xbeta_start = self.x_centre_focal + self.x_elbow_local + \
                        (safe_fraction-padding_fraction) * xdiff
        ybeta_start = self.y_centre_focal + self.y_elbow_local + \
                        (safe_fraction-padding_fraction) * ydiff
         
        # Divide the vulnerable zone into nzones rectangles
        zone_fraction = (self.length2 + 2*padding - self.bsafe - \
                         INS_POS_B1 - INS_POS_B2) / \
                        ( nzones * self.length2)
        xfraction = xdiff * zone_fraction
        yfraction = ydiff * zone_fraction
        xzone_start = xbeta_start
        yzone_start = ybeta_start
        for zone in range(0, nzones):
            xzone_end = xzone_start + xfraction
            yzone_end = yzone_start + yfraction
            x1 = xzone_end - xdelta
            y1 = yzone_end + ydelta
            x2 = xzone_end + xdelta
            y2 = yzone_end - ydelta
            x3 = xzone_start + xdelta
            y3 = yzone_start - ydelta
            x4 = xzone_start - xdelta
            y4 = yzone_start + ydelta
            rect_list.append( ((x1,y1),(x2,y2),(x3,y3),(x4,y4)) )
            xzone_start = xzone_end
            yzone_start = yzone_end
         
        return rect_list
    
    def avoidance_ellipse_beta(self, padding=0.0):
        """
        
        Return, in focal plane coordinates, the parameters of
        an ellipse covering the beta arm avoidance zone.
        
        :Parameters:
        
        padding: float (optional)
            Safety padding added around the edge of the avoidance zone.
            Defaults to 0.0
                                
        :Returns:
        
        (xcen, ycen, major, minor, tilt)
            The centre, major axis, minor axis and tilt angle of the ellipse
        
        """
        # TODO: These avoidance zone functions need to be optimised.
        # Determine the extent of the beta arm zone.
        xdiff = self.x_fibre_local - self.x_elbow_local
        ydiff = self.y_fibre_local - self.y_elbow_local
        
        # Find the start of the vulnerable zone on the beta arm.
        # Expand the zone slightly in the direction of the beta axis.
        safe_fraction = 0.85 * (self.bsafe-padding) / self.length2
        xbeta_start = self.x_centre_focal + self.x_elbow_local + \
                        safe_fraction * xdiff
        ybeta_start = self.y_centre_focal + self.y_elbow_local + \
                        safe_fraction * ydiff

        zone_length = padding + self.length2 - self.bsafe - INS_POS_B1 - INS_POS_B2
        zone_fraction = zone_length / self.length2
        xbeta_end = xbeta_start + xdiff * zone_fraction
        ybeta_end = ybeta_start + ydiff * zone_fraction
        
        # The ellipse is placed in the middle of this zone
        xcen = (xbeta_start + xbeta_end) / 2.0
        ycen = (ybeta_start + ybeta_end) / 2.0
        
        # The major and minor axes are determined by the length and width
        # of the zone. Over-size the length a little to ensure the zone
        # is mostly covered by the ellipse.
        major = 1.1 * (self.length2 - self.bsafe - INS_POS_B1 - INS_POS_B2)/2.0
        minor = INS_POS_WIDTH2 /  2.0
        tilt = math.atan2(ydiff, xdiff)
        
        return (xcen, ycen, major, minor, tilt)
   
    def avoidance_zone_datum(self, padding=0.0):
        """
         
        Return, in focal plane cooordinates, the vertices of the rectangle
        connecting the circular "datum foot" to the axis of the beta arm.
        
        Two beta arms will clash if their datum rectangles overlap.
         
        :Parameters:
         
        padding: float (optional)
            Safety padding added around the edge of the avoidance zone.
            Defaults to 0.0
                    
        :Returns:
        
        (rect, circle)
            rect =((x1, y2), (x2, y2), (x3, y3), (x4, y4))
                The rectangular attachment to the "datum foot".
            circle = (xcen, ycen, radius)
                The centre and radius of the circular "datum foot" zone.
        
        """
        # TODO: These avoidance zone functions need to be optimised.
        shape_list = []

        # Calculate X and Y increments for the corners of each rectangle
        # from the length and orientation of the beta arm.
        xdiff = self.x_fibre_local - self.x_elbow_local
        ydiff = self.y_fibre_local - self.y_elbow_local
        half_width = padding + INS_POS_DW/2.0
        xdelta = half_width * ydiff / self.length2
        ydelta = half_width * xdiff / self.length2
        
        # The zone starts at the beta axis (elbow).
        padding_fraction = padding / self.length2
        xdatum_start = self.x_centre_focal + self.x_elbow_local - \
            padding_fraction * xdiff
        ydatum_start = self.y_centre_focal + self.y_elbow_local - \
            padding_fraction * ydiff
        
        # The zone ends at the datum switch above the alpha axis, so it is
        # the alpha arm length from the beta axis.
        #xdiff = self.x_fibre_local - self.x_elbow_local
        #ydiff = self.y_fibre_local - self.y_elbow_local
        datum_fraction = INS_POS_DL / self.length2
        xdatum_circ = self.x_centre_focal + self.x_elbow_local + datum_fraction * xdiff
        ydatum_circ = self.y_centre_focal + self.y_elbow_local + datum_fraction * ydiff
        xdatum_end = xdatum_circ + padding_fraction * xdiff
        ydatum_end = ydatum_circ + padding_fraction * ydiff

        x1 = xdatum_end-xdelta
        y1 = ydatum_end+ydelta
        x2 = xdatum_end+xdelta
        y2 = ydatum_end-ydelta
        x3 = xdatum_start+xdelta
        y3 = ydatum_start-ydelta
        x4 = xdatum_start-xdelta
        y4 = ydatum_start+ydelta
        rect = ((x1,y1),(x2,y2),(x3,y3),(x4,y4))
        shape_list.append( rect )

        # Assume the radius of the "foot" is the same as the
        # half-width of the beta arm.
        # radius = padding + INS_POS_DW/2.0
        # THE DATUM CIRCLE HAS GONE.  REPLACE FOR NOW WITH A TINY CIRCLE.
        radius = 0.1
        circle = (xdatum_circ, ydatum_circ, radius)
        shape_list.append( circle )
        return shape_list
         
    # ***
    def in_conflict_with(self, other, safety=INS_POS_SAFETY, symmetrical=True):
        """
        
        Determine whether this positioner is in conflict with another,
        given their current targets.
        
        :Parameters:
        
        other: Positioner object
            The other positioner object with which conflict is to be checked.
        safety: float (optional)
            All conflict avoidance zones are padded with this safety
            tolerance. The default is the design tolerance, INS_POS_SAFETY.
        symmetrical: boolean (optional)
            Set to True to make all conflict checks symmetrical, so that when
            one positioner is in conflict both are in conflict. This only makes
            a difference to the fibre snag check.
            The default is True - when there is a fibre snag, both positioners
            are considered to be in conflict.
   
        :Returns:

        in_conflict: boolean
            True if the two positioners are in conflict
            False if the positioners are not in conflict
                
        """
        assert isinstance(other, Positioner)
        
        # Determine the distance between the fibre targets associated with
        # the two positioners.
        xdiff = self.x_fibre_focal - other.x_fibre_focal
        ydiff = self.y_fibre_focal - other.y_fibre_focal
        distsq = xdiff * xdiff + ydiff * ydiff
        
        # Are the targets sufficiently far apart (further than the maximum
        # reach of this positioner plus minimum separation) that a conflict
        # is not possible?
        if distsq > self._safedistsq:
            return False

        # Are the targets closer than the minimum separation distance
        # between fibre holders?
        mindistsq = (INS_POS_MINDIST+2*safety) * (INS_POS_MINDIST+2*safety)
        if distsq <= mindistsq:
            self.in_conflict = True
            self.conflict_type = self.CONFLICT_TOOCLOSE
            self.conflict_reason = "Conflict between %s and %s: " % \
                (self.name, other.name)
            self.conflict_reason += "Targets too close: dist^2 %.3f <= %.3f" % \
                (distsq, mindistsq)
            return True
       
        # Do the rectangular metrology zones of the two positioners overlap?
        # TODO: REPLACE WITH POLYGONAL ZONE - SEE BELOW
        this_metrology_rectangle = \
            self.avoidance_rectangle_metrology(padding=safety)
        other_metrology_rectangle = \
            other.avoidance_rectangle_metrology(padding=safety)
        logger.debug("My metrology rectangle: " + str(this_metrology_rectangle))
        logger.debug("Other metrology rectangle: " + str(other_metrology_rectangle))
        if util.quadrangles_intersect(this_metrology_rectangle,
                                      other_metrology_rectangle):
            self.in_conflict = True
            self.conflict_type = self.CONFLICT_METROLOGY
            self.conflict_reason = "Conflict between %s and %s: " % \
                (self.name, other.name)
            self.conflict_reason += "Metrology rectangles intersect"
            return True
#         # REPLACE WITH THE FOLLOWING
#         this_metrology_polygon = self.avoidance_polygon_metrology()
#         other_polygon_polygon = other.avoidance_polygon_metrology()
#         if util.polygons_intersect(this_metrology_polygon, other_polygon_polygon):
#             self.in_conflict = True
#             self.conflict_type = self.CONFLICT_METROLOGY
#             self.conflict_reason = "Metrology polygons intersect"
#             return True

        # Does one of the beta arm zones overlap with the same zone on
        # the other arm?
        if INCLUDE_BETA_ZONE:
            this_rectangles_beta = \
                self.avoidance_rectangles_beta(padding=safety)
            other_rectangles_beta = \
                other.avoidance_rectangles_beta(padding=safety)
            # Test the rectangles one at a time.
            zone = 1
            for this_rect, other_rect in zip(this_rectangles_beta,
                                         other_rectangles_beta):
                logger.debug("My beta rectangle " + str(zone) + ": " + str(this_rect))
                logger.debug("Other beta rectangle " + str(zone) + ": " + str(other_rect))
                if util.quadrangles_intersect(this_rect, other_rect):
                    self.conflict_type = self.CONFLICT_BETA
                    self.in_conflict = True
                    self.conflict_reason = "Conflict between %s and %s: " % \
                        (self.name, other.name)
                    self.conflict_reason += \
                        "Beta arms intersect at same level (zone %d)" % zone
#                     if GRAPHICAL_DEBUGGING and plotting is not None:
#                         (xlist1, ylist1) = util.polygon_to_points(this_rect)
#                         (xlist2, ylist2) = util.polygon_to_points(other_rect)
#                         plotaxis = plotting.plot_xy(xlist1, ylist1,
#                                     title="Beta zone conflict", showplot=False)
#                         plotting.plot_xy(xlist2, ylist2, plotaxis=plotaxis,
#                                     showplot=True)
                    return True
                zone = zone + 1

        if INCLUDE_DATUM_FOOT:
            (this_datum_rect, this_datum_circle) = \
                self.avoidance_zone_datum(padding=safety)
            (other_datum_rect, other_datum_circle) = \
                other.avoidance_zone_datum(padding=safety)
            logger.debug("My datum circle: " + str(this_datum_circle))
            logger.debug("Other datum circle: " + str(other_datum_circle))
            limitsq = this_datum_circle[2] * this_datum_circle[2] * 4.0
            if util.closer_than(this_datum_circle[0], this_datum_circle[1],
                                other_datum_circle[0], other_datum_circle[1],
                                limitsq):
                self.in_conflict = True
                self.conflict_type = self.CONFLICT_DATUM
                self.conflict_reason = "Conflict between %s and %s: " % \
                        (self.name, other.name)
                self.conflict_reason += \
                        "Datum actuators intersect"
#                 if GRAPHICAL_DEBUGGING and plotting is not None:
#                     xcen1 = [this_datum_circle[0]]
#                     ycen1 = [this_datum_circle[1]]
#                     rad1 = this_datum_circle[2]
#                     plotaxis = plotting.plot_circles( xcen1, ycen1,
#                         rad1, title="Datum actuators intersect", showplot=False )
#                     xcen2 = [other_datum_circle[0]]
#                     ycen2 = [other_datum_circle[1]]
#                     rad2 = other_datum_circle[2]
#                     plotaxis = plotting.plot_circles( xcen2, ycen2,
#                         rad2, showplot=True )
                return True
                
            logger.debug("My datum rect: " + str(this_datum_rect))
            logger.debug("Other datum rect: " + str(other_datum_rect))
            if util.quadrangles_intersect(this_datum_rect, other_datum_rect):
                self.in_conflict = True
                self.conflict_type = self.CONFLICT_DATUM
                self.conflict_reason = "Conflict between %s and %s: " % \
                        (self.name, other.name)
                self.conflict_reason += \
                        "Datum actuator rectangles intersect"
#                 if GRAPHICAL_DEBUGGING and plotting is not None:
#                     (xlist1, ylist1) = util.polygon_to_points(this_datum_rect)
#                     (xlist2, ylist2) = util.polygon_to_points(other_datum_rect)
#                     plotaxis = plotting.plot_xy(xlist1, ylist1,
#                                     title="Datum actuator rectangles conflict",
#                                     showplot=False)
#                     plotting.plot_xy(xlist2, ylist2, plotaxis=plotaxis,
#                                     showplot=True)
                return True

        # Does the target lie within the other positioner's fibre avoidance triangle?
        other_fibre_triangle = other.avoidance_triangle_fibre(padding=safety)
        logger.debug("Other fibre triangle: " + str(other_fibre_triangle))
        if util.point_inside_triangle(self.x_fibre_focal, self.y_fibre_focal,
                    other_fibre_triangle[0][0], other_fibre_triangle[0][1],
                    other_fibre_triangle[1][0], other_fibre_triangle[1][1],
                    other_fibre_triangle[2][0], other_fibre_triangle[2][1]):
            self.in_conflict = True
            self.conflict_type = self.CONFLICT_FIBRE_SNAG
            self.conflict_reason = "Conflict between %s and %s: " % \
                (self.name, other.name)
            self.conflict_reason += "Target inside fibre snag triangle"
            logger.debug(self.conflict_reason)
            return True
        # Does the fibre holder intersect with the other positioner's avoidance
        # triangle (within the given safety tolerance)?
        if util.triangle_intersects_circle(other_fibre_triangle,
                                           self.x_fibre_focal,
                                           self.y_fibre_focal,
                                           (FIBRE_RADIUS+safety)):
            self.in_conflict = True
            self.conflict_type = self.CONFLICT_FIBRE_SNAG
            self.conflict_reason = "Conflict between %s and %s: " % \
                (self.name, other.name)
            self.conflict_reason += "Fibre holder intersects fibre snag triangle"
            logger.debug(self.conflict_reason)
#             if GRAPHICAL_DEBUGGING and plotting is not None:
#                 (xlist, ylist) = util.polygon_to_points(other_fibre_triangle)
#                 plotaxis = plotting.plot_circles([self.x_fibre_focal],
#                     [self.y_fibre_focal], FIBRE_RADIUS,
#                     title="Fibre snag zone conflict", showplot=False)
#                 plotting.plot_xy(xlist, ylist, plotaxis=plotaxis, showplot=True)
#                 plotting.close()
            return True
        
        if symmetrical:
            # If a symmetrical check is needed, also test whether the other
            # positioner's target lies within this positioner's fibre snag zone.
            this_fibre_triangle = self.avoidance_triangle_fibre(padding=safety)
            logger.debug("This fibre triangle: " + str(this_fibre_triangle))
            if util.point_inside_triangle(other.x_fibre_focal, other.y_fibre_focal,
                    this_fibre_triangle[0][0], this_fibre_triangle[0][1],
                    this_fibre_triangle[1][0], this_fibre_triangle[1][1],
                    this_fibre_triangle[2][0], this_fibre_triangle[2][1]):
                self.in_conflict = True
                self.conflict_type = self.CONFLICT_FIBRE_SNAG
                self.conflict_reason = "Conflict between %s and %s: " % \
                    (self.name, other.name)
                self.conflict_reason += "Other's target inside fibre snag triangle"
                logger.debug(self.conflict_reason)
                return True
            if util.triangle_intersects_circle(this_fibre_triangle,
                                           other.x_fibre_focal,
                                           other.y_fibre_focal,
                                           (FIBRE_RADIUS+safety)):
                self.in_conflict = True
                self.conflict_type = self.CONFLICT_FIBRE_SNAG
                self.conflict_reason = "Conflict between %s and %s: " % \
                (self.name, other.name)
                self.conflict_reason += "Other's fibre holder inside fibre snag triangle"
                logger.debug(self.conflict_reason)
#                 if GRAPHICAL_DEBUGGING and plotting is not None:
#                     (xlist, ylist) = util.polygon_to_points(this_fibre_triangle)
#                     plotaxis = plotting.plot_circles([other.x_fibre_focal],
#                         [other.y_fibre_focal], FIBRE_RADIUS,
#                         title="Fibre snag zone conflict in reverse", showplot=False)
#                     plotting.plot_xy(xlist, ylist, plotaxis=plotaxis, showplot=True)
#                     plotting.close()
                return True

        return False

    # ***
    def in_conflict_with_neighbours(self, test_all=True, safety=INS_POS_SAFETY,
                                    symmetrical=True):
        """
        
        Determine whether this positioner is in conflict with any of its
        neighbours, given their current targets.
        
        :Parameters:
        
        test_all: bool (optional)
            Set to True (the default) to test all the neighbours.
            If False, the function returns after the first conflict
            is detected.
        safety: float (optional)
            All conflict avoidance zones are padded with this safety
            tolerance. The default is the design tolerance, INS_POS_SAFETY.
        symmetrical: boolean (optional)
            Set to True to make all conflict checks symmetrical, so that when
            one positioner is in conflict both are in conflict. This only makes
            a difference to the fibre snag check.
            The default is True - when there is a fibre snag, both positioners
            are considered to be in conflict.
   
        :Returns:

        in_conflict: boolean
            True if the positioner is in conflict with any of its neighbours
            False if the positioner is not in conflict with all neighbours.
        
        """
        in_conflict = False
        for neighbour in self.near_neighbours:
            logger.debug("Check for conflict between %s and near neighbour %s" % \
                          (self.name, neighbour.name))
            if self.in_conflict_with(neighbour, safety=safety,
                                     symmetrical=symmetrical):
                in_conflict = True
                if GRAPHICAL_DEBUGGING:
                    title = "%s in conflict with %s" % (self.name,
                                                        neighbour.name)
                    title += "\n" + self.conflict_reason
                    plotfig = self.plot(showplot=False, description=title)
                    neighbour.plot(plotfig=plotfig)
                if not test_all:
                    break
#             else:
#                 if GRAPHICAL_DEBUGGING:
#                     title = "%s NOT in conflict with %s" % (self.name,
#                                                             neighbour.name)
#                     plotfig = self.plot(showplot=False, description=title)
#                     neighbour.plot(plotfig=plotfig)
                
        for neighbour in self.far_neighbours:
            logger.debug("Check for conflict between %s and far neighbour %s" % \
                          (self.name, neighbour.name))
            # A far neighbour can only conflict when the positioner's
            # reach is outside a defined limit.
            if self.r_fibre_local >= (self.limitdist - 2.0*safety):
                if self.in_conflict_with(neighbour, safety=safety,
                                         symmetrical=symmetrical):
                    in_conflict = True
                    if GRAPHICAL_DEBUGGING:
                        title = "%s in conflict with %s" % (self.name,
                                                            neighbour.name)
                        title += "\n" + self.conflict_reason
                        plotfig = self.plot(showplot=False, description=title)
                        neighbour.plot(plotfig=plotfig)
                    if not test_all:
                        break
#                 else:
#                     if GRAPHICAL_DEBUGGING:
#                         title = "%s NOT in conflict with %s" % (self.name,
#                                                                 neighbour.name)
#                         plotfig = self.plot(showplot=False, description=title)
#                         neighbour.plot(plotfig=plotfig)
        return in_conflict

    def plot(self, description='', plotfig=None, showplot=True):
        """
        
        Plot the configuration of the positioner.
        
        NOTE: This function requires access to the plotting
        module, plotting.py. If the module can't be imported
        the function returns with an apologetic message.
        
        :Parameters:
        
        description: str, optional
            Optional description to be added to the positioner plot.
        plotfig: matplotlib Figure object, optional
            An existing matplotlib figure onto which to add the plot.
            If not given, a new figure is created.
        showplot: boolean, optional
            If True (the default) show the plot when finished.
            Set of False if you are overlaying several plots, and
            this plot is not the last.
            
        :Returns:
        
        plotfig: matplotlib Figure object
            A matplotlib figure containing the plot.
        
        """
        if plotting is None:
            logger.warning("Plotting is disabled")
            return
        
        # Create a new plotting figure, if necessary.
        if plotfig is None:
            plotfig = plotting.new_figure(1, figsize=(10,9), stitle='')

        # First plot the centres of the positioners on the grid
        # using "O" symbols.
        xcentres = [self.x_centre_focal]
        ycentres = [self.y_centre_focal]

        # Plot the centre circles showing the patrol field.
        # Inner boundary - black dashed line.
        plotaxis = plotting.plot_circles( xcentres, ycentres,
                        self.inner, plotfig=plotfig,
                        facecolor='w', circcolor='k ', linefmt='ko',
                        linestyle='--', linewidth=0.3,
                        xlabel="X (mm)", ylabel="Y (mm)", title=description,
                        grid=True, showplot=False )
        # Nominally "Safe" boundary - blue dashed line.
        plotaxis = plotting.plot_circles( xcentres, ycentres,
                        self.outer/2.0,
                        plotfig=plotfig, plotaxis=plotaxis,
                        facecolor='w', circcolor='b ', linefmt='b ',
                        linestyle='--', linewidth=0.2, grid=False,
                        showplot=False )
        # Outer boundary - black dashed line.
        plotaxis = plotting.plot_circles( xcentres, ycentres,
                        self.outer, plotfig=plotfig, plotaxis=plotaxis,
                        facecolor='w', circcolor='k ', linefmt='k ',
                        linestyle='--', linewidth=0.3, grid=False,
                        showplot=False )
        
        # Add a label giving the positioner ID
        xtxt = self.x_centre_focal - self.length1/3.0
        ytxt = self.y_centre_focal - self.length1/3.0
        plotaxis = plotting.plot_text( str(self.ident), xtxt, ytxt,
                                       plotfig=plotfig, plotaxis=plotaxis,
                                       showplot=False)

        # Change the colour of the positioner arms depending on whether there
        # is a target assigned and whether the positioner is in conflict.
        if self.target_assigned:
            target_fmt = 'b+'   # Blue + - assigned
            if self.in_conflict:
                arm_fmt = 'r '  # Red - in conflict
            else:
                arm_fmt = 'k '  # Black - target assigned 
        else:
            target_fmt = 'gx'   # Green x - not assigned
            if self.in_conflict:
                arm_fmt = 'r '  # Red - in conflict
            elif self.conflict_type == self.CONFLICT_UNREACHABLE:
                arm_fmt = 'm '  # Magenta - attempted to assign unreachable target                
            else:
                arm_fmt = 'g '  # Green - target not assigned

        # Add to the plot the location of all the fibres using the
        # symbol defined in target_fmt.
        xtargets = [self.x_fibre_focal]
        ytargets = [self.y_fibre_focal]
        plotaxis = plotting.plot_xy( xtargets, ytargets, plotfig=plotfig,
                                     plotaxis=plotaxis,
                                     linefmt=target_fmt, linestyle=' ',
                                     showplot=False )
        # Show the fibre approach tolerance as a thick red circle.
        plotaxis = plotting.plot_circles( xtargets, ytargets,
                        FIBRE_RADIUS,
                        plotfig=plotfig, plotaxis=plotaxis,
                        facecolor='w', circcolor='r ', linefmt='r ',
                        linestyle='-', linewidth=3, grid=False,
                        showplot=False )
        if PLOT_PADDING:
            plotaxis = plotting.plot_circles( xtargets, ytargets,
                            FIBRE_RADIUS+INS_POS_SAFETY,
                            plotfig=plotfig, plotaxis=plotaxis,
                            facecolor='w', circcolor='r ', linefmt='r ',
                            linestyle='--', linewidth=3, grid=False,
                            showplot=False )
        
        # Show the orientation of the alpha and beta arms as thick lines.
        x_elbow_focal = self.x_centre_focal + self.x_elbow_local
        y_elbow_focal = self.y_centre_focal + self.y_elbow_local
        xlist = [self.x_centre_focal, x_elbow_focal, self.x_fibre_focal]
        ylist = [self.y_centre_focal, y_elbow_focal, self.y_fibre_focal]
        plotaxis = plotting.plot_xy( xlist, ylist, plotfig=plotfig,
                    plotaxis=plotaxis, linefmt=arm_fmt, linestyle='-',
                    linewidth=6, showplot=False )

        # Add thick lined polygons to show the avoidance zones. Note that
        # avoidance zones are plotted without the safety tolerance.
        # (1a) The metrology target zone - red.
        poly1 = self.avoidance_polygon_metrology()
        (xlist, ylist) = util.polygon_to_points(poly1)
        plotaxis = plotting.plot_xy( xlist, ylist, plotfig=plotfig,
                    plotaxis=plotaxis, linefmt='r ', linestyle='-',
                    linewidth=3, showplot=False )
        if PLOT_PADDING:
            poly1_padded = \
                self.avoidance_polygon_metrology(padding=INS_POS_SAFETY)
            (xlist, ylist) = util.polygon_to_points(poly1_padded)
            plotaxis = plotting.plot_xy( xlist, ylist, plotfig=plotfig,
                            plotaxis=plotaxis, linefmt='r ', linestyle='--',
                            linewidth=3, showplot=False )
        
        # (1b) Cover the metrology zone with an ellipse.
        # FIXME: The ellipse becomes squashed at certain locations and tilts.
        if PLOT_ELLIPSES:
            (xcen, ycen, major, minor, tilt) = \
                self.avoidance_ellipse_metrology()
            plotaxis = plotting.plot_ellipses( [xcen], [ycen], major, minor,
                            tilt, plotfig=plotfig, plotaxis=plotaxis,
                            facecolor='w', ellipsecolor='k ', linefmt='r ',
                            linestyle='-', linewidth=3, grid=False,
                            showplot=False )

        # (2) The fibre snag avoidance zone - yellow.
        triangle2 = self.avoidance_triangle_fibre()
        (xlist, ylist) = util.polygon_to_points(triangle2)
        plotaxis = plotting.plot_xy( xlist, ylist, plotfig=plotfig,
                    plotaxis=plotaxis, linefmt='y ', linestyle='-',
                    linewidth=3, showplot=False )
#         if PLOT_PADDING:
#             triangle2_padded = \
#                 self.avoidance_triangle_fibre(padding=INS_POS_SAFETY)
#             (xlist, ylist) = util.polygon_to_points(triangle2_padded)
#             plotaxis = plotting.plot_xy( xlist, ylist, plotfig=plotfig,
#                             plotaxis=plotaxis, linefmt='y ', linestyle='--',
#                             linewidth=3, showplot=False )

        # (3a) Add a thick dotted rectangles showing the beta arm
        # collision zones - magenta.
        if INCLUDE_BETA_ZONE:
            rectangles_beta = self.avoidance_rectangles_beta()
            for rect in rectangles_beta:
        
                (xlist, ylist) = util.polygon_to_points(rect)
                plotaxis = plotting.plot_xy( xlist, ylist, plotfig=plotfig,
                                plotaxis=plotaxis, linefmt='m ', linestyle='-.',
                                linewidth=3, showplot=False )
            if PLOT_PADDING:
                rectangles_beta_padded = \
                    self.avoidance_rectangles_beta(padding=INS_POS_SAFETY)
                for rect in rectangles_beta_padded:
        
                    (xlist, ylist) = util.polygon_to_points(rect)
                    plotaxis = plotting.plot_xy( xlist, ylist, plotfig=plotfig,
                                    plotaxis=plotaxis, linefmt='m ',
                                    linestyle=':', linewidth=3,
                                    showplot=False )

        # (3b) Cover the beta zone with an ellipse.
        # FIXME: The ellipse becomes squashed at certain locations and tilts.
        if PLOT_ELLIPSES:
            (xcen, ycen, major, minor, tilt) = \
                self.avoidance_ellipse_beta()
            plotaxis = plotting.plot_ellipses( [xcen], [ycen], major, minor,
                            tilt, plotfig=plotfig, plotaxis=plotaxis,
                            facecolor='w', ellipsecolor='b ', linefmt='r ',
                            linestyle='-', linewidth=3, grid=False,
                            showplot=False )

        # (4) Add a thick dotted circle showing the location of the
        # "datum foot" of the beta arm - cyan.
        if INCLUDE_DATUM_FOOT:
            (rect4, circle4) = self.avoidance_zone_datum()
            plotaxis = plotting.plot_circles( [circle4[0]], [circle4[1]],
                            circle4[2], plotfig=plotfig, plotaxis=plotaxis, 
                            facecolor='w', circcolor='c ', linefmt='co',
                            linestyle=':', linewidth=3,
                            xlabel="X (mm)", ylabel="Y (mm)", title=description,
                            grid=True, showplot=False )
            (xlist, ylist) = util.polygon_to_points(rect4)
            plotaxis = plotting.plot_xy( xlist, ylist, plotfig=plotfig,
                            plotaxis=plotaxis, linefmt='c ', linestyle=':',
                            linewidth=3, showplot=False )
            if PLOT_PADDING:
                (rect4_pad, circle4_pad) = \
                    self.avoidance_zone_datum(padding=INS_POS_SAFETY)
                plotaxis = plotting.plot_circles( [circle4_pad[0]],
                            [circle4_pad[1]], circle4_pad[2], plotfig=plotfig,
                            plotaxis=plotaxis, 
                            facecolor='w', circcolor='c ', linefmt='co',
                            linestyle=':', linewidth=2,
                            xlabel="X (mm)", ylabel="Y (mm)", title=description,
                            grid=True, showplot=False )
                (xlist, ylist) = util.polygon_to_points(rect4_pad)
                plotaxis = plotting.plot_xy( xlist, ylist, plotfig=plotfig,
                            plotaxis=plotaxis, linefmt='c ', linestyle=':',
                            linewidth=2, showplot=False )
        if showplot:
            plotting.show_plot()
            plotting.close()
        return plotfig


    # FROM PATH ANALYSIS RELEASE - 13 JUN 2016
    # FIXME: REWORK. REPLACED BY SET_ARM_ANGLES AND GET_ARM_ANGLES?

    def motor_position_to_focal_cartesian(self, alpha_position, beta_position):
        x_fibre_focal = self.x_centre_focal + self.length1 * math.cos(alpha_position) + \
                        self.length2 * math.cos(alpha_position + beta_position + math.pi)
        y_fibre_focal = self.y_centre_focal + self.length1 * math.sin(alpha_position) + \
                        self.length2 * math.sin(alpha_position + beta_position + math.pi)

        return (x_fibre_focal, y_fibre_focal)

    # Add motors
    # TODO: REWORK. read the maximum and minimum speed from an interface with hardware
    def add_motors(self, sim_length, alfa_position_target, beta_position_target,\
                   max_speed_alpha,min_speed_alpha,max_speed_beta,min_speed_beta):
        self.motor1 = Motor( sim_length, alfa_position_target,max_speed_alpha,min_speed_alpha )
        self.motor2 = Motor( sim_length, beta_position_target,max_speed_beta,min_speed_beta)

class PositionerGrid(object):
    """
    
    Class representing a hexagonal grid of fibre positioners.
    
    Each grid is built from a collection of fibre positioners at
    different MOONS focal plane locations and orientations. Column and
    row numbers in the grid are also given explicitly to disambiguate
    the cells.
    
    The defined locations and orientations are used to acquire targets
    and check for conflicts between fibre positioners.
    
    The defined column and row numbers are used to determine which fibre
    positioners are neighbours.
    
    """
    # Class variables defining the location of information within the
    # configuration list.
    C_IDENT = 0
    C_RCEN = 1
    C_THCEN = 2
    C_ORIENT = 3
    C_COLUMN = 4
    C_ROW = 5
    C_SIMULATED = 6
    C_LOCKED = 7
    C_RFIBRE = 8
    C_THFIBRE = 9
    C_PFIBRE = 10
   
    def __init__(self, config_list, length1=INS_POS_LENGTH1,
                 length2=INS_POS_LENGTH2):
        """
        
        Constructor
        
        :Parameters:
        
        config_list: list of (ident, rcen, thcen, orient, column, row,
            [simulated, locked, rfibre, thfibre, pfibre])
            List of configuration parameters describing the IDs,
            locations, orientations, column numbers and row numbers
            of the fibre positioners. Each list may also contain
            a simulation flag, positioner locked (broken) flag and
            a default location for the fibre.
        length1: float (optional)
            Length of all alpha arms.
            Defaults to the designed length, INS_POS_LENGTH1
        length2: float (optional)
            Length of all beta arms.
            Defaults to the designed length, INS_POS_LENGTH2

        """
        assert isinstance(config_list, (tuple,list))
        assert isinstance(config_list[0], (tuple,list))
        assert len(config_list) > 1
        assert len(config_list[0]) >= 6
            
        # A list in which to store the positioners and a dictionary
        # in which to translate ident number into row and column.
        self.grid = []
        self.positioner_count = 0
        self.ident_dict = {}
        
        # The pitch size of the grid is assumed to be the same as the reach of
        # the positioners.
        self.pitch = length1 + length2
        
        # Determine the number of columns and rows in the grid.
        cmax = config_list[0][self.C_COLUMN]
        rmax = config_list[0][self.C_ROW]
        idmax = config_list[0][self.C_IDENT]
        for config in config_list:
            if config[self.C_COLUMN] > cmax:
                cmax = config[self.C_COLUMN]
            if config[self.C_ROW] > rmax:
                rmax = config[self.C_ROW]
            if config[self.C_IDENT] > idmax:
                idmax = config[self.C_IDENT]
        # Use the range of column and row numbers given
        self.columns = cmax + 1
        self.rows = rmax + 1
        pmax = idmax + 1
        strg = "Fibre positioner grid contains %d positioners " % pmax
        strg += "arranged in up to %d columns and %d rows" % \
            (self.columns, self.rows)
        logger.info(strg)

        # Create an empty list in which to store the positioners
        self.positioners = [None] * pmax
        
        # Create an empty grid of the correct size.
        # TODO: Creating a rectangular grid to cover a circular focal plane
        # leaves several blank locations. Is there a better way of addressing
        # positioners but still allowing neighbours to be quickly identified?
        for row in range(0, self.rows):
            grid_rows = []
            for column in range(0, self.columns):
                grid_rows.append(None)
            self.grid.append(grid_rows)

        # Populate the grid with positioners. column and row are assumed
        # to have been already derived correctly from the centre of each
        # positioner.
        for config in config_list:
            ident = config[self.C_IDENT]
            rcen = config[self.C_RCEN]
            thcen = config[self.C_THCEN]
            orient = config[self.C_ORIENT]
            column = config[self.C_COLUMN]
            row = config[self.C_ROW]
            if self.positioners[ident] is None:
                if self.grid[row][column] is None:
                    # TODO: Locked (broken) positioners not taken into account.
                    positioner = Positioner(ident, rcen, thcen, orient, column, row,
                                            simulated=False, locked=False,
                                            length1=length1, length2=length2)
                    self.positioners[ident] = positioner
#                     print("Adding positioner %d at (%d,%d)" % (ident, row, column))
                    self.grid[row][column] = ident
                    self.ident_dict[str(ident)] = (column,row) # NEEDED?
                    self.positioner_count += 1
                    logger.debug("Positioner %s added at row=%d, column=%d" % \
                             (positioner.name, row, column) )
                else:
                    # Two or more positioners defined at the same column and row
                    strg = "More than one positioner found at row=%d, column=%d\n" % \
                        (row, column)
                    strg += "\tAttempt to add positioner "
                    strg += "ident=%d, rcen=%f, thcen=%f\n" % (ident, rcen, thcen)
                    strg += "\tSpace occupied by "
                    ident = self.grid[row][column]
                    other = self.positioners[ident]
                    strg += "ident=%d, rcen=%f, thcen=%f" % \
                        (other.ident, other.r_centre_focal, other.theta_centre_focal)
                    raise AttributeError(strg)
            else:
                # Two or more positioners defined with the same identifier
                strg = "More than one positioner with identifier %d\n" % ident
                strg += "\tAttempt to add positioner "
                strg += "rcen=%f, thcen=%f\n" % (rcen, thcen)
                strg += "\tSpace occupied by "
                other = self.positioners[ident]
                strg += "rcen=%f, thcen=%f" % \
                    (other.r_centre_focal, other.theta_centre_focal)
                raise AttributeError(strg)
           
        # Set up the references between each positioner and its neighbours.
        self._define_neighbours()
        
        # Initialise the conflict counters.
        self.init_counters()
 
    def __del__(self):
        """
        
        Destructor
        
        """
        try:
            for positioner in self.positioners:
                if positioner is not None:
                    del positioner
        except:
            pass

    def __str__(self):
        """
        
        Return a readable string describing the positioner grid.
        
        """
        strg = '=' * 75 + '\n'
        strg += "Grid of %d positioners in %d cols x %d rows with pitch=%.3f:\n" % \
            (self.positioner_count, self.columns, self.rows, self.pitch)
        for positioner in self.positioners:
            if positioner is not None:
                strg += '-' * 75 + '\n'
                strg += str(positioner) + '\n'
        return strg


    def _define_neighbours(self):
        """
        
        Helper function which defines the neighbours surrounding
        each positioner in the grid. Called on initialisation.
                
        :Parameters:
        
        None.
        
        """
        # Near neighbours are within an inner ring adjacent to each
        # positioner. Far neighbours are within the second ring adjacent
        # to the near neighbours. Near and far neighbours both overlap
        # each positioner's patrol field, but far neighbours only need to
        # be considered when a positioner is reaching beyond a certain limit.
        for row in range(0, self.rows):
            # The column offsets are slightly different for odd and even
            # numbered rows in the grid.
            if row % 2 == 0:
                # Even numbered row
                for column in range(0, self.columns):
                    positioner = self.get_row_column(row, column)
                    if positioner is not None:
                        neighbour_offsets = [ \
                            (row-2, column-1, False),
                            (row-2, column, False),
                            (row-2, column+1, False),
                            (row-1, column-1, False),
                            (row-1, column, True),
                            (row-1, column+1, True),
                            (row-1, column+2, False),
                            (row, column-2, False),
                            (row, column-1, True),
                            (row, column+1, True),
                            (row, column+2, False),
                            (row+1, column-1, False),
                            (row+1, column, True),
                            (row+1, column+1, True),
                            (row+1, column+2, False),
                            (row+2, column-1, False),
                            (row+2, column, False),
                            (row+2, column+1, False),
                            ]
                        # Add neighbours at the specified row and column offsets.
                        for (roff, coff, near) in neighbour_offsets:
                            neighbour = self.get_row_column(roff, coff)
                            if neighbour is not None:
                                positioner.add_neighbour(neighbour, near=near)
                                    
                        if logger.getEffectiveLevel() == logging.DEBUG:
                            strg = "Even row. Positioner %s has the following neighbours." % \
                                positioner.name
                            strg += "\nNEAR neighbours: "
                            for neighbour in positioner.near_neighbours:
                                strg += "%s " % neighbour.name
                            strg += "\nFAR neighbours: "
                            for neighbour in positioner.far_neighbours:
                                strg += "%s " % neighbour.name
                            logger.debug(strg)
            else:
                # Odd numbered row
                for column in range(0, self.columns):
                    positioner = self.get_row_column(row, column)
                    if positioner is not None:
                        neighbour_offsets = [ \
                            (row-2, column-1, False),
                            (row-2, column, False),
                            (row-2, column+1, False),
                            (row-1, column-2, False),
                            (row-1, column-1, True),
                            (row-1, column, True),
                            (row-1, column+1, False),
                            (row, column-2, False),
                            (row, column-1, True),
                            (row, column+1, True),
                            (row, column+2, False),
                            (row+1, column-2, False),
                            (row+1, column-1, True),
                            (row+1, column, True),
                            (row+1, column+1, False),
                            (row+2, column-1, False),
                            (row+2, column, False),
                            (row+2, column+1, False),
                            ]
                        # Add neighbours at the specified row and column offsets.
                        for (roff, coff, near) in neighbour_offsets:
                            neighbour = self.get_row_column(roff, coff)
                            if neighbour is not None:
                                positioner.add_neighbour(neighbour, near=near)
                        if logger.getEffectiveLevel() == logging.DEBUG:
                            strg = "Odd row. Positioner %s has the following neighbours." % \
                                positioner.name
                            strg += "\nNEAR neighbours: "
                            for neighbour in positioner.near_neighbours:
                                strg += "%s " % neighbour.name
                            strg += "\nFAR neighbours: "
                            for neighbour in positioner.far_neighbours:
                                strg += "%s " % neighbour.name
                            logger.debug(strg)

    def get_positioner(self, ident):
        """
        
        Return the positioner in the grid with the given ID.
                
        :Parameters:
        
        ident: int
            Positioner ID.
            
        :Returns:
        
        positioner: Positioner
            The positioner object with the given ID, or None if the positioner
            cannot be found.
            
        """
        if ident >= 0 and ident < len(self.positioners):
            return self.positioners[ident]
        else:
            return None

    def get_row_column(self, row, column):
        """
        
        Return the positioner in the grid at the given row and column location.
                
        :Parameters:
        
        row: int
            Row number in the grid.
        column: int
            Column number in the grid.
            
        :Returns:
        
        positioner: Positioner
            The positioner object with the given location, or None if the
            positioner cannot be found.
            
        """
        if row >= 0 and row < self.rows and column >= 0 and column < self.columns:
            ident = self.grid[row][column]
            if ident is not None:
                return self.positioners[ident]
            else:
                return None
        else:
            return None

    def plot(self, description='', targetlist=[], withcircles=True, simple=True,
             trivial=False, showplot=True):
        """
        
        Plot all the positioners in the grid.
        
        :Parameters:
        
        description: str (optional)
            Optional description to be added to the plot title.
        targetlist: list of Target objects (optional)
            Optional list of targets whose positions are to be
            shown on the plot.
        withcircles: boolean (optional)
            Set to True to include circles on the plot showing
            the overlapping ranges of the positioners.
        simple: boolean (optional)
            If True, show the arms as thick lines.
            If False, show the arms as realistic rectangles.
        trivial: boolean (optional)
            If True, show all positioners as a simple +.
            If False, plot all positioners.
        showplot: boolean, optional
            If True (the default) show the plot when finished.
            Set of False if you are overlaying several plots, and
            this plot is not the last.
            
        :Returns:
        
        plotfig: matplotlib Figure object
            A matplotlib figure containing the plot.
        
        """
        if plotting is None:
            logger.warning("Plotting is disabled")
            return
        
        # Determine a figure size appropriate to the grid size,
        # and some sensible line widths.
        maxrowcol = max(self.columns, self.rows)
        if maxrowcol > 10:
            figsize = (15,12)
        else:
            figsize = (10,8)
        
        title = "Grid of %d positioners in %d cols x %d rows with pitch=%.3f:" % \
            (self.positioner_count, self.columns, self.rows, self.pitch)
        if description:
            title += "\n%s" % description
        plotfig = plotting.new_figure(1, figsize=figsize, stitle=title)

        # Plot the positioners.
        for positioner in self.positioners:
            if positioner is not None:
                if trivial:
                    # For a trivial plot, simply mark the location of
                    # each positioner with a +.
                    plotting.plot_xy( positioner.x_centre_focal,
                                      positioner.y_centre_focal,
                                      plotfig=plotfig, linefmt='b+',
                                      linestyle=' ', showplot=False )
                else:
                    # Plot the actual configuration of each positioner.
                    plotfig = positioner.plot(plotfig=plotfig,
                                              showplot=False)
               
        # If there is a list of targets, show these on the plot as well,
        # with x symbols.
        if targetlist:
            xtargets = []
            ytargets = []
            for target in targetlist:
                xtargets.append(target[1])
                ytargets.append(target[2])
                plotting.plot_xy( xtargets, ytargets,
                                  plotfig=plotfig, linefmt='cx', linestyle=' ',
                                  showplot=False )

        if showplot:
            plotting.show_plot()
            plotting.close()
        return plotfig

    def plot_subset(self, ident_list, description=''):
        """
        
        Plot a subset the positioners in the grid.
        
        :Parameters:
        
        ident_list: list of int
            A list of all the positioners to be plotted.
        description: str (optional)
            Optional description to be added to the plot title.
       
        """
        if plotting is None:
            logger.warning("Plotting is disabled")
            return

        last_positioner = None
        for ident in ident_list:
            positioner = self.get_positioner(ident)
            if positioner is not None:
                last_positioner = positioner
                positioner.plot(showplot=False)
        if last_positioner is not None:
            last_positioner.plot(showplot=True, description=description)

    def count(self):
        """
        
        Count the number of positioners in the grid.
        
        :Returns:
        
        count: int
            The total number of positioners.
        
        """
        count = 0
        for positioner in self.positioners:
            if positioner is not None:
                count += 1
        return count

    def reset(self):
        """
        
        Reset all positioners in the grid to their initial state.
        
        """
        for positioner in self.positioners:
            if positioner is not None:
                positioner.initialise()

    def init_counters(self):
        """
        
        Initialise the conflict counters
        
        """
        self.counters = [0, 0, 0, 0, 0, 0, 0, 0]

    def counters_str(self):
        """
        
        Return a summary of the conflict counters in a string.
        
        """
        strg = "Summary of conflicts:\n"
        strg += "  %d positioners ok\n" % self.counters[Positioner.CONFLICT_OK]
        strg += "  %d positioners with unreachable targets\n" % self.counters[Positioner.CONFLICT_UNREACHABLE]
        strg += "  %d locked (broken) positioners\n" % self.counters[Positioner.CONFLICT_LOCKED]
        strg += "  %d targets too close\n" % self.counters[Positioner.CONFLICT_TOOCLOSE]
        strg += "  %d conflicts between metrology zones\n" % self.counters[Positioner.CONFLICT_METROLOGY]
        strg += "  %d conflicts between beta arms\n" % self.counters[Positioner.CONFLICT_BETA]
        strg += "  %d conflicts between datum actuators\n" % self.counters[Positioner.CONFLICT_DATUM]
        strg += "  %d conflicts between fibre and beta arm" % self.counters[Positioner.CONFLICT_FIBRE_SNAG]
        return strg

    #+++
    # Class variables defining the location of information within the
    # conflict list.
    CL_ID = 0
    CL_R_OK = 1
    CL_L_OK = 2
    CL_R_DIFF = 3
    CL_L_DIFF = 4
    # Class variables defining the location of information within the
    # reply list.
    R_IDENT = 0
    R_OK = 1
    R_PARITY = 2
    R_DIFFICULTY = 3
    R_ALT_DIFFICULTY = 4
    def check_targets(self, fibretargets, all_neighbours=False,
                      setup_positioners=False):
        """
           
        Check a collection of targets, with nominated positioners, to
        determine if any will result in conflict between positioners
        in the grid when in their final configuration. The best
        parity is suggested for each positioner.
          
        Target fibre coordinates are expressed in polar coordinates
        on the MOONS focal plane.
           
        :Parameters:
   
        fibretargets: list of (ident, rfibre, thfibre)
            List of target fibre positions required for each
            positioner. 
        all_neighbours: boolean (optional)
            If True each positioner is checked for a conflict with all
            its known neighbours, regardless of whether they are included
            in the fibretargets list.
            If False (the default) each positioner is checked for a
            conflict with other positioners in the fibretarget list.
            True is more efficient when checking large collections of
            targets, for example an entire grid configuration.
        setup_positioners: boolean (optional)
            If True, two extra passes through the data are used to setup
            the fibre positioners into the recommended configuration.
            Set to True if the final configuration of the fibre
            positioners is important, or False (the default) if only the
            reply is important and performance is the priority.
              
        :Returns:
         
        reply: list of (ident, ok, parity, difficulty, alt_difficulty)
            ident is the fibre positioner identifier
            ok is True if the target can be reached without conflict
            parity is the recommended parity
            difficulty is the degree of difficulty at the recommended parity
            alt_difficulty is the degree of dificulty at the other parity.
                    
        """
        assert isinstance(fibretargets, (tuple,list))
        assert isinstance(fibretargets[0], (tuple,list))
        assert len(fibretargets) > 1
        assert len(fibretargets[0]) == 3
 
        # First assign all the targets and record whether they can be reached
        # at each parity setting. Make a list of all the positioners included
        # in the list.
        conflict_table = []
        positioner_ids = []
        # For each target
        for (ident, rfibre, thfibre) in fibretargets:
            logger.debug("Checking target %d at r=%f, theta=%f" % \
                         (ident, rfibre, math.degrees(thfibre)))
            # Locate the positioner with the given identifier.
            positioner = self.get_positioner(ident)
            if positioner is not None:
                positioner_ids.append(ident)
                # Attempt to assign the target to the positioner at RIGHT parity
                # and record the result
                if positioner.set_target(rfibre, thfibre, parity=PARITY_RIGHT):
                    # Target can be reached
                    ok_right = True
                    diff_right = 0.0
                else:
                    # Target cannot be reached.
                    ok_right = False
                    diff_right = 1.0
                # Attempt to assign the target to the positioner at LEFT parity
                # and record the result
                if positioner.set_target(rfibre, thfibre, parity=PARITY_LEFT):
                    # Target can be reached
                    ok_left = True
                    diff_left = 0.0
                else:
                    # Target cannot be reached.
                    ok_left = False
                    diff_left = 1.0
                conflict_table.append( [ident, ok_right, ok_left,
                                        diff_right, diff_left] )
            else:
                # Positioner is not known
                strg = "Positioner %d is not known." % ident
                logger.error(strg)
                #raise ValueError(strg)
                conflict_table.append( [ident, False, False, 1.0, 1.0] )
        if logger.getEffectiveLevel() >= logging.DEBUG:
            logger.debug("After first pass: ID R_OK L_OK R_DIFF L_DIFF")
            for entry in conflict_table:
                logger.debug(str(entry))
                                     
        # Now make a second pass through the target list and check each
        # positioner for a potential conflict between neighbours.
        r_record = 0
        for (ident, rfibre, thfibre) in fibretargets:
            logger.debug("Second pass for positioner %d" % ident)
            assert( conflict_table[r_record][self.CL_ID] == ident)
            positioner = self.get_positioner(ident)
            # Only make the test if the positioner was able to reach its target.
            if conflict_table[r_record][self.CL_R_OK]:   
                # Check for conflict at the RIGHT parity
                positioner.set_target(rfibre, thfibre, parity=PARITY_RIGHT)
                if all_neighbours:
                    # Check all neighbours.
#                     logging.debug("CHECKING ALL NEIGHBOURS")
                    if positioner.in_conflict_with_neighbours():
                        # Change the status of this parity to not ok.
                        conflict_table[r_record][self.CL_R_OK] = False
                        conflict_table[r_record][self.CL_R_DIFF] = 1.0
                    else:
                        # Calculate the difficulty for this parity
                        right_diff = positioner.difficulty_value()
                        conflict_table[r_record][self.CL_R_DIFF] = right_diff
                else:
                    # Only check neighbours which are part of this subset.
#                     logging.debug("ONLY CHECKING THE GIVEN SUBSET")
                    subset = []
                    neighbours = positioner.get_neighbours()
#                     logging.debug("Number of neighbours=%d" % len(neighbours))
                    for neighbour in neighbours:
#                         logging.debug("Checking neighbour %s" % neighbour.name)
                        if neighbour.ident in positioner_ids:
#                             logging.debug("Inside subset")
                            subset.append(neighbour)
                            if positioner.in_conflict_with(neighbour):
                                # Change the status of this parity to not ok.
                                conflict_table[r_record][self.CL_R_OK] = False
                                conflict_table[r_record][self.CL_R_DIFF] = 1.0
                                logger.debug("RIGHT: Positioner %s in conflict with %s" % \
                                      (positioner.name, neighbour.name))
#                                 if GRAPHICAL_DEBUGGING:
#                                     title = "%s in conflict with %s" % \
#                                         (positioner.name, neighbour.name)
#                                     title += "\n" + positioner.conflict_reason
#                                     plotfig = positioner.plot(showplot=False,
#                                                 description=title)
#                                     neighbour.plot(plotfig=plotfig,
#                                                 showplot=True)
                    if conflict_table[r_record][self.CL_R_OK]:
                        # Calculate the difficulty for this parity
                        right_diff = positioner.difficulty_value(others=subset)
                        logger.debug("RIGHT: Positioner %s ok. Difficulty %f" % \
                                (positioner.name, right_diff))
                        conflict_table[r_record][self.CL_R_DIFF] = right_diff
            # Only make the test if the positioner was able to reach its target.
            if conflict_table[r_record][self.CL_L_OK]:                            
                # Check for conflict at the LEFT parity
                positioner.of(rfibre, thfibre, parity=PARITY_LEFT)
                if all_neighbours:
                    # Check all neighbours.
#                     logging.debug("CHECKING ALL NEIGHBOURS")
                    if positioner.in_conflict_with_neighbours():
                        # Change the status of this parity to not ok.
                        conflict_table[r_record][self.CL_L_OK] = False
                        conflict_table[r_record][self.CL_L_DIFF] = 1.0
                    else:
                        # Calculate the difficulty for this parity
                        left_diff = positioner.difficulty_value()
                        conflict_table[r_record][self.CL_L_DIFF] = left_diff
                else:
                    # Only check neighbours which are part of this subset.
#                     logging.debug("ONLY CHECKING THE GIVEN SUBSET")
                    neighbours = positioner.get_neighbours()
                    subset = []
                    for neighbour in neighbours:
#                         logging.debug("Checking neighbour %s" % neighbour.name)
                        if neighbour.ident in positioner_ids:
#                             logging.debug("Inside subset")
                            subset.append(neighbour)
                            if positioner.in_conflict_with(neighbour):
                                # Change the status of this parity to not ok.
                                conflict_table[r_record][self.CL_L_OK] = False
                                conflict_table[r_record][self.CL_L_DIFF] = 1.0
                                logger.debug("LEFT: Positioner %s in conflict with %s" % \
                                      (positioner.name, neighbour.name))
#                                 if GRAPHICAL_DEBUGGING:
#                                     title = "%s in conflict with %s" % \
#                                         (positioner.name, neighbour.name)
#                                     title += "\n" + positioner.conflict_reason
#                                     plotfig = positioner.plot(showplot=False,
#                                                 description=title)
#                                     neighbour.plot(plotfig=plotfig,
#                                                 showplot=True)
                    if conflict_table[r_record][self.CL_L_OK]:
                        # Calculate the difficulty for this parity
                        left_diff = positioner.difficulty_value(others=subset)
                        logger.debug("LEFT: Positioner %s ok. Difficulty %f" % \
                                (positioner.name, left_diff))
                        conflict_table[r_record][self.CL_L_DIFF] = left_diff
            r_record = r_record + 1
        if logger.getEffectiveLevel() == logging.DEBUG:
            logger.debug("After second pass: ID R_OK L_OK R_DIFF L_DIFF")
            for entry in conflict_table:
                logger.debug(str(entry))

        # Now make a third pass through the target list, choose the best
        # configuration, and construct the reply.
        reply = []
        r_record = 0
        for (ident, rfibre, thfibre) in fibretargets:
            logger.debug("Third pass for positioner %d" % ident)
            assert( conflict_table[r_record][0] == ident)
            positioner = self.get_positioner(ident)
            if conflict_table[r_record][self.CL_R_OK] and \
               conflict_table[r_record][self.CL_L_OK]:
                # No conflict at either parity. Choose the one with the
                # smallest difficulty.
                if conflict_table[r_record][self.CL_R_DIFF] < \
                   conflict_table[r_record][self.CL_L_DIFF]:
                    # Best choice RIGHT, but LEFT also possible.
                    ok = True
                    parity_choice = PARITY_RIGHT
                    difficulty = conflict_table[r_record][self.CL_R_DIFF]
                    alt_difficulty = conflict_table[r_record][self.CL_L_DIFF]
                else:
                    # Best choice LEFT, but RIGHT also possible.
                    ok = True
                    parity_choice = PARITY_LEFT
                    difficulty = conflict_table[r_record][self.CL_L_DIFF]
                    alt_difficulty = conflict_table[r_record][self.CL_R_DIFF]
            elif conflict_table[r_record][self.CL_R_OK]:
                # No conflict at RIGHT parity only
                ok = True
                parity_choice = PARITY_RIGHT
                difficulty = conflict_table[r_record][self.CL_R_DIFF]
                alt_difficulty = 1.0
            elif conflict_table[r_record][self.CL_L_OK]:
                # No conflict at LEFT parity only
                ok = True
                parity_choice = PARITY_RIGHT
                difficulty = conflict_table[r_record][self.CL_L_DIFF]
                alt_difficulty = 1.0
            else:
                # Conflict
                ok = False
                parity_choice = PARITY_RIGHT
                difficulty = 1.0
                alt_difficulty = 1.0
            reply.append( [ident, ok, parity_choice, difficulty, alt_difficulty] )
            r_record = r_record + 1

        if logger.getEffectiveLevel() == logging.DEBUG:
            logger.debug("After third pass: ID OK PARITY DIFF ALT_DIFF")
            for entry in reply:
                logger.debug(str(entry))

        # Now make a fourth and fifth pass through the target list, in the
        # fourth pass moving the positioners to their final configurations
        # and in the fifth pass checking for any unresolved conflict.
        # These last two passes are mainly for simulation purposes, since
        # the reply has already been constructed.
        if setup_positioners:
            r_record = 0
            for (ident, rfibre, thfibre) in fibretargets:
                assert( conflict_table[r_record][self.R_IDENT] == ident)
                positioner = self.get_positioner(ident)
                positioner.initialise()
                parity_choice = reply[r_record][self.R_PARITY]
                logger.debug("Fourth pass for positioner %d at parity %s" % \
                      (ident, util.elbow_parity_str(parity_choice)))
                positioner.in_conflict = False
                positioner.set_target(rfibre, thfibre, parity=parity_choice)
                r_record = r_record + 1
            r_record = 0
            for (ident, rfibre, thfibre) in fibretargets:
                assert( conflict_table[r_record][self.R_IDENT] == ident)
                positioner = self.get_positioner(ident)
                parity_choice = reply[r_record][self.R_PARITY]
                logger.debug("Fifth pass for positioner %d at parity %s" % \
                      (ident, util.elbow_parity_str(parity_choice)))
                # Go to the chosen configuration and log any unresolved conflict.
                if positioner.set_target(rfibre, thfibre, parity=parity_choice):
                    if all_neighbours:
                        if positioner.in_conflict_with_neighbours():
                            logger.warning(positioner.conflict_reason)
                    else:
                        # Only check neighbours which are part of this subset.
                        neighbours = positioner.get_neighbours()
                        for neighbour in neighbours:
                            if neighbour.ident in positioner_ids:
                                if positioner.in_conflict_with(neighbour):
                                    logger.debug("Positioner %s in conflict with %s" % \
                                            (positioner.name, neighbour.name))
                                    logger.warning(positioner.conflict_reason)
                                    if GRAPHICAL_DEBUGGING:
                                        title = "%s in conflict with %s" % \
                                            (positioner.name, neighbour.name)
                                        title += "\n" + positioner.conflict_reason
                                        plotfig = positioner.plot(showplot=False,
                                                    description=title)
                                        neighbour.plot(plotfig=plotfig, showplot=True)
                else:
                    logger.warning(positioner.conflict_reason)
                self.counters[positioner.conflict_type] += 1
                r_record = r_record + 1
        return reply

    def test_pair(self, target1, parity1, target2, parity2):
        """
           
        Test for a conflict between two targets assigned to specific
        positioners for a particular combination of parities.
          
        Target fibre coordinates are expressed in polar coordinates
        on the MOONS focal plane.
           
        :Parameters:

        target1: tuple of (ident, rfibre, thfibre)
            The positioner identifier and proposed target R and theta
            coordinates for the first target.
        parity1: int
            The proposed elbow parity for the first target.
        target2: tuple of (ident, rfibre, thfibre)
            The positioner identifier and proposed target R and theta
            coordinates for the second target.
        parity2: int
            The proposed elbow parity for the second target.
              
        :Returns:
         
        (in_conflict, difficulty1, difficulty2)
            in_conflict is True if the two positioners are in conflict,
            or if any of the targets cannot be reached, and False if both
            targets are ok and the positioners are not in conflict.
            difficulty1 and difficulty2 give the degree of difficulty
            estimates for each positioner in the proposed configuration
            (in the range 0.0-1.0).
                    
        """
        assert isinstance(target1, (tuple,list))
        assert isinstance(target2, (tuple,list))
        assert len(target1) >= 3
        assert len(target2) >= 3
        # Unpack the target parameters
        (ident1, rfibre1, thfibre1) = target1
        (ident2, rfibre2, thfibre2) = target2
        
        # Verify that both positioners exist.
        positioner1 = self.positioners[ident1]
        if positioner1 is None:
            # Positioner is not known
            strg = "Positioner %d is not known." % ident1
            logger.error(strg)
            return (True, 1.0, 1.0)
        positioner2 = self.positioners[ident2]
        if positioner2 is None:
            # Positioner is not known
            strg = "Positioner %d is not known." % ident2
            logger.error(strg)
            return (True, 1.0, 1.0)
        
        # Attempt to assign the given targets to the positioners.
        if not positioner1.set_target(rfibre1, thfibre1, parity=parity1):
            # Target cannot be reached.
            logger.warning(positioner1.conflict_reason)
            return (True, 1.0, 1.0)
        if not positioner2.set_target(rfibre2, thfibre2, parity=parity2):
            # Target cannot be reached.
            logger.warning(positioner2.conflict_reason)
            return (True, 1.0, 1.0)

        # Test for conflict.
        in_conflict = positioner1.in_conflict_with( positioner2 )
        if in_conflict:
            difficulty1 = 1.0
            difficulty2 = 1.0
        else:
            difficulty1 = positioner1.difficulty_value(others=[positioner2])
            difficulty2 = positioner2.difficulty_value(others=[positioner1])
        return (in_conflict, difficulty1, difficulty2)

    CONFLICT_RR = 0
    CONFLICT_RL = 1
    CONFLICT_LR = 2
    CONFLICT_LL = 3
    def check_pair(self, target1, target2):
        """
           
        Check for a conflict between two targets assigned to specific
        positioners, trying all four combinations of parities (RR, RL,
        LR, LL) and ranking the ones that don't conflict.
          
        Target fibre coordinates are expressed in polar coordinates
        on the MOONS focal plane.
           
        :Paramete,
        rs:

        target1: tuple of (ident, rfibre, thfibre)
            The positioner identifier and proposed target R and theta
            coordinates for the first target.
        target2: tuple of (ident, rfibre, thfibre)
            The positioner identifier and proposed target R and theta
            coordinates for the second target.
              
        :Returns:
         
        (rr, rl, lr, ll)
            A conflict list giving the result of the test at each
            parity combination (e.g. rl means target1 at RIGHT parity
            and target2 at LEFT parity). Each entry in the list
            contains a difficulty value, in the range 0.0 to 1.0,
            giving the likelihood of conflict. Values greater than or 
            equal to 1.0 indicate conflict.
                    
        """
        # Try each combination of target parities.
        conflict_list = []
        parity_combinations = [(PARITY_RIGHT, PARITY_RIGHT),
                               (PARITY_RIGHT, PARITY_LEFT),
                               (PARITY_LEFT, PARITY_RIGHT),
                               (PARITY_LEFT, PARITY_LEFT),
                               ]
        for (parity1, parity2) in parity_combinations:
            (in_conflict, diff1, diff2) = self.test_pair(target1, parity1,
                                                         target2, parity2)
            if in_conflict:
                conflict_list.append( 1.0 )
            else:
                combined_diff = (diff1 + diff2)/2.0
                conflict_list.append( combined_diff )
                
        return conflict_list


def count_conflicts( conflict_reply ):
    """
    
    Count the number of conflicts in a conflict check reply list.
    
    :Parameters:
    
    conflict_reply: list of (ident, ok, parity, difficulty, alt_difficulty)
        The reply list returned by the PositionerGrid  "check_targets"
        method.
        
    :Returns:
    
    (ngood, nconflicts): list of int
        The number of good replies and the number of conflicts.
    
    """
    ngood = 0
    nconflicts = 0
    for entry in conflict_reply:
        if entry[PositionerGrid.R_OK]:
            ngood += 1
        else:
            nconflicts += 1
    return (ngood, nconflicts)
        

if __name__ == "__main__":
    print("\nMOONS OPS/FPS shared library test...")

    PLOTTING = True
    if plotting is None:
        PLOTTING = False
        
    # Decide which tests to run.
    POSITIONER_TESTS = True
    MOTOR_TESTS = True
    POSITIONER_GRID_TESTS = True
    POSITIONER_PAIR_TESTS = True

    if POSITIONER_TESTS:
        def set_target_angles( positioner, r, theta, parity, description):
            """
             
            Test function to set a positioner to a given target, display
            and then plot the configuration, including the given description
            in the title.
             
            """
            if r is not None and theta is not None and parity is not None:
                if not positioner.set_target(r, math.radians(theta), parity):
                    print("***Failed to set target (%f,%f)" % (r,theta))
            print("")
            print(positioner)
            strg = "Arm angles are alpha=%.3f, beta=%.3f" % positioner.get_arm_angles()
            strg += " (orient=%.3f)" % positioner.orient
            print(description + " " + strg)
            if PLOTTING:
                positioner.plot(showplot=True,
                                description=description + "\n" + strg)

        def set_motor_angles( positioner, alpha, beta, description):
            """
             
            Test function to set a positioner's arms to given orientations,
            display and then plot the configuration, including the given 
            description in the title.
             
            """
            if alpha is not None and beta is not None:
                positioner.set_arm_angles(alpha, beta)
            print("")
            print(positioner)
            strg = " Arm angles are alpha=%.3f, beta=%.3f" % positioner.get_arm_angles()
            strg += " (orient=%.3f)" % math.degrees(positioner.orient)
            print(description + " " + strg)
            if PLOTTING:
                positioner.plot(showplot=True,
                                description=description + "\n" + strg)

        # Create a positioner and move its arms to different orientations.
        if MOTOR_TESTS:
            print("\n----------------------------------------------------------------------")
            print("Positioner motor tests...")
            import numpy as np
            positioner = Positioner(1, 0.0, 0.0, math.radians(0.0), 17, 19)
            beta = BETA_DEFAULT
            for alpha in np.arange(ALPHA_LIMITS[0], ALPHA_LIMITS[1], 60.0):
                set_motor_angles( positioner, alpha, beta, "Moving alpha.")
            set_motor_angles( positioner, ALPHA_LIMITS[1], beta, "Moving alpha.")
            positioner.initialise()
            alpha = ALPHA_DEFAULT
            for beta in np.arange(BETA_LIMITS[0], BETA_LIMITS[1], 60.0):
                set_motor_angles( positioner, alpha, beta, "Moving beta.")
            set_motor_angles( positioner, alpha, BETA_LIMITS[1], "Moving beta.")
            del positioner
            for orient in [0.0, 45.0, 225.0]:
                positioner = Positioner(1, 0.0, 0.0, math.radians(orient), 17, 19)
                set_motor_angles( positioner, 42.0, -120.0, "Changing orientation.")
                del positioner

        # Create a positioner, assign a target to it and then plot it,
        # showing the avoidance zones used by the conflict detection
        # function.
        # ident, r_centre_focal, theta_centre_focal, orient, column, row
        print("\n----------------------------------------------------------------------")
        print("Positioner tests...")
        positioner = Positioner(1, 0.0, 0.0, math.radians(0.0), 17, 19)
        strg = "Inner value=%f, " % positioner.inner
        strg += "outer value=%f " % positioner.outer
        strg += "and critical length=%f" % positioner.lengthcrit
        print(strg)
        set_target_angles( positioner, None, None, None, "Starting position.")
        set_target_angles( positioner, positioner.outer, 90.0, PARITY_RIGHT,
                           "Fully stretched position.")
        set_target_angles( positioner, positioner.lengthcrit, 90.0, PARITY_LEFT,
                           "Opposite fully stretched position.")
        set_target_angles( positioner, 20.0, 15.0, PARITY_RIGHT,
                           "r=20.0 theta=15.0 with right-armed parity.")
        set_target_angles( positioner, 20.0, 15.0, PARITY_LEFT,
                           "r=20.0 theta=15.0 with left-armed parity.")
        set_target_angles( positioner, 12.0, 180.0, PARITY_RIGHT,
                           "r=12.0 theta=180.0 with right-armed parity.")
        set_target_angles( positioner, 12.0, 180.0, PARITY_LEFT,
                           "r=12.0 theta=180.0 with left-armed parity.")
        set_target_angles( positioner, 16.0, -65.0, PARITY_RIGHT,
                           "r=16.0 theta=-65.0 with right-armed parity.")
        set_target_angles( positioner, 16.0, -65.0, PARITY_LEFT,
                           "r=16.0 theta=-65.0 with left-armed parity.")
         
        print("\nAdd a second positioner and check for conflict")
        positioner2 = Positioner(2, 20, math.radians(45.0), 0.0, 16, 21)
        set_target_angles(positioner2, 30.0, 30.0, PARITY_RIGHT,
                        "r=30.0 theta=30.0 with right-armed parity");
     
        # Test for conflict
        print("\nin_conflict_with function test.")
        if (positioner2.in_conflict_with( positioner ) ):
            print("Positioners are in conflict (" + \
                  positioner.conflict_reason + ")")
        else:
            print("Positioners are NOT in conflict")
     
        print("\ncheck_for_conflict function test.")
        in_conflict = check_for_conflict(0.0, 0.0, 0.0,
                            16.0, math.radians(-65.0), PARITY_LEFT,
                            20.0, math.radians(45.0), 0.0,
                            30.0, math.radians(30.0), PARITY_RIGHT)
        if ( in_conflict ):
            print("Positioners are in conflict (" + positioner.conflict_reason + ")")
        else:
            print("Positioners are NOT in conflict")
    
    if POSITIONER_GRID_TESTS:
        print("\n----------------------------------------------------------------------")
        print("PositionerGrid test...\n")
        # Define a fibre positioner grid in Cartesian coordinates
        positioner_configs_xy = [
                # id,     r,       theta,       orient,      col, row
                ( 1,   0.00000,   0.00000, math.radians(0.0), 17, 19),
                ( 2,  12.62935,  21.87468, math.radians(0.0), 17, 20),
                ( 3,  25.25870,   0.00000, math.radians(0.0), 18, 19),
                ( 4,  12.62935, -21.87470, math.radians(0.0), 17, 18),
                ( 5, -12.62940, -21.87470, math.radians(0.0), 16, 18),
                ( 6, -25.25870,   0.00000, math.radians(0.0), 16, 19),
                ( 7, -12.62940,  21.87468, math.radians(0.0), 16, 20),
                ( 8,   0.00000,  43.74782, math.radians(0.0), 17, 21),
                ( 9,  25.25737,  43.74705, math.radians(0.0), 18, 21),
                (10,  37.88672,  21.87391, math.radians(0.0), 18, 20),
                (11,  50.51474,   0.00000, math.radians(0.0), 19, 19),
                (12,  37.88672, -21.87390, math.radians(0.0), 18, 18),
                (13,  25.25737, -43.74700, math.radians(0.0), 18, 17),
                (14,   0.00000, -43.74780, math.radians(0.0), 17, 17),
                (15, -25.25740, -43.74700, math.radians(0.0), 16, 17),
                (16, -37.88670, -21.87390, math.radians(0.0), 15, 18),
                (17, -50.51470,   0.00000, math.radians(0.0), 15, 19),
                (18, -37.88670,  21.87391, math.radians(0.0), 15, 20),
                (19, -25.25740,  43.74705, math.radians(0.0), 16, 21)
            ]
        # Convert to a grid configuration in polar coordinates.
        positioner_configs = []
        for config in positioner_configs_xy:
            newconfig = []
            newconfig.append(config[0])
            (r,theta) = util.cartesian_to_polar(config[1], config[2])
            newconfig.append(r)
            newconfig.append(theta)
            newconfig.append(config[3])
            newconfig.append(config[4])
            newconfig.append(config[5])
            positioner_configs.append(newconfig)
        # Create a fibre positioner grid object.
        positioner_grid = PositionerGrid( positioner_configs )
        print( positioner_grid )
         
        print("\nConflict check test 1 (check_targets with 19 positioners)...\n")
        # Define a set of targets in Cartesian coordinates (for convenience),
        # then convert to polar coordinates to match the agreed interface.
        fibrecoords = [(1, 0.87, 10),
                       (2, 26.50, 42),
                       (3, 22.94, -16),
                       (4, 8.80, -33),
                       (5, 6.86, -32),
                       (6, -11.57, 14.0),
                       (7, -16.00,  0.0),
                       (8, -20.0, 60.0),
                       (9,  16.2, 43.75),
                       (10, 48.0, 19.0),
                       (11, 45.0, 16.0),
                       (13, 10.0, -38.0),
                       (14, -24.9, -43.75),
                       (15, -25.0, -60.0 ),
                       (16, -50.7, -10.0),
                       (17, -54.0, -18.0),
                       (18, -37.0, 40.0),
                       (19, -20.0, 35.0)
                       ]
        fibretargets = []
        for (ident, x, y) in fibrecoords:
            (r, theta) = util.cartesian_to_polar(x, y)
            fibretargets.append( (ident, r, theta) )
     
        time.sleep(1) # Ensure any logging text appears before the reply string.
        strg = "Target suggestions [ident, r, theta] (focal plane coordinates):"
        for entry in fibretargets:
            theta_degrees = math.degrees(entry[2])
            strg += "\n  [%d, %.5f, %.5f (%.5f deg)]" % \
                (entry[0], entry[1], entry[2], theta_degrees)
        print(strg)
         
        # Now assign these targets to the fibre positioner grid and
        # check for conflicts.
        replies = positioner_grid.check_targets(fibretargets,
                                                setup_positioners=True)
        print(positioner_grid.counters_str())
         
        time.sleep(1) # Ensure any logging text appears before the reply string.
        strg = "Reply string 1 [ident, ok, parity, difficulty, alt_difficulty]: "
        for entry in replies:
            strg += "\n  [%d, %s, %s, %.4f, %.4f]" % \
                (entry[0], str(entry[1]), util.elbow_parity_str(entry[2]),
                 entry[3], entry[4])
        print(strg)
         
        targetlist = []
        for ident, rtarget, thtarget in fibretargets:
            xtarget, ytarget = util.polar_to_cartesian(rtarget, thtarget)
            targetlist.append( [ident, xtarget, ytarget] )
        if PLOTTING:
            plotfig = positioner_grid.plot(targetlist=targetlist,
                            description="(Green=no target; " + \
                            "Magenta=target unreachable; Black=target assigned; " + \
                            "Red=In conflict; x=targets)")

        print("\nConflict check test 2 (check_targets with just 2 positioners)...\n")
        # Define two simple targets for 16 and 17
        (r1, theta1) = util.cartesian_to_polar(-50.7, -10.0)
        (r2, theta2) = util.cartesian_to_polar(-54.0, -18.0)
        fibretargets = [(16, r1, theta1),
                        (17, r2, theta2)]
     
        # Now assign these targets to the fibre positioner grid and
        # check for conflicts.
        replies = positioner_grid.check_targets(fibretargets)
          
        time.sleep(1) # Ensure any logging text appears before the reply string.
        strg = "Reply string 2 [ident, ok, parity, difficulty, alt_difficulty]: "
        for entry in replies:
            strg += "\n  [%d, %s, %s, %.4f, %.4f]" % \
                (entry[0], str(entry[1]), util.elbow_parity_str(entry[2]),
                 entry[3], entry[4])
        print(strg)
     
    if POSITIONER_PAIR_TESTS:
        print("\nConflict check test 3 (test_pair)...\n")
        # Try the same pair of targets in 4 different parity configurations.
        (r1, theta1) = util.cartesian_to_polar(-50.7, -10.0)
        (r2, theta2) = util.cartesian_to_polar(-54.0, -18.0)
        target1 = (16, r1, theta1)
        target2 = (17, r2, theta2)
        parity_checks = ([PARITY_RIGHT, PARITY_RIGHT],
                         [PARITY_RIGHT, PARITY_LEFT],
                         [PARITY_LEFT,  PARITY_RIGHT],
                         [PARITY_LEFT,  PARITY_LEFT]
                         )
        for (parity1, parity2) in parity_checks:
            # Construct a string from the first letter of each parity name.
            pstrg1 = util.elbow_parity_str(parity1)
            pstrg2 = util.elbow_parity_str(parity2)
            paritystrg = "%.1s%.1s" % (pstrg1, pstrg2)
            (conflict, diff1, diff2) = positioner_grid.test_pair( target1, parity1,
                                                                  target2, parity2 )
            if conflict:
                strg = "Positioners 16 and 17 are in conflict (%s)." % paritystrg
            else:
                strg = "Positioners 16 and 17 are not in conflict (%s)." % paritystrg
                strg += " Difficulty estimates: %.3f, %.3f." % (diff1, diff2)
            print( strg )
            if PLOTTING:
                positioner_grid.plot_subset( [16, 17], description=strg )
     
        print("\nConflict check test 4 (check_pair)...\n")
        conflict_list = positioner_grid.check_pair(target1, target2)
        print("Conflict list: rr=%.3f, rl=%.3f, lr=%.3f, ll=%.3f" % tuple(conflict_list))

    print("Test finished.")
