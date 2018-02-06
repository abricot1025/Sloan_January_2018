try:
    import mocpath.conflict.fps_shared as fps
except ImportError:
    import fps_shared as fps

import numpy as np
from numpy import linalg as LA
import math

import random

import json as simplejson
import matplotlib.pyplot as plt

import parameters as param

__author__ = "Tao Dominique, Laleh Makarem"
__copyright__ = "Belong to EPFL and MOONS"
__credits__ = ["Rob Knight", "Peter Maxwell", "Gavin Huttley",
                    "Matthew Wakefield"]
__version__ = "Developed under python 3.4.3, but should work with python 2.7.6 too"
__email__ = "dominique.tao@epfl.ch"
__status__ = "in Progress (I guess)"

"""
MOONS: Potential field and finite state machine(decision layer)

A prototype module written in python containing function used in MOONS
for moving the positioners (collision avoidance with potential field
decentralized navigation function) and coordinating them with priorities
(priority coordination).

It is expected that the final version will be translated into C or C++
for greater efficiency and compatibility with ESO/VLT software.

The original version with the potential field (decentralized navigation function) was written
by Laleh Makarem. Then the extension of the decision layer (finite state machine)
was then added on top of it by Dominique Tao.

The list of modification concerns the implementation of
finite state machine (decision layer) implemented by Dominique Tao:

27.02.2017: Read priority from OPS file
            Repulsive and attractive force depending on priority
03.03.2017: Variable zone of influence of the force
17.03.2017: Mathematical modification of the DNF
            Goal extension of DNF of MOONS of Laleh work for MOONS
            with beautiful mathematical function
            but did not succeed
24.03.2017: Added first attribute of finite state machine:
            ON  & OFF (deactivate the attractive force)
28.03.2017: Added second and third attribute of finite state machine:
            ID and PR (remember the positioner ID and priority (or
            priority of the guy who put this neigh pos off if it was OFF)
            that put it OFF)
03.04.2017: Change of parameters in force, remove variability (removing the
            adaptive zone of influence of the force )and thus
            reduction of complexity
18.04.2017: Finish correction of small bug with same priority (use
            random priority assignement and study curve of success depending on
            how priority is assigned to have a model (but not working...)),
            Correction of projection mistake,
            weird behaviors of positioners
            due to attributes of finite state machine,
            creation and test on new constructed test cases
27.04.2017: Correction of chain reaction,
            Creation of the fourth attribute of the finite state machine:
             INFL (is OFF when positioner reached its target !,otherwise ON,
             when OFF, this pos not able to influence state of other pos within
             its vicinity)
03.05.2017: Correction of chaine reaction of positioner with different
            priority influencing each others (somehow correction of deadlock)
09.05.2017: Change of tactic for same priority positioner creating of deadlock:
            with random priority assignement, introduction of a distance
            based heuristic to have more determine result,
            added noise to try solving deadlock with same priority
            Began presentation with reveal.js of my work (cool stuff with video)
xx.05.2017 : presentation
08.06.2017: creation of fourth variable:
            LAYER (used to manage, isolate two positioners with same priority
            that is creating a deadlock),
            creation of two priority: OPS_priority used in normal condition
                                      pseudo_priority used when fourth attribute > 0,
            Constant margin of status zone
16.06.2017: Remove random assignement of priority for deadlock, used a vector-distance
            strategy (fully deterministic strategy
             improvement but not really working well (not complete))
21.06.2017: change of stopping criteria for when solving deadlock of pos with same priorities
07.07.2017: change strategy to solve deadlock with same priority:
            creation of phantom positioner covering larger space ahead
18.07.2017:  Change strategy to solve deadlock with same priority:
            Previous one created uncessary movement and deadlocks...
            use a vector-distance based heuristic that identify if target and
            neigh positioner on same side or not in order to know if possibe deadlock
            (kind of point of view strategy)
31.07.2017: Although kind of working, still some deadlock...
            Change strategy to solve deadlock with same priority:
            Same as before but with different point of view !!! , since our problem is in 2D
            need several point of view !
03.08.2017: change of stopping criteria
08.08.2017: correction of bug related to collision and deadlock problem with strategy
12.08.2017: Final version of the code for the paper... (not time left to improve it...)

@author: Laleh Makarem
@author: Dominique Tao

"""


"""
####### TIPS TIPS TIPS TIPS TIPS TIPS TIPS TIPS#######

# With change in the velociy of the attraction force: CONST_ATTRACTIVE_FORCE 2 -> 1.8
            M1attractiveF  = -self.k11 * self.ggoal1 * ((self.agent.pseudo_priority+CONST_ATTRACTIVE_FORCE)/2)
            M2attractiveF  = -self.k21 * self.ggoal2 * ((self.agent.pseudo_priority+CONST_ATTRACTIVE_FORCE)/2)

# because otherwise the positioners was too fast and high risk of collision and bad behavior
# ==> BUT STILL NEED the attraction force to be big enough so the positioner can bypass the positioner that creates deadlock
# or reach target that are otherwise difficult to reach beause they are near other positioners that create a repulsive force!

## When changed what is above, need to change also
   REPULSIVE_TARGET_FACTOR :
#    used for when deadlock situation, the target of the other positioner creating the deadlock
#   create a isotropic force to this positioner!
#  if velocity of the posistioners above is smaller, the force created by the target need to be bigger, resulting in a bigger
# velocity of the positioner affected by this force ==> in order to bypass the other positioner
# that creates the deadlock or its repulsive force


"""

OFF=0
ON=1

#################################################################################################################
 ################################################
### THESE ARE THE PARAMETERS THAT I USED TO MODIFY###
 ################################################

# zone where the positioner is affected by the repulsive force of other positioners
STATUS_OFF_ZONE_INFLUENCE = param.STATUS_OFF_ZONE_INFLUENCE # 3.5 # 3.5
STATUS_ON_ZONE_INFLUENCE = param.STATUS_ON_ZONE_INFLUENCE # 9.0 # 8.8

LIMITE_FOR_NOISE = param.LIMITE_FOR_NOISE # was measured by running several time on several test cases

CSTE_FORCE = param.CSTE_FORCE #6.2# 5.8 # 5.9

# zone where the positioner can influence the status of the other positioner!
STATUSZONE = param.STATUSZONE  #6 # 6 # 6.5

# adding this cst_margin: avoid the unnecessary oscillation due to the too strict condition "if distance < statusZone"
# so less oscillation in the "ON/OFF statusForceAttr" thus takes less time to reach target TODO ?
CSTE_MARGIN_STATUSZONE = param.CSTE_MARGIN_STATUSZONE # for the status change zone
CSTE_MARGIN_STATUSZONE_STOPPING_DEADLOCK = param.CSTE_MARGIN_STATUSZONE_STOPPING_DEADLOCK # for the stopping criteria of whe not deadlock anymore
ARRIVED_PRIORITY = param.ARRIVED_PRIORITY

# the smaller it gets, the stronger the repulsice force from the repulsive target it becomes
REPULSIVE_TARGET_FACTOR = param.REPULSIVE_TARGET_FACTOR

# a constant added to the attrative force to better modulate it
CONST_ATTRACTIVE_FORCE = param.CONST_ATTRACTIVE_FORCE

#################################################################################################################

class Motion_Planning(object):
    """
    Class representing the algorithm of the motion planning (decentralize potential field) that makes
    the positioners move to their target
    """

    _params_defined = False
    @classmethod
    def define_common_params_Motion_Planning(cls):
        """

        Class method which defines the parameters which are common to all
        fibre positioners. These parameters only need to be calculated once.

        """

        cls.l1= fps.INS_POS_LENGTH1
        cls.l2w= fps.INS_POS_WIDTH2

        # Parameters of the motor velocity
        cls.k11 = 0.2
        cls.k12 = 0.00001 # 0.001
        cls.k13 = 0.2

        cls.k21 = 0.2
        cls.k22 = 0.00001 # 0.001
        cls.k23 = 0.2

        cls._params_defined=True


    def __init__(self, positioner):

        self.agent = positioner
        self.m = 0 # simulation time (value assigned at pas_config_dnf)

        if not Motion_Planning._params_defined:
            Motion_Planning.define_common_params_Motion_Planning()

        self.d = 1           # weighted distance is already normalized
        self.d2 = self.l2w/2 + 1          # half-width of the beta arm + 1 mm of margin

        # initialized the three terms of the potential function
        self.ggoal1 = 0
        self.ggoal2 = 0

        self.gcol1 = 0
        self.gcol2 = 0

        self.gcol_end1 = 0   # added for protection of fiber
        self.gcol_end2 = 0

        self.motor1VelList=[]
        self.motor1Addition=0

        self.motor2VelList=[]
        self.motor2Addition=0

    def weigthed_distance (self,q, q_center, rotation,l1,l2):
        """

        :param q, q_center:
        :param rotation:
        :param l1, l2:
        :return:
        """

        a = l1/2
        b = l2/2
        distance = math.sqrt((((q[0]-q_center[0])*math.cos(rotation) + (q[1]-q_center[1])*math.sin(rotation))**2)/(a**2) +
                             (((q[1]-q_center[1])*math.cos(rotation) - (q[0]-q_center[0])*math.sin(rotation))**2)/(b**2))
        return distance

    def get_angle_length_data(self):
        """
        Get the position data in variables with shorter names...
        cartesian position of the end-effector of current positioner in focal plane

        :return:
         r1: float
            alpha angle in radian
         r2: float
            beta angle in radian

         l1: float
            length of first arm
         l2: float
            length of the second arm

         r_beta:  float
            half a turn back to compensate for beta rotating from tacked in position

        """
        r1 = self.agent.motor1.position_array[self.m]
        r2 = self.agent.motor2.position_array[self.m]

        l1 = self.agent.length1
        l2 = self.agent.length2

        r_beta = r1 + r2 + math.pi # half a turn back to compensate for beta rotating from tacked in position

        return r1,r2,l1,l2,r_beta


    # add the attractive force to the target
    def add_attractive_force(self):
        """
        Calculate the coefficient for the attractive force of the positioner to its target
        :return:
        """
        if self.agent.motor1.position_array[-1] < 10000 and\
                math.fabs(self.agent.motor1.position_array[self.m] - self.agent.motor1.position_array[-1]) > 0.019 :
            self.ggoal1 += ((self.agent.motor1.position_array[self.m] - self.agent.motor1.position_array[-1])\
                / LA.norm(self.agent.motor1.position_array[self.m] - self.agent.motor1.position_array[-1]))

        if self.agent.motor2.position_array[-1] < 10000 and\
                math.fabs(self.agent.motor2.position_array[self.m] - self.agent.motor2.position_array[-1]) > 0.025:
            self.ggoal2 += ((self.agent.motor2.position_array[self.m] - self.agent.motor2.position_array[-1])\
                / LA.norm(self.agent.motor2.position_array[self.m] - self.agent.motor2.position_array[-1]))

    # Saturate on maximum
    def saturate_max(self,value, max_value):
        """
        :param value: float
        :param max_value: float
        :return:
        """

        re_value = value
        if value > max_value:
            re_value = max_value
        elif value < -max_value:
            re_value = - max_value

        return re_value

    def saturate_min(self,value,min_value):
        """
        :param value: float
        :param min_value: float
        :return:
        """
        re_value = value
        if 0 <= value < min_value/2:
            re_value = 0
        elif min_value/2 < value <= min_value:
            re_value = min_value


        if -min_value < value <= 0:
            re_value = 0
        elif -min_value < value <= -min_value/2:
            re_value = -min_value
        return re_value

    def change_direction_case(self,motor,motor_velocity):
        """
        :param motor: float
        :param motor_velocity: float
        :return:
        """

        if (motor.speed_array[self.m-1]<0 and motor_velocity>0) or (motor.speed_array[self.m-1]>0 and motor_velocity<0):
            motor_velocity_result1 = 0
        else:
            motor_velocity_result1=motor_velocity

        return motor_velocity_result1


    def next_speed_RiseFall_constraint(self,motor,motor_velocity):
        """
        :param motor: Motor object
        :param motor_velocity: float
        :return:
        """
        if (motor.speed_array[self.m-1] >= motor.minimum_speed) or (motor.speed_array[self.m-1] <= -motor.minimum_speed) :

            bool=True
            if abs(motor_velocity) > abs(motor.speed_array[self.m-1]) \
                    and (motor_velocity - motor.speed_array[self.m-1])/motor.speed_array[self.m-1]>0.2:
                motor_velocity_result2 = motor.speed_array[self.m-1] * 1.2
                bool=False

            if motor.minimum_speed/2 < abs(motor_velocity) < abs(motor.speed_array[self.m-1]) \
                    and (motor.speed_array[self.m-1] - motor_velocity)/motor.speed_array[self.m-1]>0.4:
                motor_velocity_result2 = motor.speed_array[self.m-1] * 0.8
                bool=False

            if bool==True:
                motor_velocity_result2=motor_velocity
        else:
            motor_velocity_result2 = motor_velocity

        return motor_velocity_result2

    def from_stop_motor_velocity(self, motor,motor_velocity):
        """

        :param motor: Motor object
        :param motor_velocity: float
        :return:
        """
        if motor.speed_array[self.m-1] == 0:
            if motor_velocity > 0:
                motor_velocity_result3 = motor.minimum_speed
            elif motor_velocity < 0:
                motor_velocity_result3 = - motor.minimum_speed
            else:
                motor_velocity_result3 = motor_velocity
        else:
            motor_velocity_result3 = motor_velocity

        return motor_velocity_result3

    def assign_minimum_velocity_speed(self,motor,motor_velocity):
        """

        :param motor: Motor object
        :param motor_velocity: float
        :return:
        """
        if motor_velocity > 0:
            motor_velocity_result4 = motor.minimum_speed
        elif motor_velocity < 0:
            motor_velocity_result4 = -motor.minimum_speed
        else:
            motor_velocity_result4 = motor_velocity

        return motor_velocity_result4

    def DeadlockLookUpTable_endEff(self, own_targ, neig_pos,  neig_targ):
        """
        Based on calculated distance between the two part of the two considered positioners, create a list of order
        distance (distance of target and pos), and compare this list to a lookup table containing configuration of if
        this positioners might be a deadlock or not

        NB: the look up table was created by pure logic, determining if 2 points placed in one axis (1D) with their
            corresponding targets (so 4 points in one axis) are in deadlock or not
            Since the positioners are in a focal plane (2D), then this approximation of a deadlock in 1D is used more
            than one time in the code (need two axis for 2D right...)

        :param own_targ: numpy array
            vector " agent positioner end-effector to target "
        :param neig_pos: numpy array
            vector " agent positioner end-effector to neigh pos end-effector "projected  on "own_targ" vector
        :param neig_targ: numpy array
            vector " agent positioner end-effector to neigh pos target" projected  on "own_targ" vector
        :return:

            flag_goodConf: boolean
                If it is a deadlock (possibly), then return false !
        """

        # evaluate distance norm:
        dist_own_targ = LA.norm(own_targ)

        ## CREATION OF THE ORDER DISTANCE TABLE BETWEEN TARGETS AND POSITIONER
        ## EVALUATION OF the different possible configuration (where the neighbour end-eff, target and agent
        # target are w.r.t to agent joint part) based on the measured distance

        # Assignement:
        # the method here is to attribute a "number" to them(target and part of the joint that are considered),
        #  from 1 to 4 and as hypothese that 1 "match" with  3 and 2 with 4
        # positioner being 1 (agent) and 2 (neigh)
        # respective target 3 (agent's target) and 4 (neigh's target)

         # don't need to put the "own target" since it is already present in the order_table below
        vec_table=[neig_pos,neig_targ]

        order_table=[1,3]
        order_distance=[0,dist_own_targ]

        U_own_targ=own_targ/LA.norm(own_targ)

        j=0
        for v_eph in vec_table:

            Uv_eph=v_eph/LA.norm(v_eph)

            # take its own target (U_own_targ) as reference
            if((np.sign(U_own_targ)[0] == np.sign(Uv_eph)[0]) and np.sign(U_own_targ)[1] == np.sign(Uv_eph)[1]):
                Eph_sign = 1
            else:
                Eph_sign = -1

            i=0
            biggest=False
            for dist_eph in order_distance:

                varDistEph = Eph_sign*LA.norm(v_eph)
                if(varDistEph < dist_eph):

                    biggest=True

                    if j==0: # means that v_eph is the neig_pos, so we add 2
                        order_table.insert(i,2)
                        order_distance.insert(i,varDistEph)
                        break
                    if j==1: # means that v_eph is the neig_targ, so we add 4
                        order_table.insert(i,4)
                        order_distance.insert(i,varDistEph)
                        break

                i=i+1

            if biggest==False:
                if j==0: # means that v_eph is the neig_pos, so we add 2
                    order_table.insert(i,2)
                    order_distance.insert(i,varDistEph)
                if j==1: # means that v_eph is the neig_targ, so we add 4
                    order_table.insert(i,4)
                    order_distance.insert(i,varDistEph)

            j=j+1

        ## COMPARISON with a look up table !
        #TODO: think of a better way to charaterize a deadlock  ?
        # positions where there is no deadlock
        position_GoodfRef_table = [[1,2,3,4],\
                                  [2,1,4,3],\
                                  [3,1,2,4],\
                                  [4,2,1,3],\
                                  [4,3,2,1],\
                                  [3,4,1,2],\
                                  [4,2,3,1],\
                                  [3,1,4,2],\
                                  [1,3,4,2],\
                                  [2,4,3,1],\
                                  [2,4,1,3],\
                                  [1,3,2,4] ]

        # positions where there is a deadlock
        position_BadfRef_table = [[1,2,4,3],\
                                  [2,1,3,4],\
                                  [3,2,1,4],\
                                  [4,1,2,3],\
                                  [3,4,2,1],\
                                  [4,3,1,2],\
                                  [3,2,4,1],\
                                  [4,1,3,2],\
                                  [1,4,3,2],\
                                  [2,3,4,1],\
                                  [1,4,2,3],\
                                  [2,3,1,4]]
        flag_goodConf=None

        for conf in position_GoodfRef_table:
            if (conf ==order_table):
                flag_goodConf=True
                break

        if flag_goodConf!=True:
            for conf in position_BadfRef_table:
                if (conf ==order_table):
                    flag_goodConf=False
                    break

        # if return false then it is a dealock !
        return flag_goodConf

    def DeadlockNotLookUpTable_notEndEff(self,qa,qn,qn_end,qt):
        """

        create a vector "qa to qn" AND based on this vector, see if end-effector of neigh positioner and its target are
        on same side or other side w.r.t to "qn" point

        ## CALCULATE FOR NON END-EFFECTOR TEST
        # if the part of the positioner which is possibly responsible of the deadlock is not the end-effector
        # then we create a another "point of view" from this joint to the other joint of the other positioner responsible
        # of the possible deadlock.
        # What is called a Point of view Here is basically constructing a vector to verify if the end-effector of the
        # other positioner is same side as its target w.r.t to the q_agent_ephClosest or q_neigh_ephClosest.
        # Then we first test this point of view before testing the point of view of from
        # the end-effector as calculated above

        :param qa, qn: numpy array
            vector with the coordinate in focal plane of Joint of the second arm of the agent and
            neighbour positioner respectively
        :param qn_end: numpy array
            vector with the coordinate in focal plane of end-effector positoiner neigh
        :param qt: numpy array
            vector with the target coordinate in focal plane of end-effector ppositioner neigh
        :return:
            flag_same_side: boolean
                if false, then deadlock possibly
        """

        sameSide = None

        v_qa_qn = qn-qa
        v_qa_qn_end = qn_end-qa
        v_qa_qt = qt-qa

        # projection on reference vector taken as v_qa_qn
        v_qn_projQn_end = (np.dot(v_qa_qn_end,v_qa_qn)/(LA.norm(v_qa_qn)**2))*v_qa_qn
        v_qn_projQn_t = (np.dot(v_qa_qt,v_qa_qn)/(LA.norm(v_qa_qn)**2))*v_qa_qn

        # Unit vector
        Uv_qn_projQn_end = v_qn_projQn_end/LA.norm(v_qn_projQn_end)
        Uv_qn_projQn_t = v_qn_projQn_t/LA.norm(v_qn_projQn_t)

        flag_same_side=None

        if np.sign(Uv_qn_projQn_end[0])!=np.sign(Uv_qn_projQn_t[0]) and \
                np.sign(Uv_qn_projQn_end[1])!=np.sign(Uv_qn_projQn_t[1]):
            flag_same_side=False
        else:
            flag_same_side=True

        # if return false, then possibly a deadlock
        return flag_same_side

    def deadlock_part_positioners(self, q_a_joint,q_a_middle, q_a_end,q_n_joint,q_n_middle, q_n_end):

        """
        Determine which part of the of two positioners (q_end, q_middle, q_joint) are the closest

        :param q_a_joint, q_a_middle, q_a_end: numpy array
             vector containing the coordinate on focal plane of the agent positioner's second arm
        :param q_n_joint, q_n_middle, q_n_end: numpy array
            vector containing the coordinate on focal plane of the neighbour positioner's second arm
        :return: numpy array
            return the two part of the positiones that are the closest

        """

        q_Table_agent = [ q_a_joint,q_a_middle, q_a_end]
        q_Table_neig = [q_n_joint,q_n_middle, q_n_end]

        min_dist_posPart = LA.norm(q_a_end-q_n_end)

        q_agent_eph1=q_a_end
        q_neigh_eph1=q_n_end

        for qa_Eph in q_Table_agent:
            for qn_Eph in  q_Table_neig:
                if( LA.norm(qa_Eph-qn_Eph) < min_dist_posPart):
                    min_dist_posPart = LA.norm(qa_Eph-qn_Eph)
                    q_agent_eph1 = qa_Eph
                    q_neigh_eph1 = qn_Eph

        return q_agent_eph1,q_neigh_eph1

    def deadlock_scenario(self,q_agent_ephClosest,q_neigh_ephClosest,q_t_a,q_t_n,q_a_end,q_n_end):

        """
        For the evalution of if two positioners are in a deadlock situtation, look for cases where it is possible that
        deadlocks are generated by the endeffector or not the end-effector of the positioners (if the closest part of the
        two positioners that is possibly generating the deadlock  is either a end-effector or not )

        :param q_agent_ephClosest, q_neigh_ephClosest:   numpy array

        :param q_t_a, q_t_n:  numpy array
            vector containing the coordinate on focal plane of the target of the positioner
            and its neighbours respectively
        :param q_a_end, q_n_end:  numpy array
            vector containing the coordinate on focal plane of the agent and neigh positioner's
            end effector of second arm

        :return: flag_sameSide_wrtToagent_notEndEff,flag_ptw_from_agent_end: boolean

        """
        flag_ptw_from_agent_end=None
        flag_sameSide_wrtToagent_notEndEff=None

        if(q_agent_ephClosest[0]==q_a_end[0] and q_agent_ephClosest[1]==q_a_end[1]) :
            # test: point of view from the agent: create vector agent end-eff to its target
            #       and project on it the end-eff of neigh and target neig

            #vector:
            v_End1_t1 = q_t_a-q_a_end

            v_End1_End2 = q_n_end - q_a_end
            v_End1_t2 =  q_t_n - q_a_end

            # projection on v_End1_t1
            V_End1_Proj_q_n_end = (np.dot(v_End1_End2,v_End1_t1)/(LA.norm(v_End1_t1)**2))*v_End1_t1
            V_End1_Proj_q_t_n = (np.dot(v_End1_t2,v_End1_t1)/(LA.norm(v_End1_t1)**2))*v_End1_t1

            flag_ptw_from_agent_end=self.DeadlockLookUpTable_endEff( v_End1_t1, V_End1_Proj_q_n_end,  V_End1_Proj_q_t_n)

        else:
            flag_sameSide_wrtToagent_notEndEff = self.DeadlockNotLookUpTable_notEndEff(q_agent_ephClosest,\
                                                                                       q_neigh_ephClosest,q_n_end,q_t_n)

        return flag_sameSide_wrtToagent_notEndEff,flag_ptw_from_agent_end

    def deadlock_evaluation(self,q_t_a,q_t_n, q_a_joint,q_a_middle, q_a_end,q_n_joint,q_n_middle, q_n_end, neighbour):
        """
        Evaluation if the configuration between two positioners is a deadlock or not.
            - Call a function "deadlock scenario" which will look if the conflicting situation is created by the end-eff
              then will look at a list of predefined and simplified situation if the situation is a deadlock
            - Then evaluate if deadlock or not depending on the result (of the point above)for the target
              and its neighbour

        :param q_t_a, q_t_n: numpy array
                vector containing the coordinate on focal plane of the target of the positioner
                and its neighbours respectively
        :param q_a_joint, q_a_middle, q_a_end: numpy array
                vector containing the coordinate on focal plane of the agent positioner's second
                arm (q_i_joint, q_i_middle, q_i_end)
        :param q_n_joint, q_n_middle, q_n_end: numpy array
                vector containing the coordinate on focal plane of the neighbour positioner's
                second arm (q_i_joint, q_i_middle, q_i_end)
        :param neighbour: positioner object
            positioner that is the neighbour of the current positioner (agent positioner)

        :return:
         flag_dealock: boolean
            if true, then deadlock, otherwise it is not consider as a deadlock by the algorithm
        """
        # Todo: probably not the best way of identifying a deadlock and correcting it...

        ## EVALUATE the distance between which part of the positioner are responsible of the deadlock (the two closest)
        q_agent_ephClosest,q_neigh_ephClosest =self.deadlock_part_positioners(q_a_joint,q_a_middle, q_a_end,\
                                                                              q_n_joint,q_n_middle, q_n_end)

        ## ESTIMATE FOR END-EFFECTOR OR NON END-EFFECTOR if they are possibly in a deadlock configuration or not !
        #from target point of view
        flag_ptw_from_agent_end=None
        flag_sameSide_wrtToagent_notEndEff=None
        flag_sameSide_wrtToagent_notEndEff,flag_ptw_from_agent_end = self.deadlock_scenario(q_agent_ephClosest,\
                                                                                            q_neigh_ephClosest,
                                                                                       q_t_a,q_t_n,q_a_end,q_n_end)
        #from neigh point of view
        flag_ptw_from_neigh_end=None
        flag_sameSide_wrtToneigh_notEndEff=None
        flag_sameSide_wrtToneigh_notEndEff,flag_ptw_from_neigh_end = self.deadlock_scenario(q_neigh_ephClosest,\
                                                                                            q_agent_ephClosest,
                                                                                         q_t_n,q_t_a,q_n_end,q_a_end)

        ## TAKE THE RESULTING BOOLEAN VALUE FROM ABOVE TO TEST IF DEADLOCK OR NOT
        flag_dealock=False

        if(flag_ptw_from_agent_end==None and flag_ptw_from_neigh_end==None):
            if (flag_sameSide_wrtToneigh_notEndEff==False or flag_sameSide_wrtToagent_notEndEff==False):
                flag_dealock=True

        elif(flag_ptw_from_agent_end==None and flag_sameSide_wrtToneigh_notEndEff==None):
            if(flag_sameSide_wrtToagent_notEndEff==False or flag_ptw_from_neigh_end==False):
                flag_dealock=True

        elif(flag_sameSide_wrtToagent_notEndEff == None and  flag_ptw_from_neigh_end ==None):
            if(flag_ptw_from_agent_end==False or flag_sameSide_wrtToneigh_notEndEff==False):
                flag_dealock=True

        elif(flag_sameSide_wrtToagent_notEndEff == None and  flag_sameSide_wrtToneigh_notEndEff ==None):
            if(flag_ptw_from_agent_end==False or flag_sameSide_wrtToneigh_notEndEff==False):
                flag_dealock=True

        return flag_dealock


    def reset_status_FromDeadlockSituationToNone(self):
        """
        When one positioner status is reset to its defaut, reset also the neighbour positioner that this agent
        positioner has affected !

        :return:
        """
        if self.agent.statusDeadlock_sp>0:
            if self.agent.pseudo_priority>self.agent.OPS_priority: # means that it was the one who is influencing the other
                for neighbour in self.agent.get_neighbours():
                    if neighbour.statusDeadlock_sp>0 :
                        neighbour.statusDeadlock_sp=0
                        break
            else: # means that the agent is the one that is being influenced !
                for neighbour in self.agent.get_neighbours():
                    if neighbour.statusDeadlock_sp>0 :
                        neighbour.statusDeadlock_sp=0
                        neighbour.pseudo_priority=neighbour.OPS_priority
                        break

    def agent_status_change(self, q_a_joint,q_a_middle, q_a_end,q_n_joint,q_n_middle, q_n_end, neighbour):
        """
        As positioners's behavior are determined by their status and its neighbours status, this function takes into
        account different scenario to determine when and what to change on the AGENT positioner status

        NB: this list of condition was created by logic, envisionning possibles scenarios that can happen if positioners
            with this set of status configuration were to encounter a positioners with another set of status

        :param q_a_joint, param q_a_middle, q_a_end: numpy array
            vector containing the coordinate on focal plane of the agent positioner's second
            arm (q_i_joint, q_i_middle, q_i_end)
        :param q_n_joint, q_n_middle, q_n_end:  numpy array
            vector containing the coordinate on focal plane of the neighbour positioner's
            second arm (q_i_joint, q_i_middle, q_i_end)
        :param neighbour: Positioner object

        :return:

            - changed in the positioner status
            - flagC: boolean
                if true, then

        """

        #TODO: all this condition can be reduced to fewer lines... In addition, possible a better strategy can be found
        #TODO: In addition, probably not the best way of interaction between positioner ?!?

        # use in order to avoid the self.agent positioner to be influenced by other positioner that are within the agent
        # positioner and are not supposed to in
        flagC = True

        ## Some explanation
        # the second condition after the "or" "neighbour.statusForceAttr==OFF":
        # in order to allow the reaction in chain process, when a positioner
        # A already reached its target, so A.statusArr==OFF
        # and so it cannot influence the status of other neighbour positioner !
        # now if a neighbour positioner B of positioner A  with higher priority
        # come nearby, positioner A will move to let it pass
        # and thus A.statusForceAttr==OFF...
        # However, if by moving, the positioner A has a chance to collide against
        # another neighbour C, then  positioner C will also need to have its status
        # change according to positioner A even if A.statusArr==OFF,
        # thus we have the second condition "neighbour.statusForceAttr==OFF" !
        if neighbour.statusArr==ON or neighbour.statusForceAttr==OFF:

            if self.agent.statusForceAttr==OFF and neighbour.statusForceAttr==OFF:

              if self.agent.statusDeadlock_sp>0 and neighbour.statusDeadlock_sp > 0:
                if self.agent.statusMemPr < neighbour.statusMemPr:

                    self.agent.statusMemID=neighbour.ident
                    self.agent.statusMemPr=neighbour.statusMemPr

                    flagC = False

              if self.agent.statusDeadlock_sp>0 and neighbour.statusDeadlock_sp == 0:
                if self.agent.OPS_priority < neighbour.statusMemPr:

                    self.reset_status_FromDeadlockSituationToNone()
                    self.agent.statusDeadlock_sp=0
                    self.agent.pseudo_priority=self.agent.OPS_priority

                    self.agent.statusMemID=neighbour.ident
                    self.agent.statusMemPr=neighbour.statusMemPr

                    flagC = False

              if self.agent.statusDeadlock_sp==0 and neighbour.statusDeadlock_sp > 0:
                if self.agent.statusMemPr < neighbour.OPS_priority:

                    self.agent.statusMemID=neighbour.ident
                    self.agent.statusMemPr=neighbour.OPS_priority

                    flagC = False

              if self.agent.statusDeadlock_sp==0 and neighbour.statusDeadlock_sp == 0:
                  if self.agent.statusMemPr < neighbour.statusMemPr:

                    self.agent.statusMemID=neighbour.ident
                    self.agent.statusMemPr=neighbour.statusMemPr

                    flagC = False

            elif self.agent.statusForceAttr==ON and neighbour.statusForceAttr==ON:

                if self.agent.statusDeadlock_sp>0 and neighbour.statusDeadlock_sp >0:
                  if self.agent.pseudo_priority < neighbour.pseudo_priority:

                    self.agent.statusForceAttr=OFF
                    self.agent.statusMemID=neighbour.ident
                    self.agent.statusMemPr=neighbour.pseudo_priority

                    flagC = False

                if self.agent.statusDeadlock_sp>0 and neighbour.statusDeadlock_sp == 0:
                  if  self.agent.OPS_priority < neighbour.OPS_priority:

                    self.reset_status_FromDeadlockSituationToNone()
                    self.agent.statusDeadlock_sp=0
                    self.agent.pseudo_priority=self.agent.OPS_priority

                    self.agent.statusForceAttr=OFF
                    self.agent.statusMemID=neighbour.ident
                    self.agent.statusMemPr=neighbour.pseudo_priority

                    flagC = False

                if self.agent.statusDeadlock_sp==0 and neighbour.statusDeadlock_sp>0:
                  if self.agent.OPS_priority < neighbour.OPS_priority :

                    self.agent.statusForceAttr=OFF
                    self.agent.statusMemID=neighbour.ident
                    self.agent.statusMemPr=neighbour.OPS_priority

                    flagC = False

                if self.agent.statusDeadlock_sp==0 and neighbour.statusDeadlock_sp == 0:
                  if self.agent.OPS_priority < neighbour.OPS_priority:

                    self.agent.statusForceAttr=OFF
                    self.agent.statusMemID=neighbour.ident
                    self.agent.statusMemPr=neighbour.pseudo_priority

                    flagC = False

            elif self.agent.statusForceAttr==ON and neighbour.statusForceAttr==OFF:

                if self.agent.statusDeadlock_sp>0 and neighbour.statusDeadlock_sp> 0:
                  if self.agent.pseudo_priority < neighbour.statusMemPr :

                    self.agent.statusForceAttr=OFF
                    self.agent.statusMemID=neighbour.ident
                    self.agent.statusMemPr=neighbour.statusMemPr

                    flagC = False

                if self.agent.statusDeadlock_sp>0 and neighbour.statusDeadlock_sp == 0:
                  if self.agent.OPS_priority < neighbour.statusMemPr :

                    self.reset_status_FromDeadlockSituationToNone()
                    self.agent.statusDeadlock_sp=0
                    self.agent.pseudo_priority=self.agent.OPS_priority

                    self.agent.statusForceAttr=OFF
                    self.agent.statusMemID=neighbour.ident
                    self.agent.statusMemPr=neighbour.statusMemPr

                    flagC = False

                if self.agent.statusDeadlock_sp==0 and neighbour.statusDeadlock_sp >0:
                  if self.agent.OPS_priority < neighbour.OPS_priority :

                    self.agent.statusForceAttr=OFF
                    self.agent.statusMemID=neighbour.ident
                    self.agent.statusMemPr=neighbour.OPS_priority

                    flagC = False

                if self.agent.statusDeadlock_sp==0 and neighbour.statusDeadlock_sp ==0:
                  if self.agent.OPS_priority < neighbour.statusMemPr :

                    self.agent.statusForceAttr=OFF
                    self.agent.statusMemID=neighbour.ident
                    self.agent.statusMemPr=neighbour.statusMemPr

                    flagC = False


            elif self.agent.statusForceAttr==OFF and neighbour.statusForceAttr==ON:

              if self.agent.statusDeadlock_sp>0 and neighbour.statusDeadlock_sp >0:
                  if self.agent.statusMemPr < neighbour.pseudo_priority :

                    self.agent.statusMemID=neighbour.ident
                    self.agent.statusMemPr=neighbour.pseudo_priority

                    flagC = False

              if self.agent.statusDeadlock_sp>0 and neighbour.statusDeadlock_sp ==0:
                  if self.agent.OPS_priority < neighbour.OPS_priority:

                    self.reset_status_FromDeadlockSituationToNone()
                    self.agent.statusDeadlock_sp=0
                    self.agent.pseudo_priority=self.agent.OPS_priority

                    self.agent.statusMemID = neighbour.ident
                    self.agent.statusMemPr = neighbour.OPS_priority

                    flagC = False

              if self.agent.statusDeadlock_sp==0 and neighbour.statusDeadlock_sp >0:
                  if self.agent.statusMemPr < neighbour.OPS_priority:

                    self.agent.statusMemID = neighbour.ident
                    self.agent.statusMemPr = neighbour.OPS_priority

                    flagC = False

              if self.agent.statusDeadlock_sp==0 and neighbour.statusDeadlock_sp ==0:
                  if self.agent.statusMemPr < neighbour.OPS_priority:

                    self.agent.statusMemID = neighbour.ident
                    self.agent.statusMemPr = neighbour.OPS_priority

                    flagC = False

            # When agent is inside the zone of influence of the neighbour that affected its statusForceAttr,
            # then if this agent come back into this loop, it won't go into any of the condition if
            #  above, thus the flacC will be True and thus it's status will be changed to ON even though
            #  it is still affected by its neighbour (creation of an oscillation ON/OFF)...
            #  which is not we want so we add this condition
            if self.agent.statusMemID == neighbour.ident:
                flagC = False

        if flagC == False:
            return False
        else:
            return True

    def repulsive_force_between(self,neighbour,q1,q2,distance,end,round1,round2,round1_end,round2_end,h):
        """
        Calculate the coefficient of the repulsive force

        :param neighbour: Positioner object
        :param q1,q2: numpy array
        :param distance: float
        :param end: boolean
            used to distinguish when using the weighted distance and the normal one
        :param round1, round2, round1_end, round2_end: numpy array
        :param h: int
        :return:
        """

        global i_Joint
        r_init = self.agent.ident % 7 + 1
        r = r_init/7 + 1
        bigD = (neighbour.OPS_priority+6.5)/2
        #prim = (neighbour.priority/2)+2
        rprim = 30

        if h==1:
            prim2= ((neighbour.OPS_priority)+CSTE_FORCE)
        if h==2:
            prim2= ((neighbour.OPS_priority)+CSTE_FORCE/2)

        # when using weight-distance
        if end == False:
            if distance < rprim:
                delta_sai = (((bigD)**2 - self.d**2) / (distance**2 - self.d**2)**2 * (q2-q1))

                self.gcol1 += np.dot(delta_sai, round1)
                self.gcol2 += np.dot(delta_sai, round2)

        # when using normal distance (euclidian)
        if end == True:
            if neighbour.statusArr==OFF :

                # added this additional condition "if" in order to make the chain reaction of positioner possible
                # explanation:
                # when the positioner A reached its target, A.statusArr==OFF, then if a higher priority positioner B
                # comes near by the positioner A, A.statusForceAttr==OFF and then positioner A will move to let B pass
                # now what if when positioner A moves and goes against a neighbour positioner C,
                # then because A.statusArr==OFF
                # A NORMALLY won't be able to affect the status of C (but changes in the code in
                # def "agent_status_change"were made such that it is possible). anyway the problem is that,
                # because now it affect the positioner C by chain reaction because of B through A, the zone of
                # influence of the force should be the same as when a positioner could affect the status of an
                #  another positioner (even though statusArr==OFF)
                if neighbour.statusForceAttr==OFF:
                    condition_zone= STATUS_ON_ZONE_INFLUENCE
                else:
                    condition_zone=STATUS_OFF_ZONE_INFLUENCE
            else:
                condition_zone= STATUS_ON_ZONE_INFLUENCE


            if distance < condition_zone :
                delta_sai = (prim2**2 - self.d2**2) / (distance**2 - self.d**2)**2 * (q2-q1)

                self.gcol_end1 += np.dot(delta_sai,round1_end)

                if i_Joint ==0:
                    self.gcol_end2 += np.dot(delta_sai,round2_end)


    def firstDepth_StatusChange_ChainPositioner(self,pseudo_agent):
        """

        first depth search algorithm in order to reset all the other positioners states ON, ID and Pr,
        whose states were affected by the current positioner that has the highest priority

        :param pseudo_agent: Positioner object
            current positioner that affected the state of an another positioner
        :return:


        #NB:
        # This function is to avoid endless status change through the chain reaction of positioners:
        # If the agent positioner A is inside its neighbour positioner B that changed its status, and positioner B
        # status was changed by another one. If positioner B.statusForceAttr becomes ON but has lower priority than A
        # then it can cause problem: positioner 3 can be [0,1,4,1] (influenced by positioner 1) and
        # positioner 2 [0,3,4,1] which was influenced by positioner 3 (chain reaction influence),
        # then when positioner 1 is out, we have for 3 [1,0,0,1] and 2 [0,3,4,1] but positioner 3
        # will be influenced by 2 thus 3: [0,2,4,1]... so 2 influence 3 and 3 influence 2...

        """

        # depth first search is used in order to :
        for neighbourDFS in pseudo_agent.get_neighbours():
            if neighbourDFS.statusMemID==pseudo_agent.ident and neighbourDFS.statusForceAttr==OFF:
                neighbourDFS.statusForceAttr= ON
                neighbourDFS.statusMemID= 0
                neighbourDFS.statusMemPr= 0
                self.firstDepth_StatusChange_ChainPositioner(neighbourDFS)
                break

    def localize_deadlock(self,q_t_a,q_t_n,q_i_joint,q_i_middle, q_i_end,q_j_joint,q_j_middle, q_j_end, neighbour):

        """

        Localise deadlock produced by two positioners with the same priority:
            - look if the problematic configuration between two postioners is really a deadlock or
              just a situation where they are temporarly blocked with the function "scenario_deadlock" listing
              all sort of possible deadlock situations
            - If deadlock and depending on distance:
                - change their statusDeadlock_sp, isolate the two problematics positioners where additional action/behavior
                  are used to solve the conflict
                - update their pseudo_priority, ie their pseudo-priority only used to solve the conflicting situation
                - update the repulsiveTarget of the positioners, ie the target of neighbour positioner
                  that will be use to create a isotropic and constant repulsive force on the agent positioners
                  (also toadd additional force to bypass the local minima)

        :param q_t_a, q_t_n:  numpy array
            vector containing the coordinate on focal plane of the target of the
            positioner and its neighbours respectively
        :param q_i_joint, q_i_middle, q_i_end: numpy arrays
            vector containing the coordinate on focal plane of the positioner's
            second arm (q_i_joint, q_i_middle, q_i_end)
        :param q_j_joint,  q_j_middle, q_j_end:  numpy arrays
            vector containing the coordinate on focal plane of the neighbour positioner's
            second arm (q_i_joint, q_i_middle, q_i_end)

        :param neighbour: positioner object
            positioner that is the neighbour of the current positioner (agent positioner)
        :return:
        """



        if self.agent.statusForceAttr==ON and neighbour.statusForceAttr==ON :
            if self.agent.OPS_priority == neighbour.OPS_priority: # if have the same priority
                if self.agent.statusDeadlock_sp==0 and neighbour.statusDeadlock_sp==0: # on the original layer

                    distance1=LA.norm(q_j_joint-q_i_end)
                    distance2=LA.norm(q_j_middle-q_i_end)

                    if (math.fabs(self.agent.MeanVelMotor[1])< LIMITE_FOR_NOISE and\
                      math.fabs(neighbour.MeanVelMotor[1])< LIMITE_FOR_NOISE) or\
                           (distance1< STATUS_ON_ZONE_INFLUENCE and distance2<STATUS_ON_ZONE_INFLUENCE) :

                        #agent_priority=changed=False
                        flag_deadlock = self.deadlock_evaluation(q_t_a,q_t_n,q_i_joint,q_i_middle, q_i_end,\
                                                                 q_j_joint,q_j_middle, q_j_end, neighbour)

                        if(flag_deadlock==True):

                            # distance from end to target
                            dist_aEnd_at = LA.norm(q_i_end - q_t_a)
                            dist_nEnd_nt = LA.norm(q_j_end - q_t_n)

                            agent_priority=None
                            if (dist_aEnd_at<dist_nEnd_nt):
                                agent_priority=False
                            else:
                                agent_priority=True

                            if agent_priority != None:
                                self.agent.statusDeadlock_sp=self.agent.ident
                                neighbour.statusDeadlock_sp=self.agent.ident

                                if agent_priority==True:
                                    self.agent.pseudo_priority=self.agent.OPS_priority+1
                                    neighbour.repulsiveTarget = np.array((self.agent.x_fibre_focal,\
                                                                          self.agent.y_fibre_focal))
                                else:
                                    neighbour.pseudo_priority=neighbour.OPS_priority+1
                                    self.agent.repulsiveTarget = np.array((neighbour.x_fibre_focal,\
                                                                           neighbour.y_fibre_focal))


    def add_repulsive_force(self,q_i_joint,q_i_middle,q_i_end,r_beta,round1,round2,round1_end,round2_end,l1,l2):

        """
        Call function to calculate the repulsive force and the change of its status w.r.t to the positioner's neighbour

        :param q_i_joint: numpy array
            vector on the focal plane of q_i_joint (second arm)
        :param q_i_middle: numpy array
            vector on the focal plane of q_i_middle (second arm)
        :param q_i_end: numpy array
            vector on the focal plane of q_i_end (second arm)
        :param r_beta: float
             angle:half a turn back to compensate for beta rotating from tacked in position
        :param round1, round2:  numpy array
             chain derivatives of position to rotation of the middle joint of the beta arm
        :param round1_end , round2_end:  numpy array
             chain derivatives of position to rotation of the joints (middle of beta arm, fiber position)
        :param l1: float
             length of first arm
        :param l2: float
             length of the second arm

        :return:
            - calculate the coefficient of the repulsive force (assigned to the motion planning's class attribute )
            - change the status of the positioners (assigned to the positioner's class attribute)
        """

        # Repulsive forces from the neighbours (19 neighbours)
        flagStatus=True

        for neighbour in self.agent.get_neighbours():

            rj1 = neighbour.motor1.position_array[self.m]
            rj2 = neighbour.motor2.position_array[self.m]

            # angle of the beta arm
            r_beta_j = rj1 + rj2 + math.pi

            q_j_joint,q_j_middle,q_j_end  =self.get_coordinate_second_arm_focal_plan(neighbour,r_beta_j,l2)

            q_t_a= np.array((self.agent.x_fibre_focal, self.agent.y_fibre_focal))
            q_t_n= np.array((neighbour.x_fibre_focal,neighbour.y_fibre_focal))

            # Used of weighted_distance
            """
            # NB:
            # Use a rectangular shape to represent the zone of influence of the positioner
            # I use instead a circle with q_joint, q_end and q_middle as its center !

            self.weigthed_distance(q_j_end,q_i_middle,r_beta,l1,l2)

            # add the repulsive force from the other positioners:
            # It is computationally more efficient to take three points of the other
            # positioners and calculate the repulsive force
            end=False

            distance= self.weigthed_distance(q_j_end,q_i_middle,r_beta,l1,l2)
            #self.repulsive_force_between(neighbour,q_j_end,q_i_middle, distance,end,\
                                            round1,round2,round1_end,round2_end)

            distance= self.weigthed_distance(q_j_middle,q_i_middle,r_beta,l1,l2)
            #self.repulsive_force_between(neighbour,q_j_middle,q_i_middle, distance,end,\
                                            round1,round2,round1_end,round2_end)

            distance= self.weigthed_distance(q_j_joint,q_i_middle,r_beta,l1,l2)
            #self.repulsive_force_between(neighbour,q_j_joint,q_i_middle, distance,end,\
                                            round1,round2,round1_end,round2_end)
            """

            # Use of normal distance
            end=True

            # make it start but not at the beginning of the program,
            # because the positioner are not moving at the beginning of the program,
            #  it could be interpreted as a deadlock...
            if self.m>10 and self.agent.statusArr==ON:
               self.localize_deadlock(q_t_a,q_t_n,q_i_joint,q_i_middle, q_i_end,q_j_joint,q_j_middle, q_j_end, neighbour)

            ##STOPPING CRITERIA
            ## if the two positionres are in deadlock situtation and the positioner consider is the one having priority
            if self.agent.statusDeadlock_sp==neighbour.statusDeadlock_sp and self.agent.pseudo_priority>self.agent.OPS_priority:

                # TODO: INSTEAD OF ALWAYS RE-SEARCHING WHICH PART OF THE POS CREATE THE DEADLOCK,
                # todo: THE FIRST TIME WE SEARCHED FOR IT IN "deadlock_evaluation", SAVE IT ?!?

                flag_deadlock = self.deadlock_evaluation(q_t_a,q_t_n,q_i_joint,q_i_middle, q_i_end,\
                                                         q_j_joint,q_j_middle, q_j_end, neighbour)
                q_agent_ephClosest,q_neigh_ephClosest =self.deadlock_part_positioners(q_i_joint,q_i_middle, q_i_end,\
                                                                                      q_j_joint,q_j_middle, q_j_end)
                distanceNearest=LA.norm(q_agent_ephClosest-q_neigh_ephClosest)

                # if situation is not consider deadlock anymore "or"
                # if the two positioners already arrived at their target
                if(flag_deadlock==False or (self.agent.statusArr==OFF and neighbour.statusArr==OFF) \
                   or distanceNearest > CSTE_MARGIN_STATUSZONE_STOPPING_DEADLOCK+STATUSZONE):

                    self.agent.statusDeadlock_sp=0
                    neighbour.statusDeadlock_sp=0
                    self.agent.pseudo_priority=self.agent.OPS_priority

                    self.agent.repulsiveTarget=np.array((0,0))
                    neighbour.repulsiveTarget=np.array((0,0))

                    if(neighbour.statusMemID==self.agent.ident and neighbour.statusForceAttr==OFF):
                        neighbour.statusForceAttr=ON

            ## CALCULATE IF POSITIONERS ARE WITHIN "STATUSZONE" SO THEY CAN INFLUENCE EACH OTHER STATUS

            qTable=[q_i_joint,q_i_middle,q_i_end]
            qTableNeig=[q_j_joint,q_j_middle,q_j_end]

            # distance from middle to all 3 points of the positioners: use to put the changes of status of the agent
            for q_Eph in qTable:

                # distance of each agent positioner to end-eff neigh positioner
                distance=LA.norm(q_Eph-q_j_end)
                if self.agent.statusForceAttr==OFF and distance <STATUSZONE+CSTE_MARGIN_STATUSZONE :
                    flagStatus = self.agent_status_change(q_i_joint,q_i_middle,q_i_end,\
                                                          q_j_joint,q_j_middle,q_j_end,neighbour)
                elif distance < STATUSZONE:
                    flagStatus = self.agent_status_change(q_i_joint,q_i_middle,q_i_end,\
                                                          q_j_joint,q_j_middle,q_j_end,neighbour)

                # distance of each agent positioner to middle neigh positioner
                distance=LA.norm(q_Eph-q_j_middle)
                if self.agent.statusForceAttr==OFF and distance < STATUSZONE+CSTE_MARGIN_STATUSZONE:
                    flagStatus = self.agent_status_change(q_i_joint,q_i_middle,q_i_end,\
                                                          q_j_joint,q_j_middle,q_j_end,neighbour)
                elif distance < STATUSZONE:
                    flagStatus = self.agent_status_change(q_i_joint,q_i_middle,q_i_end,\
                                                          q_j_joint,q_j_middle,q_j_end,neighbour)

               # distance of each agent positioner to joint neigh positioner
               # distance= LA.norm(q_Eph-q_j_joint)
               # if self.agent.statusForceAttr==OFF and distance < STATUSZONE+CSTE_MARGIN_STATUSZONE:
               #     flagStatus = self.agent_status_change(q_i_joint,q_i_middle,q_i_end,\
               #                                           q_j_joint,q_j_middle,q_j_end,neighbour)
               # if distance < STATUSZONE:
               #     flagStatus = self.agent_status_change(q_i_joint,q_i_middle,q_i_end,\
               #                                            q_j_joint,q_j_middle,q_j_end,neighbour)


            # use to allow calculation of repulsive force on second and first arm in the
            # "repulsive_force_between" function
            #TODO: get rid of this global and find a better name ?!?...
            global i_Joint
            i_Joint=1
            #distance= LA.norm(q_i_joint - q_j_end)
            #self.repulsive_force_between(neighbour,q_j_end,q_i_joint, distance,end,\
                                        # round1,round2,round1_end,round2_end,2)

            #distance= LA.norm(q_i_joint - q_j_middle)
            #self.repulsive_force_between(neighbour,q_j_middle,q_i_joint, distance,end,\
                                        # round1,round2,round1_end,round2_end,2)

            #distance= LA.norm(q_i_joint - q_j_joint)
            #self.repulsive_force_between(neighbour,q_j_joint,q_i_joint, distance,end,\
                                        # round1,round2,round1_end,round2_end)

            i_Joint=0
            distance= LA.norm(q_i_end - q_j_end)
            self.repulsive_force_between(neighbour,q_j_end,q_i_end, distance,end,\
                                         round1,round2,round1_end,round2_end,1)

            distance= LA.norm(q_i_end - q_j_middle)
            self.repulsive_force_between(neighbour,q_j_middle,q_i_end, distance,end,\
                                         round1,round2,round1_end,round2_end,1)

            #distance= LA.norm(q_i_end - q_j_joint)
            #self.repulsive_force_between(neighbour,q_j_joint,q_i_end, distance,end,\
                                        # round1,round2,round1_end,round2_end,2)

            distance= LA.norm(q_i_middle - q_j_middle)
            self.repulsive_force_between(neighbour,q_j_middle,q_i_middle, distance,end,\
                                         round1,round2,round1_end,round2_end,1)


            ## REACH TARGET STATUS
            # condition in order to stop the agent.positioner to influence the ON/OFF of the other neighbour positioners
            # only useful to do it on the second arm (when reaches the target)
            # the range of value when the positioner arrives at the target (by measurement): [0.006 to 0.021]

            # NB: the second condition "elif" with value 0.035 is used so the change between statusArr= ON/OFF
            # is not so abrupt in order to avoid unnecessary oscillation in statusArr, thus no influence in statusForceAttr
            # The value of 0.065 is chosen by trial and error (normally 0.35 is ok too but let's have a good margin ;) )
            if  math.fabs(self.agent.motor2.position_array[self.m] - self.agent.motor2.position_array[-1]) < 0.025:
                self.agent.statusArr=OFF

            elif math.fabs(self.agent.motor2.position_array[self.m] - self.agent.motor2.position_array[-1]) > 0.035 \
                and self.agent.statusForceAttr==ON:
                self.agent.statusArr=ON

            ##RESET STATUS TO DEFAUT
            # Now if the agent is not inside the neighbour that previously changed its statusForceAttr before,
            # then agent positioner is ON
            if flagStatus==True and self.agent.statusMemID==neighbour.ident:

                # see description on fps_shared.py
                #([ON/OFF status, ident of neigbourh influencer, priority of neighbour influencer,
                # can this agent positioner influenced the other])
                self.agent.statusForceAttr=ON
                self.agent.statusMemID=0
                self.agent.statusMemPr=0

                self.firstDepth_StatusChange_ChainPositioner(self.agent)

    def get_coordinate_second_arm_focal_plan(self,positioner,r_beta,l2):
        """
        :param r_beta: float
            half a turn back to compensate for beta rotating from tacked in position
        :param l2: float
            length of the second arm
        :return:  numpy array vector
            Coordinate of the second arm(q_joint,q_middle and q_end_effector) on focal plan/cartesian
        """

        q_i_end = np.array(positioner.motor_position_to_focal_cartesian(positioner.motor1.position_array[self.m],\
                                                                        positioner.motor2.position_array[self.m]))
        q_i_middle = q_i_end + (l2/2) * np.array((-math.cos(r_beta),-math.sin(r_beta)))
        q_i_joint = q_i_end + l2 * np.array((-math.cos(r_beta),-math.sin(r_beta)))

        return q_i_joint,q_i_middle,q_i_end

    def pas_config_dnf(self,time):
        """

        Calculate the potential field and the state of the finite state machine for the current positioner
        with recpect to its neighbouring positioners

        return the velocity of its two arms

        :param time: int
            Is the time of the application
        :return: floats
            the velocity of the motor 1 and 2
        """
        assert isinstance(self.agent, fps.Positioner)

        self.m = time # simulation time

        # initialized the three terms of the potential function
        self.ggoal1 = 0
        self.ggoal2 = 0

        self.gcol1 = 0
        self.gcol2 = 0

        self.gcol_end1 = 0   # added for protection of fiber
        self.gcol_end2 = 0


        r1,r2,l1,l2,r_beta = self.get_angle_length_data()

        q_i_joint,q_i_middle,q_i_end  =self.get_coordinate_second_arm_focal_plan(self.agent,r_beta,l2)

        # chain derivatives of position to rotation of the joints (middle of beta arm, fiber position)
        round1 = np.array((-l1*math.sin(r1) - l2/2*math.sin(r_beta), l1*math.cos(r1) + l2/2*math.cos(r_beta)))
        round2 = np.array((-l2/2*math.sin(r_beta), l2/2*math.cos(r_beta)))

        round1_end = np.array((-l1*math.sin(r1) - l2*math.sin(r_beta), l1*math.cos(r1) + l2*math.cos(r_beta)))
        round2_end = np.array((-l2*math.sin(r_beta), l2*math.cos(r_beta)))

        # Calculate potential field forces
        self.add_repulsive_force(q_i_joint,q_i_middle,q_i_end,r_beta,round1,round2,round1_end,round2_end,l1,l2)

        # Calculate the attractive force
        if self.agent.statusForceAttr==ON :
            self.add_attractive_force()

        # when agent and neigh positioners are in deadlock situation(statusDeadlock_sp) because they have originally the same
        # priority
        if self.agent.statusDeadlock_sp > 0 and self.agent.pseudo_priority == self.agent.OPS_priority :
            M1attractiveF  = 0
            M2attractiveF  = 0

            tmpVec1=(self.agent.repulsiveTarget - q_i_end)
            tmpVec2=(self.agent.repulsiveTarget - q_i_end)

            tmpVec1=tmpVec1/LA.norm(tmpVec1)
            tmpVec2=tmpVec2/LA.norm(tmpVec2)

            # add to agent positioner an another isotropic and constant repulsive force from neigh positioner target
            M1repulsiveF_NeighTarget = -((self.agent.motor1.maximum_speed+\
                                          self.agent.motor1.minimum_speed)/REPULSIVE_TARGET_FACTOR) * \
                                         np.dot(tmpVec1,round1_end)
            M2repulsiveF_NeighTarget = -((self.agent.motor2.maximum_speed+\
                                          self.agent.motor2.minimum_speed)/REPULSIVE_TARGET_FACTOR) * \
                                         np.dot(tmpVec2,round2_end)

        else:
            M1attractiveF  = -self.k11 * self.ggoal1 * ((self.agent.pseudo_priority+CONST_ATTRACTIVE_FORCE)/2)
            M2attractiveF  = -self.k21 * self.ggoal2 * ((self.agent.pseudo_priority+CONST_ATTRACTIVE_FORCE)/2)

            M1repulsiveF_NeighTarget = 0
            M2repulsiveF_NeighTarget = 0

        M1repulsiveFAll= self.k12 * self.gcol1
        M1repulsiveFEnd= self.k13 * self.gcol_end1

        M2repulsiveFAll= self.k22 * self.gcol2
        M2repulsiveFEnd= self.k23 * self.gcol_end2

        motor1_velocity = M1attractiveF + M1repulsiveFAll + M1repulsiveFEnd + M1repulsiveF_NeighTarget
        motor2_velocity = M2attractiveF + M2repulsiveFAll + M2repulsiveFEnd + M2repulsiveF_NeighTarget


        ### NOISE:
        # In order to avoid deadlock where one positioner is completely orthogonal to an other one, making the forces
        # of repulsion completely equal and thus deadlocks (and since the strategies of priority does not work here we
        # added some "noise" !)
        #  value chosen randomly and kind of validated with error and trial

        # if not in target
        if math.fabs(self.agent.motor2.position_array[self.m] - self.agent.motor2.position_array[-1]) > 0.05:
            if math.fabs(self.agent.MeanVelMotor[0])< LIMITE_FOR_NOISE:
                motor1_velocity = motor1_velocity+self.motor1Addition
            else:
                self.motor1Addition=self.agent.MeanVelMotor[0]+(np.sign(self.agent.MeanVelMotor[0])*LIMITE_FOR_NOISE)

            if math.fabs(self.agent.MeanVelMotor[1])< LIMITE_FOR_NOISE:
                motor2_velocity = motor2_velocity+self.motor2Addition
            else:
                self.motor2Addition=self.agent.MeanVelMotor[1] + (np.sign(self.agent.MeanVelMotor[1])*LIMITE_FOR_NOISE)


        if self.m!=0:
            #in case of change in the direction, stop for this timestep

            motor1_velocity = self.change_direction_case(self.agent.motor1,motor1_velocity)
            motor2_velocity = self.change_direction_case(self.agent.motor2,motor2_velocity)

            # check if the next speed rises/ fall more than 20%
            # TODO: After tests by David Atkinson, this might change to 40%

            motor1_velocity = self.next_speed_RiseFall_constraint(self.agent.motor1,motor1_velocity)
            motor2_velocity = self.next_speed_RiseFall_constraint(self.agent.motor2,motor2_velocity)

            #From stop motors can't begin with a sudden change
            motor1_velocity = self.from_stop_motor_velocity(self.agent.motor1,motor1_velocity)
            motor2_velocity = self.from_stop_motor_velocity(self.agent.motor2,motor2_velocity)

        #First time_step rules; From David attkinson on May 17, 2016

        elif self.m==0: # if first step of the simulation, give the motor the minimum speed

            motor1_velocity=self.assign_minimum_velocity_speed(self.agent.motor1,motor1_velocity)
            motor2_velocity=self.assign_minimum_velocity_speed(self.agent.motor2,motor2_velocity)


        motor1_velocity = self.saturate_max(motor1_velocity, self.agent.motor1.maximum_speed)
        motor1_velocity = self.saturate_min(motor1_velocity, self.agent.motor1.minimum_speed)

        motor2_velocity = self.saturate_max(motor2_velocity, self.agent.motor2.maximum_speed)
        motor2_velocity = self.saturate_min(motor2_velocity, self.agent.motor2.minimum_speed)


        # CALCULATE THE AVERAGE VELOCITY OVER A "wind_time" WINDOWS
        # calculation of the average is done in this pa_dnf file, motion planning class but the average
        # velocity of the two motors belongs to the positioner class in

        wind_time=6 # threshold in order to consider that two positioners are in deadlock situation 0.03

        if time>wind_time:

            self.agent.MeanVelMotor[0]=sum(self.motor1VelList)/wind_time
            self.agent.MeanVelMotor[1]=sum(self.motor2VelList)/wind_time

            self.motor1VelList.pop(0)
            self.motor2VelList.pop(0)
            self.motor1VelList.append(motor1_velocity)
            self.motor2VelList.append(motor2_velocity)

        else:
            self.motor1VelList.append(motor1_velocity)
            self.motor2VelList.append(motor2_velocity)


        return (motor1_velocity,motor2_velocity)
