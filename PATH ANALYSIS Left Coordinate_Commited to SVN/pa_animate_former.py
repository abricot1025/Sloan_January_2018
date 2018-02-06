try:
    import mocpath.conflict.fps_shared as fps
except ImportError:
    import fps_shared as fps

try:
    import mocpath.common.util as util
except ImportError:
    import util as util

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math

import pa_dnf as pa
import parameters as param

INS_POS_LENGTH1 = fps.INS_POS_LENGTH1     # Length of alpha arm in mm.
INS_POS_WIDTH1 = fps.INS_POS_WIDTH1       # Width of alpha arm in mm.
INS_POS_LENGTH2 = fps.INS_POS_LENGTH2       # Length of beta arm in mm (including metrology zone)
INS_POS_WIDTH2 = fps.INS_POS_WIDTH2         # Width of beta arm in mm (before metrology zone).
PITCH = INS_POS_LENGTH1 + INS_POS_LENGTH2

STATUS_ON_ZONE_INFLUENCE = pa.STATUS_ON_ZONE_INFLUENCE
STATUSZONE = pa.STATUSZONE

dt = param.dt
sim_length = param.sim_length

# Counters for making replay and delays in the simulation for recording
global Counter
Counter = -1

global counter2
counter2 = 0

def run_animation(agents):

    # Plottig variables
    x_focal_min = agents[0].x_centre_focal
    x_focal_max = agents[0].x_centre_focal
    y_focal_min = agents[0].y_centre_focal
    y_focal_max = agents[0].y_centre_focal

    for agent in agents:

        if agent.x_centre_focal < x_focal_min:
            x_focal_min = agent.x_centre_focal
        elif agent.x_centre_focal > x_focal_max:
            x_focal_max = agent.x_centre_focal

        if agent.y_centre_focal < y_focal_min:
            y_focal_min = agent.y_centre_focal
        elif agent.y_centre_focal > y_focal_max:
            y_focal_max = agent.y_centre_focal

    pad = 1.5 * PITCH
    fig = plt.figure(1)
    #Find the size of the plot
    plt.axis([x_focal_min-pad, x_focal_max+pad, y_focal_min-pad, y_focal_max+pad])

    plt.gca().set_aspect('equal', adjustable='box')
    subfig = fig.add_subplot(1,1,1)


    beta_polygon = []
    beta_circle1 = []
    beta_circle2 = []
    fiber = []

    end_effector_forceZone  = []
    end_effector_statusZone = []

    for agent in agents:
        # Draw the pich and the alpha motor and the target point
        firstcircle_p = np.array([agent.x_centre_focal,agent.y_centre_focal])

        pitchcircle = plt.Circle(firstcircle_p, radius= PITCH/2, color=[0.8,0.8,0.8], fill=False)

        motorcircle = plt.Circle(firstcircle_p, radius= INS_POS_LENGTH1, color=[0.4,0.4,0.4], fill=False)

        # make the circle for the targets
        thetag = agent.motor1.position_array[-1]
        phig = agent.motor2.position_array[-1]


        # make the circle for the targets from the target list
        targetcircle = plt.Circle(util.polar_to_cartesian(agent.r_fibre_focal,agent.theta_fibre_focal),radius=1, color='r', fill=True)
        subfig.add_patch(targetcircle)

        subfig.add_patch(motorcircle)
        subfig.add_patch(pitchcircle)


        # Get the theta and phi angels
        theta = agent.motor1.position_array[0]
        phi = agent.motor2.position_array[0] + math.pi

        # calculate the positions of the four edges of the polygon:
        # this polygon will represent the beta arm along with two circles:
        # the shape at the end is like an extended ellipse

        joint = firstcircle_p + INS_POS_LENGTH1 * np.array([math.cos(theta),math.sin(theta)])\
        + INS_POS_WIDTH2/2*np.array([math.cos(theta+phi),math.sin(theta+phi)])

        point1 = joint + INS_POS_WIDTH2/2*np.array([math.sin(theta+phi),-math.cos(theta+phi)])
        point2 = joint + INS_POS_WIDTH2/2*np.array([-math.sin(theta+phi),math.cos(theta+phi)])
        point3 = point2 + (INS_POS_LENGTH2 - INS_POS_WIDTH2/2)*np.array([math.cos(theta+phi),math.sin(theta+phi)])
        point4 = point1 + (INS_POS_LENGTH2 - INS_POS_WIDTH2/2)*np.array([math.cos(theta+phi),math.sin(theta+phi)])

        joint_polygon = plt.Polygon([point1,point2,point3,point4], edgecolor=[0,0.6,1], facecolor=[0,0.6,1])
        subfig.add_patch(joint_polygon)

        joint_circle1 = plt.Circle(joint,radius=INS_POS_WIDTH2/2,fill=True,edgecolor=[0,0.6,1], facecolor=[0,0.6,1])
        joint_circle2 = plt.Circle(joint+(INS_POS_LENGTH2 - INS_POS_WIDTH2/2)*np.array([math.cos(theta+phi),math.sin(theta+phi)])\
                                   ,radius=INS_POS_WIDTH2/2,fill=True,edgecolor=[0,0.6,1], facecolor=[0,0.6,1])
        fiber_circle = plt.Circle(joint+(INS_POS_LENGTH2 - INS_POS_WIDTH2/2)*np.array([math.cos(theta+phi),math.sin(theta+phi)])\
                                   ,radius=INS_POS_WIDTH2/4,fill=True,edgecolor=[0,1,1], facecolor=[0,1,1])

        # create circle to have a visual on the status zone (zone where positioner's status can be changed by others)
        # and on zone where repulsive force from other positioners is applied
        #fiber_circleForce = plt.Circle(joint+(INS_POS_LENGTH2 - INS_POS_WIDTH2/2)*np.array([math.cos(theta+phi),math.sin(theta+phi)])\
        #                           ,radius=STATUS_ON_ZONE_INFLUENCE,fill=False,edgecolor=[0,0,0.5], facecolor=[0,0,1])
        #fiber_circleStatus = plt.Circle(joint+(INS_POS_LENGTH2 - INS_POS_WIDTH2/2)*np.array([math.cos(theta+phi),math.sin(theta+phi)])\
        #                           ,radius=STATUSZONE,fill=False,edgecolor=[1,0,1], facecolor=[0,0,1])

        subfig.add_patch(joint_circle1)
        subfig.add_patch(joint_circle2)
        subfig.add_patch(fiber_circle)
        #subfig.add_patch(fiber_circleForce)
        #subfig.add_patch(fiber_circleStatus)


        beta_polygon.append(joint_polygon)
        beta_circle1.append(joint_circle1)
        beta_circle2.append(joint_circle2)
        fiber.append(fiber_circle)

        #end_effector_forceZone.append(fiber_circleForce)
        #end_effector_statusZone.append(fiber_circleStatus)

    def update(data):
        x,y,z = data
        for i in range(len(x)):
            beta_polygon[i].set_xy(x[i])
            beta_circle1[i].center = y[i]
            beta_circle2[i].center = z[i]
            fiber[i].center = z[i]

            #end_effector_forceZone[i].center=z[i]
            #end_effector_statusZone[i].center=z[i]

    def data_gen():
        x = []
        y = []
        z = []
        global Counter
        global counter2
        pause_count = 5
        if counter2 < pause_count: #this counter implements a pause diruation in the beginning of the simulation
            counter2 += 1
            Counter = 0
        else:
            Counter +=1


        if Counter >= agents[0].motor1.position_array.__len__()-1:
            Counter = Counter - agents[0].motor1.position_array.__len__()-2

        for agent in agents:
            theta = agent.motor1.position_array[Counter]
            phi = agent.motor2.position_array[Counter] + math.pi

            firstcircle_p = np.array([agent.x_centre_focal,agent.y_centre_focal])

            joint = firstcircle_p + INS_POS_LENGTH1 * np.array([math.cos(theta),math.sin(theta)])\
            + INS_POS_WIDTH2/2*np.array([math.cos(theta+phi),math.sin(theta+phi)])

            point1 = joint + INS_POS_WIDTH2/2*np.array([math.sin(theta+phi),-math.cos(theta+phi)])
            point2 = joint + INS_POS_WIDTH2/2*np.array([-math.sin(theta+phi),math.cos(theta+phi)])
            point3 = point2 + (INS_POS_LENGTH2 - INS_POS_WIDTH2/2)*np.array([math.cos(theta+phi),math.sin(theta+phi)])
            point4 = point1 + (INS_POS_LENGTH2 - INS_POS_WIDTH2/2)*np.array([math.cos(theta+phi),math.sin(theta+phi)])

            joint2 = joint+(INS_POS_LENGTH2 - INS_POS_WIDTH2/2)*np.array([math.cos(theta+phi),math.sin(theta+phi)])

            x.append(np.array([point1,point2,point3,point4]))
            y.append(joint)
            z.append(joint2)

        yield x,y,z

    ani = animation.FuncAnimation(fig, update,data_gen, interval=sim_length*dt)
    plt.show()

