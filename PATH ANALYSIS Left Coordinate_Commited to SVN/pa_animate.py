try:
    import mocpath.conflict.fps_shared as fps
except ImportError:
    import fps_shared as fps

try:
    import mocpath.util as util
except ImportError:
    import util as util

import path_generator as pg

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")
import matplotlib.animation as animation
import math
import circles
import parameters as param

INS_POS_LENGTH1 = fps.INS_POS_LENGTH1     # Length of alpha arm in mm.
INS_POS_WIDTH1 = fps.INS_POS_WIDTH1       # Width of alpha arm in mm.
INS_POS_LENGTH2 = fps.INS_POS_LENGTH2       # Length of beta arm in mm (including metrology zone)
INS_POS_WIDTH2 = fps.INS_POS_WIDTH2         # Width of beta arm in mm (before metrology zone).
PITCH = INS_POS_LENGTH1 + INS_POS_LENGTH2


# Counters for making replay and delays in the simulation for recording
global Counter
Counter = -1

global counter2
counter2 = 0

def run_animation(agents, target_file):

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

    for agent in agents:
        # Draw the pich and the alpha motor and the target point

        color_alpha_pitch=''
        color_alpha_motor=''
        fill = False
        if agent.ident in circles.empty:                # if agent (positioner) is empty
            color_alpha_pitch = [0.8,0.8,0.8,0]         # pitch and motor are transparent
            color_alpha_motor = [0.4,0.4,0.4,0]

        elif agent.ident in circles.fiducial:           # if agent (positioner) is fiducial
            color_alpha_pitch = [0.8,0.8,0.8,0]         # pitch is transparent and 'motor' is red and filled
            color_alpha_motor = 'red'
            fill = True

        else:                                           # Otherwise, same color as MOONS' ones
            color_alpha_pitch = [0.8,0.8,0.8]
            color_alpha_motor = [0.4,0.4,0.4]

        firstcircle_p = np.array([agent.x_centre_focal,agent.y_centre_focal])
        pitchcircle = plt.Circle(firstcircle_p, radius= PITCH/2, color=color_alpha_pitch, fill=False)
        motorcircle = plt.Circle(firstcircle_p, radius= INS_POS_LENGTH1, color=color_alpha_motor, fill=fill)

        # make the circle for the targets
        thetag = agent.motor1.position_array[-1]
        phig = agent.motor2.position_array[-1]

        # make the circle for the targets from the target list

        #targetcircle = plt.Circle(util.polar_to_cartesian(agent.r_fibre_focal,agent.theta_fibre_focal), radius= 1, color='r', fill=True)
        #subfig.add_patch(targetcircle)

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

        color_beta = ''
        if agent.ident in circles.empty:                # if agent (positioner) is empty
            color_beta = [0,0,0,0]                      # Beta arm is transparent
        if agent.ident in circles.fiducial:             # if agent (positioner) is fiducial
            color_beta = [0,0,0,0]                      # Beta arm is transparent
        if agent.ident in circles.IR:                   # if agent (positioner) has an IR fiber
            color_beta = 'blue'                         # Beta arm is blue
        if agent.ident in circles.visible:              # if agent (positioner) has no IR fiber
            color_beta = 'green'                        # Beta arm is green

        #joint_polygon = plt.Polygon([point1,point2,point3,point4], edgecolor=[0,0.6,1], facecolor=[0,0.6,1])

        joint_polygon = plt.Polygon([point1,point2,point3,point4], edgecolor=color_beta, facecolor=color_beta)
        subfig.add_patch(joint_polygon)

        joint_circle1 = plt.Circle(joint,radius=INS_POS_WIDTH2/2,fill=True,edgecolor=color_beta, facecolor=color_beta)
        joint_circle2 = plt.Circle(joint+(INS_POS_LENGTH2 - INS_POS_WIDTH2/2)*np.array([math.cos(theta+phi),math.sin(theta+phi)])\
                                   ,radius=INS_POS_WIDTH2/2,fill=True,edgecolor=color_beta, facecolor=color_beta)

        color_fiber = ''
        if agent.ident in circles.empty:                # if agent (positioner) is empty
            color_fiber = [0.8,0.8,0.8,0]               # fiber is transparent
        elif agent.ident in circles.fiducial:           # if agent (positioner) is fiducial
            color_fiber = [0.8,0.8,0.8,0]               # fiber is transparent
        else:
            color_fiber=[0,1,1]                         # Otherwise, same color as MOONS' ones
        fiber_circle = plt.Circle(joint+(INS_POS_LENGTH2 - INS_POS_WIDTH2/2)*np.array([math.cos(theta+phi),math.sin(theta+phi)])\
                                   ,radius=INS_POS_WIDTH2/4,fill=True,edgecolor=color_fiber, facecolor=color_fiber)
        subfig.add_patch(joint_circle1)
        subfig.add_patch(joint_circle2)
        subfig.add_patch(fiber_circle)


        beta_polygon.append(joint_polygon)
        beta_circle1.append(joint_circle1)
        beta_circle2.append(joint_circle2)
        fiber.append(fiber_circle)

    fr2 = open(target_file, 'r')

    for i, line in enumerate(fr2):
        x1, x2, x3, x4, x5, x6, x7, x8 = map(float, line.split())
        color_target = ''                               # Target color depends on the fiber necessary
        if x8 == 1: color_target = 'orange'             # Orange for IR
        elif x8 == 0: color_target = 'red'              # Red for Visible
        elif x8 == -1: color_target = 'pink'            # Pink for Calibration

        targetcircle = plt.Circle(util.polar_to_cartesian(x3,x4), radius=1, color=color_target, fill=True)
        subfig.add_patch(targetcircle)

    rectangle_width = 130
    rectangle_height = 170

    rectangle = plt.Polygon([(rectangle_height/2, -rectangle_width/2 ), (rectangle_height/2, rectangle_width/2), (-rectangle_height/2, rectangle_width/2), (-rectangle_height/2, -rectangle_width/2)], edgecolor='grey', fill=False,)

    subfig.add_patch(rectangle)

    circle_main = plt.Circle((0, 0), radius=param.FOCAL_PLANE_RADIUS, fill=False, edgecolor='grey', fc=None)
    subfig.add_patch(circle_main)

    def update(data):
        x,y,z = data
        for i in range(len(x)):
            beta_polygon[i].set_xy(x[i])
            beta_circle1[i].center = y[i]
            beta_circle2[i].center = z[i]
            fiber[i].center = z[i]


    def data_gen():
        x = []
        y = []
        z = []
        global Counter
        global counter2
        pause_count = 5
        if counter2 < pause_count: #this counter implements a pause diruation in the begining of the simulation
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

    ani = animation.FuncAnimation(fig, update, data_gen(), interval=param.sim_length*param.dt)

    ani.save('ani.html')
    plt.show()



