import numpy as np
import os
#import circles as circles
#import parameters as param
import math
import matplotlib.pyplot as plt
import util

r_positioner = 22.4/2
r_focal_plane = 325
alpha_arm = 7.4
beta_arm = 15


def cartesian_to_polar(x, y, xorigin=0.0, yorigin=0.0):
    """

    Helper function to convert Cartesian coordinates to polar coordinates
    (centred at a defined origin). In the polar coordinates, theta is an
    angle measured clockwise from the Y axis.

    :Parameters:

    x: float
        X coordinate of point
    y: float
        Y coordinate of point
    xorigin: float (optional)
        X coordinate of origin (if not zero)
    yorigin: float (optional)
        Y coordinate of origin (if not zero)

    :Returns:

    (r, theta): tuple of 2 floats
        Polar coordinates of point. NOTE: theta is in radians.

    """
    PI2 = 2.0 * math.pi
    PIBY2 = math.pi / 2.0
    xdiff = float(x) - float(xorigin)
    ydiff = float(y) - float(yorigin)
    distsq = (xdiff * xdiff) + (ydiff * ydiff)
    r = math.sqrt(distsq)
    theta = PIBY2 - math.atan2(ydiff, xdiff)
    # Adjust theta to be in the range 0 - 2*PI
    while theta < 0.0:
        theta += PI2
    while theta > PI2:
        theta -= PI2
    return (r, theta)


def polar_to_cartesian(r, theta, r_ref=0.0, theta_ref=0.0):
    """

    Helper function to convert polar coordinates to Cartesian
    coordinates (relative to a defined reference point).

    :Parameters:

    r: float
        Radial distance of point from origin.
    theta: float
        Angular bearing of point, clockwise from Y axis (in radians)
    r_ref: float (optional)
        Radial distance of reference point from origin
        (if reference point is not the origin).
    theta_ref: float (optional)
        Angular bearing of reference point, clockwise from Y axis
        (in radians)(if reference point is not the origin).

    :Returns:

    (x, y): tuple of 2 floats
        Cartesian coordinates of point

    """
    if float(r_ref) > 0.0:
        x = r * math.sin(theta) - r_ref * math.sin(theta_ref)
        y = r * math.cos(theta) - r_ref * math.cos(theta_ref)
    # Old anticlockwise from the X axis code
    #         x = r * math.cos(theta) - r_ref * math.cos(theta_ref)
    #         y = r * math.sin(theta) - r_ref * math.sin(theta_ref)
    else:
        x = r * math.sin(theta)
        y = r * math.cos(theta)
    # Old anticlockwise from the X axis code
    #         x = r * math.cos(theta)
    #         y = r * math.sin(theta)

    return (x, y)

def add_fibre(name_text_file):

    text_file = open(name_text_file, 'r')
    text_file = text_file.readlines()

    output = ["%s \t %s" % (item.strip(), np.random.randint(3)-1) for item in text_file]

    name = os.path.splitext(name_text_file)[0]
    name = str(name) + str('with_fibres.txt')

    f = open(name, "w")
    f.write("\n".join(output))

    open(name_text_file, 'r').close()
    f.close()

    return

def target_position_random(name_text_file):

    text_file = open(name_text_file, 'r')
    text_file = text_file.readlines()

    priority_positioner = 1

    name_global = os.path.splitext(name_text_file)[0]

    for i in range(1,11):
        name = str(name_global) + str('_random_') + str(i) + str('.txt')
        f = open(name, "w")
        fig, ax = plt.subplots()

        for item in text_file:
            r = np.random.rand()*r_focal_plane
            theta = np.random.rand()*360 - 180
            parity = np.random.randint(2)
            output= ["%s \t %s \t %s \t %s \t %s \t %s \t %s \n" % (item.strip(), r, theta, parity, 1-parity, priority_positioner, np.random.randint(3)-1)]
            f.write("\n".join(output))
            x, y = polar_to_cartesian(r, theta)
            plt.plot(x, y, ms=5, marker='.', lw=0, color='red')

        x_positioner, y_positioner = positioner_cartesian(name_text_file)
        plt.plot(x_positioner, y_positioner, ms=5, marker='.', lw=0, color='blue')

        title_image = 'random_' + 'plot_' + str(i) + '.png'
        plt.savefig(title_image, format='png', bbox_inches='tight')
        print('Plot ' + str(i) + ' ready')
        plt.close(fig)
        f.close()

    open(name_text_file, 'r').close()

    return

def target_position_random_normal(name_text_file):

    text_file = open(name_text_file, 'r')
    text_file = text_file.readlines()

    priority_positioner = 1

    name_global = os.path.splitext(name_text_file)[0]

    for i in range(1,11):
        name = str(name_global) + str('_random_normal_') + str(i) + str('.txt')
        f = open(name, "w")
        fig, ax = plt.subplots()

        for item in text_file:
            r = np.random.randn()*r_focal_plane/2
            theta = np.random.rand()*360 - 180
            parity = np.random.randint(2)
            output= ["%s \t %s \t %s \t %s \t %s \t %s \t %s \n" % (item.strip(), r, theta, parity, 1-parity, priority_positioner, np.random.randint(3)-1)]
            f.write("\n".join(output))
            x, y = polar_to_cartesian(r, theta)
            plt.plot(x, y, ms=5, marker='.', lw=0, color='red')

        x_positioner, y_positioner = positioner_cartesian(name_text_file)
        plt.plot(x_positioner, y_positioner, ms=5, marker='.', lw=0, color='blue')

        title_image = 'random_normal_' + 'plot_' + str(i) + '.png'
        plt.savefig(title_image, format='png', bbox_inches='tight')
        print('Plot ' + str(i) + ' ready')
        plt.close(fig)
        f.close()

    open(name_text_file, 'r').close()

    return

def target_position_focus(name_text_file):

    text_file = open(name_text_file, 'r')
    text_file = text_file.readlines()

    priority_positioner = 1

    name_global = os.path.splitext(name_text_file)[0]

    for i in range(1,11):
        r_focus = r_focal_plane * np.random.rand()
        theta_focus = np.random.rand() * 360 - 180
        x_focus, y_focus = polar_to_cartesian(r_focus, math.radians(theta_focus))
        name = str(name_global) + str('_focus_') + str(i) + str('.txt')
        f = open(name, "w")

        fig, ax = plt.subplots()

        for item in text_file:
            r_target_random = np.random.rand()*50
            theta_target_random = np.random.rand() * 360 - 180
            x_target_random, y_target_random = polar_to_cartesian(r_target_random, math.radians(theta_target_random))
            x = x_target_random + x_focus
            y = y_target_random + y_focus
            r, theta = cartesian_to_polar(x, y)
            theta = math.degrees(theta)
            parity = np.random.randint(2)
            output= ["%s \t %s \t %s \t %s \t %s \t %s \t %s \n" % (item.strip(), r, theta, parity, 1-parity, priority_positioner, np.random.randint(3)-1)]
            f.write("\n".join(output))
            x, y = polar_to_cartesian(r, math.radians(theta))

            plt.plot(x, y, ms=5, marker='.', lw=0, color='red')

        x_positioner, y_positioner = positioner_cartesian(name_text_file)
        plt.plot(x_positioner, y_positioner, ms=5, marker='.', lw=0, color='blue')

        title_image = 'focus_' + 'plot_' + str(i) + '.png'
        plt.savefig(title_image, format='png', bbox_inches='tight')
        print('Plot ' + str(i) + ' ready')
        plt.close(fig)

        f.close()

    open(name_text_file, 'r').close()

    return

def target_position_focus_normal(name_text_file):

    text_file = open(name_text_file, 'r')
    text_file = text_file.readlines()

    priority_positioner = 1

    name_global = os.path.splitext(name_text_file)[0]

    for i in range(1, 11):
        r_focus = r_focal_plane * np.random.randn()
        theta_focus = np.random.rand() * 360 - 180
        x_focus, y_focus = polar_to_cartesian(r_focus, math.radians(theta_focus))
        name = str(name_global) + str('_focus_normal_') + str(i) + str('.txt')
        f = open(name, "w")

        fig, ax = plt.subplots()

        for item in text_file:
            r_target_random = np.random.randn() * 50
            theta_target_random = np.random.rand() * 360 - 180
            x_target_random, y_target_random = polar_to_cartesian(r_target_random, math.radians(theta_target_random))
            x = x_target_random + x_focus
            y = y_target_random + y_focus
            r, theta = cartesian_to_polar(x, y)
            theta = math.degrees(theta)
            parity = np.random.randint(2)
            output = ["%s \t %s \t %s \t %s \t %s \t %s \t %s \n" % (item.strip(), r, theta, parity, 1-parity, priority_positioner,
            np.random.randint(3) - 1)]
            f.write("\n".join(output))
            x, y = polar_to_cartesian(r, math.radians(theta))

            plt.plot(x, y, ms=5, marker='.', lw=0, color='red')

        x_positioner, y_positioner = positioner_cartesian(name_text_file)
        plt.plot(x_positioner, y_positioner, ms=5, marker='.', lw=0, color='blue')

        title_image = 'focus_normal_' + 'plot_' + str(i) + '.png'
        plt.savefig(title_image, format='png', bbox_inches='tight')
        print('Plot ' + str(i) + ' ready')
        plt.close(fig)

        f.close()

    open(name_text_file, 'r').close()

    return

def target_position_homogeneous(name_text_file):

    text_file = open(name_text_file, 'r')
    text_file_readlines = text_file.readlines()
    text_lines = text_file.read().split()
    length_file = len(text_lines)

    #r_target_random = np.random.rand()*r_focal_plane
    #theta_target_random = np.random.rand()*360 - 180
    #x_target_random =
    #arg_parity_1 = np.random.randint(2)
    #arg_parity_2 = np.random.randint(2)
    priority_positioner = 1
    #fibre = np.random.randint(3)-1

    name_global = os.path.splitext(name_text_file)[0]

    for i in range(1,11):
        r_focus = r_focal_plane * np.random.rand()
        theta_focus = np.random.rand() * 360 - 180
        # Output : r_target, theta_target, arg_parity_1, arg_parity_2, priority_positioner, fibre
        name = str(name_global) + str('_homogeneous_') + str(i) + str('.txt')
        f = open(name, "w")
        j = 0

        fig, ax = plt.subplots()

        for item in text_file_readlines:
            r_positioner = eval(item.split()[0])
            theta_positioner = eval(item.split()[1])
            x_positioner = r_positioner * np.cos(math.radians(theta_positioner))
            y_positioner = r_positioner * np.sin(math.radians(theta_positioner))
            r_target_random = alpha_arm + np.random.rand() * beta_arm
            theta_target_random = np.random.rand() * 360 - 180
            x_target_random = r_target_random * np.cos(math.radians(theta_target_random))
            y_target_random = r_target_random * np.sin(math.radians(theta_target_random))
            x_target = x_target_random + x_positioner
            y_target = y_target_random + y_positioner
            r_target = np.sqrt(x_target ** 2 + y_target ** 2)
            theta_target = math.degrees(math.atan2(y_target,x_target))
            parity = np.random.randint(2)
            output= ["%s \t %s \t %s \t %s \t %s \t %s \t %s \n" % (item.strip(), r_target, theta_target, parity, 1-parity, priority_positioner, np.random.randint(3)-1)]
            f.write("\n".join(output))
            j = j + 1
            x, y = polar_to_cartesian(r_target, theta_target)
            plt.plot(x, y, ms=5, marker='.', lw=0, color='red')

        x_positioner, y_positioner = positioner_cartesian(name_text_file)
        plt.plot(x_positioner, y_positioner, ms=5, marker='.', lw=0, color='blue')

        title_image = 'homogeneous_' + 'plot_' + str(i) + '.png'
        plt.savefig(title_image, format='png', bbox_inches='tight')
        print('Plot ' + str(i) + ' ready')
        plt.close(fig)

        f.close()

    open(name_text_file, 'r').close()

    return

def positioner_cartesian(name_text_file):

    text_file = open(name_text_file, 'r')
    text_file_readlines = text_file.readlines()
    x_positioner =[]
    y_positioner = []
    for item in text_file_readlines:
        r_positioner = eval(item.split()[0])
        theta_positioner = eval(item.split()[1])
        x_positioner.append(r_positioner * np.cos(math.radians(theta_positioner)))
        y_positioner.append(r_positioner * np.sin(math.radians(theta_positioner)))

    return x_positioner, y_positioner

target_position_random('650output_generator.txt')
#target_position_random_normal('650output_generator.txt')
#target_position_focus('650output_generator.txt')
#target_position_focus_normal('650output_generator.txt')
#target_position_homogeneous('650output_generator.txt')
