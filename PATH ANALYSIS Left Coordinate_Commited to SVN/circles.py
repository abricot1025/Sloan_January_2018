import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import operator
import parameters as param

# If number_of_layer = 3, number_of_positioner = 37
# If number_of_layer = 4, number_of_positioner = 61
# If number_of_layer = 5, number_of_positioner = 91
# If number_of_layer = 6, number_of_positioner = 127
# If number_of layer = 7, number_of_positioner = 169
# If number_of layer = 8, number_of_positioner = 217
# If number_of layer = 9, number_of_positioner = 271
# If number_of layer = 10, number_of_positioner = 331
# If number_of layer = 11, number_of_positioner = 397
# If number_of layer = 12, number_of_positioner = 469
# If number_of layer = 13, number_of_positioner = 547
# If number_of layer = 14, number_of_positioner = 631
# If number_of layer = 15, number_of_positioner = 721
# If number_of layer = 16, number_of_positioner = 817
# If number_of layer = 17, number_of_positioner = 919
# If number_of layer = 18, number_of_positioner = 1027

position = []
number_of_layer = param.NUMBER_OF_LAYER

R = param.POSITIONER_SEPARATION
r = R

def initialisation(layer):
    position.append((2*R*layer, 0))
    return


def addition_1():
    position.append(tuple((map(operator.add, position[(len(position)-1)], (-R, R * 3 ** 0.5)))))
    return


def addition_2():
    position.append(tuple((map(operator.add, position[(len(position)-1)], (-2*R, 0)))))
    return


def addition_3():
    position.append(tuple((map(operator.add, position[(len(position)-1)], (-R, -R * 3 ** 0.5)))))
    return


def addition_4():
    position.append(tuple((map(operator.add, position[(len(position)-1)], (R, - R * 3 ** 0.5)))))
    return


def addition_5():
    position.append(tuple((map(operator.add, position[(len(position)-1)], (2*R, 0)))))
    return


def addition_6():
    position.append(tuple((map(operator.add, position[(len(position)-1)], (R, R * 3 ** 0.5)))))
    return


def repeat_function(layer, function):                               # Used to repeat each addition, depending on the layer number
    for _ in range(layer):function()
    return


def layer_prod(layer):                                              # Used to produce every layer after the two first ones
    initialisation(layer)
    repeat_function(layer, addition_1)
    repeat_function(layer, addition_2)
    repeat_function(layer, addition_3)
    repeat_function(layer, addition_4)
    repeat_function(layer, addition_5)
    repeat_function((layer - 1), addition_6)
    return

# CIRCLE 0

position.append((0, 0))

# CIRCLE 1

initialisation(1)
addition_1()
addition_2()
addition_3()
addition_4()
addition_5()

# CIRCLE 2

initialisation(2)
repeat_function(2, addition_1)
repeat_function(2, addition_2)
repeat_function(2, addition_3)
repeat_function(2, addition_4)
repeat_function(2, addition_5)
repeat_function(1, addition_6)

# UPPER CIRCLES

for i in range(3, number_of_layer+1):layer_prod(i)

# PROPERTIES

            # Next lines are used to categorise each positioner in empty, fiducial, IR or visible
            # For empty and fiducial, the program simply takes the indexes written in empty.txt and fiducial.txt
empty = open('empty', 'r')
empty = empty.read().split()

for i in range(0, len(empty)):
    empty[i] = int(empty[i])
empty = sorted(empty)

fiducial = open('fiducial', 'r')
fiducial = fiducial.read().split()

for i in range(0, len(fiducial)):
    fiducial[i] = int(fiducial[i])
fiducial = sorted(fiducial)

            # IR's positioner are not well categorised since January 2018's modifications
IR = []

y_IR = np.sqrt(3)*R*2
y_IR_list = [0, y_IR, 2*y_IR, 3*y_IR, 4*y_IR, 5*y_IR, 6*y_IR, 7*y_IR]

for i in range(0, len(y_IR_list)):
    y_IR_list[i] = int(y_IR_list[i])

for i in range(1, len(position)+1):
    if not i in empty:
        if not i in fiducial:
            if abs(int(position[i-1][1])) in y_IR_list:
                IR.append(i)

IR.remove(576)
IR.remove(562)

            # Every positioner not categorised in empty, fiducial or IR is categorised as visible
visible = []

for i in range(1, len(position)+1):
    if not i in empty:
        if not i in fiducial:
            if not i in IR:
                visible.append(i)

# FIGURES
            # Used to visualize the structure and compare it with the desired one prior to running the full program

fig, ax = plt.subplots()

font = {'color':  'black',
        'weight': 'light',
        'size': 8,
        }

circle = []
text = []

for i in range(0, len(position)):
    circle.append(0)
    text.append(0)

for i in range(0, len(position)):
    text[i] = plt.text(position[i][0], position[i][1], i+1, fontdict=font, verticalalignment='center', horizontalalignment='center',size=4)
    ax.add_artist(text[i])

    circle_fiducial = plt.Circle(position[i], radius=r, fill=True, fc='red', edgecolor='red')
    circle_empty = plt.Circle(position[i], radius=r, fill=False, edgecolor='Black')
    circle_IR = plt.Circle(position[i], radius=r, fill=True, fc='blue', edgecolor='blue')
    circle_visible = plt.Circle(position[i], radius=r, fill=True, fc='green', edgecolor='green')

    if i+1 in fiducial:
        ax.add_artist(circle_fiducial)

    if i+1 in visible:
        ax.add_artist(circle_visible)

    #if i+1 in empty:
    #    ax.add_artist(circle_empty)

    if i+1 in IR:
        ax.add_artist(circle_IR)

radius = []
for i in range(0, len(position)):
    radius.append(np.sqrt(position[i][0]**2 + position[i][1]**2))

#limit = R*len(position)/5
focal_plane_radius = param.FOCAL_PLANE_RADIUS
limit = focal_plane_radius*1.4
plt.axis('equal')
plt.axis([-limit, limit, -limit, limit])

circle_main = plt.Circle((0, 0), radius=focal_plane_radius, fill=False, edgecolor='grey', fc=None)
ax.add_artist(circle_main)

rectangle_width = 130
rectangle_height = 170

rectangle = patches.Rectangle((-rectangle_height/2, -rectangle_width/2), rectangle_height, rectangle_width, angle=0.0, fill=False, edgecolor='grey')
#ax.add_artist(rectangle)

ax.set_xlabel('x [mm]')
ax.set_ylabel('y [mm]')

plt.axis('off')

title_image = 'image_' + 'plot_circles' + '.png'
plt.savefig(title_image, format='png', bbox_inches='tight')
print('Plot ready')
plt.close(fig)
