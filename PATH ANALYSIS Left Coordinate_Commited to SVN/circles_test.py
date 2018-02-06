import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import operator
import parameters as param


position = []
number_of_circle = 650


#r = (param.INS_POS_LENGTH1 + param.INS_POS_LENGTH2)/2
#R = 1.0*r
R = 11
r = R

def initialisation(num_circle):
    position.append((2*R*num_circle, 0))
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


def repeat_function(num_circle, function):
    for _ in range(num_circle):function()
    return


def circle_prod(num_circle):
    initialisation(num_circle)
    repeat_function(num_circle, addition_1)
    repeat_function(num_circle, addition_2)
    repeat_function(num_circle, addition_3)
    repeat_function(num_circle, addition_4)
    repeat_function(num_circle, addition_5)
    repeat_function((num_circle - 1), addition_6)
    #print('Au cercle ' + str(num_circle) + ', il y a ' + str(len(position)) + ' disques')
    return

# CIRCLE 0

position.append((0, 0))
#print('Au cercle 0, il y a 1 cercle')

# CIRCLE 1

initialisation(1)
addition_1()
addition_2()
addition_3()
addition_4()
addition_5()

#print('Au cercle ' + str(1) + ', il y a ' + str(len(position)) + ' disques')

# CIRCLE 2

initialisation(2)
repeat_function(2, addition_1)
repeat_function(2, addition_2)
repeat_function(2, addition_3)
repeat_function(2, addition_4)
repeat_function(2, addition_5)
repeat_function(1, addition_6)

#print('Au cercle ' + str(2) + ', il y a ' + str(len(position)) + ' disques')

# CIRCLE ++

for i in range(3, 400):
    circle_prod(i)
    if len(position) > number_of_circle:
        break
    else:
        continue


'''
circle_prod(3) # 37
circle_prod(4) # 61
circle_prod(5) # 91
circle_prod(6) # 127
circle_prod(7) # 169
circle_prod(8) # 217
circle_prod(9) # 271
circle_prod(10) # 331
circle_prod(11) # 397
circle_prod(12) # 469
circle_prod(13) # 547
circle_prod(14) # 631
circle_prod(15) # 721
circle_prod(16) # 817
circle_prod(17) # 919
circle_prod(18) # 1027

'''

# TARGETS
r_target = open('r_target', 'r')
r_target = r_target.read().split()

for i in range(0, len(r_target)):
    r_target[i] = float(r_target[i])

#print('r_target')
#print(r_target)

theta_target = open('theta_target', 'r')
theta_target = theta_target.read().split()

for i in range(0, len(theta_target)):
    theta_target[i] = float(theta_target[i])

#print('theta_target')
#print(theta_target)


x_target = []
y_target = []

for i in range(0, len(theta_target)):
    x_target.append(r_target[i] * np.cos(theta_target[i]))
    y_target.append(r_target[i] * np.sin(theta_target[i]))

# PROPERTIES

#empty = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 39, 40, 41, 47, 48, 49, 51, 52, 53, 59, 60, 61, 65, 74, 80, 89]

#empty = [6, 7, 8]
empty = []


#fiducial = open('fiducial', 'r')
#fiducial = fiducial.read().split()

#for i in range(0, len(fiducial)):
#   fiducial[i] = int(fiducial[i])
#fiducial = sorted(fiducial)

#fiducial = [28, 29, 30]
fiducial = []

IR = [6, 7, 8]

'''
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
'''


visible = []

for i in range(1, len(position)+1):
    if not i in empty:
        if not i in fiducial:
            if not i in IR:
                visible.append(i)

# FIGURES

'''
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
    circle_IR = plt.Circle(position[i], radius=r, fill=True, fc='blue', edgecolor='blue')
    circle_visible = plt.Circle(position[i], radius=r, fill=True, fc='green', edgecolor='green')

    if i+1 in fiducial:
        ax.add_artist(circle_fiducial)

    if i+1 in visible:
        ax.add_artist(circle_visible)

    if i+1 in IR:
        ax.add_artist(circle_IR)

#plt.plot(x_target, y_target, lw=0, ls='-', marker='o', color='red')

radius = []
for i in range(0, len(position)):
    radius.append(np.sqrt(position[i][0]**2 + position[i][1]**2))

#limit = R*len(position)/5
focal_plane_radius = 325
limit = focal_plane_radius*1.2
#axes = plt.gca()
#axes.set_xlim([-limit, limit])
#axes.set_ylim([-limit, limit])
plt.axis('equal')
plt.axis([-limit*1.2, limit, -limit, limit])

circle_main = plt.Circle((0, 0), radius=focal_plane_radius, fill=False, edgecolor='grey', fc=None)
ax.add_artist(circle_main)

rectangle_width = 130
rectangle_height = 170

rectangle = patches.Rectangle((-rectangle_height/2, -rectangle_width/2), rectangle_height, rectangle_width, angle=0.0, fill=False, edgecolor='grey')
ax.add_artist(rectangle)

#for i in range(len(position)):
    #print(str(i+1)+ ' ' + str(position[i]))

title_image = 'image_' + 'plot_circles' + '.png'
plt.savefig(title_image, format='png', bbox_inches='tight')
print('Plot ready')
plt.close(fig)


#print(position)
'''
