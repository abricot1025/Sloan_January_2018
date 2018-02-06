import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots()

R = 1
pi = np.pi

r_3 = 3 ** -0.5

def polar_to_cartesian((r, theta)):
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return (x, y)

position_polar_1 = [
    (0, 0)
]

d_2 = 2*R
position_polar_2 = [
    (d_2, pi),
    (d_2, 2*pi/3),
    (d_2, pi/3),
    (d_2, 0),
    (d_2, -pi/3),
    (d_2, -2*pi/3)
]

d_3a = R * 4
d_3b = R * 12**0.5
position_polar_3 = [
    (d_3a, 0),
    (d_3b, pi/6),
    (d_3a, 2*pi/6),
    (d_3b, 3*pi/6),
    (d_3a, 4*pi/6),
    (d_3b, 5*pi/6),
    (d_3a, 6*pi/6),
    (d_3b, 7*pi/6),
    (d_3a, 8*pi/6),
    (d_3b, 9*pi/6),
    (d_3a, 10*pi/6),
    (d_3b, 11*pi/6)
]

d_4a = 6*R
d_4b = R * (3.5**2 + 4**2)**0.5
position_polar_4 = [
    (d_4a, 0),
    (d_4b, np.arctan(2*np.sin(pi/3)/5)),
    (d_4b, np.arctan(d_3b/d_3a)),
    (d_4a, 3*pi/9),
    (d_4b, np.arctan(3**1.5)),
    (d_4b, pi - np.arctan(3**1.5)),
    (d_4a, 6*pi/9),
    (d_4b, pi - np.arctan(d_3b/d_3a)),
    (d_4b, pi - np.arctan(2*np.sin(pi/3)/5)),
    (d_4a, 9*pi/9),
    (d_4b, pi + np.arctan(2*np.sin(pi/3)/5)),
    (d_4b, pi + np.arctan(d_3b/d_3a)),
    (d_4a, 12*pi/9),
    (d_4b, - np.arctan(3**1.5)),
    (d_4b, pi + np.arctan(3**1.5)),
    (d_4a, 15*pi/9),
    (d_4b, -np.arctan(d_3b/d_3a)),
    (d_4b, -np.arctan(2*np.sin(pi/3)/5))
        ]

d_5a = 8 * R
d_5b = R * 2 * 12 ** 0.5
d_5c = R * 2 * 13 ** 0.5
position_polar_5 = [
    (d_5a, 0),
    (d_5c, np.arctan(2 * 3 ** 0.5) - (pi/3)),
    (d_5b, (0*pi/3) + (pi/6)),
    (d_5c, pi + np.arctan(-2 * 3 ** 0.5) - (pi/3)),
    (d_5a, pi/3),
    (d_5c, np.arctan(2 * 3 ** 0.5)),
    (d_5b, (1*pi/3) + (pi/6)),
    (d_5c, pi + np.arctan(-2 * 3 ** 0.5)),
    (d_5a, 2*pi/3),
    (d_5c, np.arctan(2 * 3 ** 0.5) + (pi/3)),
    (d_5b, (2*pi/3) + (pi/6)),
    (d_5c, pi + np.arctan(-2 * 3 ** 0.5) + (pi/3)),
    (d_5a, pi),
    (d_5c, np.arctan(2 * 3 ** 0.5) + (2*pi/3)),
    (d_5b, (3*pi/3) + (pi/6)),
    (d_5c, pi + np.arctan(-2 * 3 ** 0.5) + (2*pi/3)),
    (d_5a, 4*pi/3),
    (d_5c, np.arctan(2 * 3 ** 0.5) + (3*pi/3)),
    (d_5b, (4*pi/3) + (pi/6)),
    (d_5c, pi + np.arctan(-2 * 3 ** 0.5) + (3*pi/3)),
    (d_5a, 5*pi/3),
    (d_5c, np.arctan(2 * 3 ** 0.5) + (4*pi/3)),
    (d_5b, (5*pi/3) + (pi/6)),
    (d_5c, pi + np.arctan(-2 * 3 ** 0.5) + (4*pi/3)),
]

position_polar = []

for i in range(0, len(position_polar_1)):
    position_polar.append(position_polar_1[i])

for i in range(0, len(position_polar_2)):
    position_polar.append(position_polar_2[i])

for i in range(0, len(position_polar_3)):
        position_polar.append(position_polar_3[i])

for i in range(0, len(position_polar_4)):
    position_polar.append(position_polar_4[i])

for i in range(0, len(position_polar_5)):
    position_polar.append(position_polar_5[i])

circle = []
positions = []
for i in range(0, len(position_polar)):
    circle.append(0)
    positions.append(0)

for i in range(0, len(position_polar)):
    positions[i] = polar_to_cartesian(position_polar[i])
    circle[i] = plt.Circle(positions[i], radius=R, fill=True, edgecolor=[0, 0.6, 1], facecolor=[1, 1, 1])
    ax.add_artist(circle[i])


limit = 12
axes = plt.gca()
axes.set_xlim([-limit, limit])
axes.set_ylim([-limit, limit])

print('Il y a ' + str(len(position_polar)) + ' points/')

plt.show()