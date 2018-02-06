"""

MOONS Fibre Positioning System Utility Module

Contains a collection of shared functions and utilities.

24 Jul 2014: Created using functions extracted from fps_classes.py
06 Sep 2014: Added generator function for grid coordinates.
             Allow positioner grids to define the origin of their
             coordinate system at the centre.
20 Mar 2015: Polar coordinate system definition changed so that theta
             is measured clockwise from the focal plane Y axis
             (the MOONS convention).
14 Jan 2016: Additional tests. Added elbow_parity_str, polygon_to_points,
             line_intersects_circle and triangle_intersects_circle
             functions.
07 Jun 2016: Temporarily restored old coordinate conversion functions
             until the path analysis code can be brought up to date.
21 Jul 2016: Corrected a typo in a log message.

@author: Steven Beard (UKATC)

"""

from __future__ import division
import sys, math
import warnings

import numpy as np

# Define the limits for floating point values.
EPS = 10 * sys.float_info.epsilon  # Smallest possible increment (x 10)
MAXFLOAT = sys.float_info.max      # Largest possible value.

# Calculate important constants once.
ROOT3 = math.sqrt(3.0)
ROOT3BY2 = ROOT3 / 2.0
PI2 = 2.0 * math.pi
PIBY2 = math.pi / 2.0
COS30 = math.cos(math.radians(30.0))

PARITY_RIGHT = 1       # Constant meaning right-armed parity
PARITY_LEFT = -1       # Constant meaning left-armed parity

#
# A collection of global helper functions
#
def elbow_parity_str( parity ):
    """
    
    Convert the numerical code for elbow parity into a descriptive string.
    
    :Parameters:
    
    parity: int
        The elbow parity:
            
        * 1 means elbow right armed
        * -1 means elbow left armed

    """
    if parity == PARITY_RIGHT:
        return "RIGHT-armed"
    else:
        return "LEFT-armed"

# FIXME: DELETE THIS FUNCTION.
def cartesian_to_polar_OLD(x, y, xorigin=0.0, yorigin=0.0):
    """

    Helper function to convert Cartesian coordinates to polar coordinates
    (centred at a defined origin).
    Angle measured anti-clockwise from X axis.

    OLD OBSOLETE VERSION USED BY PATH ANALYSIS - 13 JUN 2016
    DELETE WHEN COMPATIBILITY PROBLEMS SORTED OUT

    """
    #FIXME: OBSOLETE COORDINATE SYSTEM
    xdiff = float(x) - float(xorigin)
    ydiff = float(y) - float(yorigin)
    distsq = (xdiff * xdiff) + (ydiff * ydiff)
    r = math.sqrt(distsq)
    theta = math.atan2(ydiff, xdiff)
    return (r, theta)

# FIXME: DELETE THIS FUNCTION.
def polar_to_cartesian_OLD(r, theta, r_ref=0.0, theta_ref=0.0):
    """

    Helper function to convert polar coordinates to Cartesian
    coordinates (relative to a defined reference point).
    Angle measured anti-clockwise from X axis.

    OLD OBSOLETE VERSION USED BY PATH ANALYSIS - 13 JUN 2016
    DELETE WHEN COMPATIBILITY PROBLEMS SORTED OUT

    """
    #FIXME: OBSOLETE COORDINATE SYSTEM
    if float(r_ref) > 0.0:
        x = r * math.cos(theta) - r_ref * math.cos(theta_ref)
        y = r * math.sin(theta) - r_ref * math.sin(theta_ref)
    else:
        x = r * math.cos(theta)
        y = r * math.sin(theta)

    return (x, y)

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

def hexagonal_to_cartesian(column, row, pitch, xzero=0.0, yzero=0.0):
    """
        
    Calculate the X, Y location of a given column and row
    within a hexagonal grid.
        
    :Parameters:
        
    column: int
        Column number in hexgonal grid.
    row: int
        Row number in hexagonal grid.
    pitch: float
        Distance between neighbouring grid points (in mm).
    xzero: float (optional)
        X coordinate of zero point. Default 0.0.
    yzero: float (optional)
        Y coordinate of zero point. Default 0.0.
        
    :Returns:
    
    (x, y): tuple of 2 floats
        Cartesian coordinates of centre of grid point.
        
    """
    # The location of the centre of the grid point depends on
    # its position in the hexagonal grid.
    if row % 2 == 0:
        # An even row - columns are aligned with the origin
        xpos = (column * pitch) - xzero
    else:
        # An odd row - columns are offset by pitch/2 from the origin
        xpos = ((column + 0.5) * pitch) - xzero
    ypos = (row * ROOT3BY2 * pitch) - yzero
    return (xpos, ypos)

def cartesian_to_hexagonal(x, y, pitch, xzero=0.0, yzero=0.0):
    """
        
    Calculate the column and row within a hexagonal grid of a
    given X, Y location.
        
    :Parameters:
        
    x: int
        Column number in hexgonal grid.
    y: int
        Row number in hexagonal grid.
    pitch: float
        Distance between neighbouring grid points (in mm).
    xzero: float (optional)
        Zero point offset for X coordinate. Default 0.0.
    yzero: float (optional)
        Zero point offset for Y coordinate. Default 0.0.
        
    :Returns:
    
    (column, row): tuple of 2 ints
        Location within hexagonal grid.
        
    """
    # First determine the nearest row number.
    row = int(round( ((yzero + y) / (ROOT3BY2 * pitch)) ))
        
    # The location of the centre of the positioner depends on
    # its position in the hexagonal grid.
    if row % 2 == 0:
        # An even row - columns are aligned with the origin
        column = int(round( ((xzero + x) / pitch) ))
    else:
        # An odd row - columns are offset by pitch/2 from the origin
        column = int(round( ((xzero + x) / pitch) - 0.5 ))
    return (column, row)

def flat_to_curved_r(flat_r, curvrad):
    """
    
    Convert a radial distance from the centre of a circular flat plane
    into the radial distance projected onto a curved focal plane of
    given radius of curvature.
    
    :Parameters:
    
    flat_r: float
        The radial distance from the centre of the flat plane.
    curvrad: float
        The radius of curvature of the curved plane.
        
    :Returns:
    
    curved_r: float
        The radial distance projected onto the curved plane.
        Note: curved_r approaches flat_r when curvrad tends to infinity.
        None is returned when flat_r > curvrad
    
    """
    if curvrad is not None:
        ratio = abs(flat_r / curvrad)
        if ratio <= 1.0:
            curved_r = curvrad * math.asin(ratio)
            return curved_r
        else:
            # Solution is not valid when flat_r > curvrad
            return None
    else:
        # If the radius of curvature is infinite or not a number
        # the flat and curved radial distances are the same
        return flat_r

def rotate_coordinates(x, y, theta):
    """
    
    Rotate a pair of Cartesian coordinates (x,y) onto a new
    coordinate grid which is rotated by angle theta.

    :Parameters:

    x: float
        X coordinate of point in the original grid
    y: float
        Y coordinate of point in the original grid
    theta: float
        Rotation angle between first grid and second grid (in radians)

    :Returns:
    
    (xnew, ynew): tuple of 2 floats
        The new coordinates
    
    """
    xnew = x * math.cos(theta) + y * math.sin(theta)
    ynew = x * math.sin(theta) + y * math.cos(theta)
    return (xnew, ynew)

def distance_squared(x1, y1, x2, y2):
    """
    
    Helper function to return the square of the distance between
    points (x1,y1) and (x2,y2).

    :Parameters:

    x1: float
        X coordinate of first point
    y1: float
        Y coordinate of first point
    x2: float
        X coordinate of second point
    y2: float
        Y coordinate of second point

    :Returns:
    
    distsq: float
        The square of the distance, (x1-x2)**2 + (y1-y2)**2
    
    """
    xdiff = float(x1) - float(x2)
    ydiff = float(y1) - float(y2)
    distsq = (xdiff * xdiff) + (ydiff * ydiff)
    return distsq

def closer_than(x1, y1, x2, y2, limitsq):
    """
        
    Helper function to determine whether the square of the distance
    between points (x1,y1) and (x2,y2) is less than or equal to limitsq.
    
    The squares of the distances are compared to save computation time.

    :Parameters:

    x1: float
        X coordinate of first point
    y1: float
        Y coordinate of first point
    x2: float
        X coordinate of second point
    y2: float
        Y coordinate of second point
    limitsq: float
        The square of the limiting distance.

    :Returns:

    True or False
        
    """
    distsq = distance_squared(x1, y1, x2, y2)
    if distsq <= limitsq:
        return True
    else:
        return False

def arms_intersect(line1, line2, half_width, maxangle=None, thinfactor=100.0,
                   quick=True):
    """
        
    Helper function to determine whether two arm
    segments intersect.
    
    NOTE: This function decides whether to call the lines_intersect
    function or the quadrangles_intersect depending on the "quick"
    parameter and the arm width given.
        
    :Parameters:

    line1: tuple of (xstart, ystart, xend, yend)
        First line segment
    line2: tuple of (xstart, ystart, xend, yend)
        Second line segment
    half_width: float
        Half the arm width in the same units as the line coordinates.
        This is used to define the tolerance for a simple line intersection
        check, or to define the corner coordinates for a full quadrangle
        intersection check.
    maxangle: float (optional)
        The maximum angle (in degrees) at which a conflict is possible.
        Two arms only conflict if the relative angle between their
        orientations is less than this limit.
        The default value of None means that all angles can conflict.
    thinfactor: float (optional)
        The definition of a thin line. A line is thin if its length
        is more than thinfactor times its width (default 100).
        Thin lines are always tested using the "lines_intersect"
        function. Thick lines are tested using the much slower
        "quadrangles_intersect" function unless the "quick" parameter
        (below) is True.
    quick: bool (optional)
        Set True for a quick but simple intersection test
        where arms are always considered lines with a tolerance.
        Set False for an accurate but very slow intersection test
        where all thick arms are treated as 2-D quadrangles.

    :Returns:

    True or False
        
    """
    ydiff1 = line1[3] - line1[1]
    xdiff1 = line1[2] - line1[0]
    length1sq = xdiff1 * xdiff1 + ydiff1 * ydiff1
    ydiff2 = line2[3] - line2[1]
    xdiff2 = line2[2] - line2[0]
    length2sq = xdiff2 * xdiff2 + ydiff2 * ydiff2
    
    # Check for relative orientations greater than the conflict limit.
    if maxangle is not None:
        angle1 = math.atan2(ydiff1, xdiff1)
        angle2 = math.atan2(ydiff2, xdiff2)
        if abs(angle2 - angle1) > math.radians(maxangle):
            return False

    thinfactsq = thinfactor * thinfactor
    maxlengthsq = max(length1sq, length2sq)
    half_widthsq = half_width * half_width
    if quick or ((half_widthsq * thinfactsq) < maxlengthsq):
        # Treat thin lines as lines with a tolerance.
        # This will tend to round off the ends of the line.
        return lines_intersect(line1, line2, tolerance=half_width)
    else:
        # Treat thick lines as quadrangles. Useful for lines with sharp
        # boundaries, but a much slower option.
        angle1 = math.atan2(ydiff1, xdiff1)
        xdelta1 = half_width * math.sin(angle1)
        ydelta1 = half_width * math.cos(angle1)
        quad1 = [(line1[0]-xdelta1,line1[1]+ydelta1),
                 (line1[0]+xdelta1,line1[1]-ydelta1),
                 (line1[2]-xdelta1,line1[3]+ydelta1),
                 (line1[2]+xdelta1,line1[3]-ydelta1)]
        angle2 = math.atan2(ydiff2, xdiff2)
        xdelta2 = half_width * math.sin(angle2)
        ydelta2 = half_width * math.cos(angle2)
        quad2 = [(line2[0]-xdelta2,line2[1]+ydelta2),
                 (line2[0]+xdelta2,line2[1]-ydelta2),
                 (line2[2]-xdelta2,line2[3]+ydelta2),
                 (line2[2]+xdelta2,line2[3]-ydelta2)]
        return quadrangles_intersect(quad1, quad2)

def lines_intersect(line1, line2, tolerance=0.0):
    """
        
    Helper function to determine whether two line segments intersect.
    
    This function is most suitable for thin lines. Intersection is
    tested by surrounding each point with a tolerance of the given
    radius. Thick lines requiring a large tolerance are effectively
    given rounded ends, which causes a slight overestimate in the
    number of intersections. A more accurate, but much slower, way
    of testing thick lines is to use the "quadrangles_intersect"
    function.
        
    :Parameters:

    line1: tuple of (xstart, ystart, xend, yend)
        First line segment
    line2: tuple of (xstart, ystart, xend, yend)
        Second line segment
    tolerance: float (optional)
        The tolerance of the intersection, to account for floating
        point rounding errors and a finite width of the lines.
        The tolerance should be roughly half the line width.
        Defaults to the smallest possible floating point value.
        NOTE: The tolerance should be much less than the
        line length - otherwise use the quadrangles_intersect
        function.

    :Returns:

    True or False
        
    """    
    assert isinstance(line1, (list,tuple))
    assert isinstance(line1, (list,tuple))
    assert len(line1) > 3
    assert len(line2) > 3
#     print("Lines intersect? Line 1:", str(line1), "Line 2:", str(line2))
    tolerance = float(tolerance) + EPS
#     print("Tolerance = ", tolerance)
    xlength1 = line1[2]-line1[0]
    ylength1 = line1[3]-line1[1]
    xlength2 = line2[2]-line2[0]
    ylength2 = line2[3]-line2[1]
    
    if line1 == line2:
        # Identical lines intersect by definition
#         print("Identical lines")
        return True
        
    if abs(xlength1) < EPS and abs(xlength2) < EPS:
        # Both slopes are infinite. The lines intersect only
        # if they have exactly the same X coordinate, to within
        # the tolerance.
        if abs(line1[0] - line2[0]) < tolerance:
#             print("Infinite slope and same X")
            return True
        else:
            return False
        
    elif abs(xlength1) < EPS:
        # Infinite line1 slope special case.
        # The lines intersect where the X coordinate of line1 is
        # substituted into the equation for line2
        slope2 = ylength2/xlength2
        const2 = line2[1] - slope2 * line2[0]
        xintersect = line1[0]
        yintersect = slope2 * xintersect + const2

    elif abs(xlength2) < EPS:
        # Infinite line2 slope special case
        # The lines intersect where the X coordinate of line2 is
        # substituted into the equation for line1
        slope1 = ylength1/xlength1
        const1 = line1[1] - slope1 * line1[0]
        xintersect = line2[0]
        yintersect = slope1 * xintersect + const1
        
    else:
        # Determine the slopes and constants of the two lines
        slope1 = ylength1/xlength1
        slope2 = ylength2/xlength2
        const1 = line1[1] - slope1 * line1[0]
        const2 = line2[1] - slope2 * line2[0]

        # If the slopes are the same, the lines are parallel and cannot
        # intersect
        slopediff = slope1 - slope2
        if abs(slopediff) < EPS:
            # Parallel lines
#             print("Parallel lines")
            return False

        # Find the point of intersection between the two lines
        # when they are extrapolated to infinity.
        xintersect = (const2 - const1) / (slope1 - slope2)
        yintersect = slope1 * xintersect + const1

#     print("Infinite lines intersect at", xintersect, yintersect)
        
    # The line segments intersect if this intersection point
    # occurs inside the range of both segments.
    left1 = min(line1[0],line1[2]) - tolerance
    right1 = max(line1[0],line1[2]) + tolerance
    left2 = min(line2[0],line2[2]) - tolerance
    right2 = max(line2[0],line2[2]) + tolerance
    bottom1 = min(line1[1],line1[3]) - tolerance
    top1 = max(line1[1],line1[3]) + tolerance
    bottom2 = min(line2[1],line2[3]) - tolerance
    top2 = max(line2[1],line2[3]) + tolerance
#     print("Line 1 extent is", left1, right1, bottom1, top1)
#     print("Line 2 extent is", left2, right2, bottom2, top2)
    if xintersect > left1 and xintersect < right1 and \
       xintersect > left2 and xintersect < right2 and \
       yintersect > bottom1 and yintersect < top1 and \
       yintersect > bottom2 and yintersect < top2:
#         print("Infinite intersection point is within this range")
        return True
        
    return False

def line_intersects_circle(x1, y1, x2, y2, xcen, ycen, radius):
    """
        
    Helper function to determine whether a specified line intersects
    with the circle described by (xcen,ycen,radius).
        
    :Parameters:

    x1: float
        X coordinate of first point
    y1: float
        Y coordinate of first point
    x2: float
        X coordinate of second point
    y2: float
        Y coordinate of second point
    xcen: float
        X coordinate of centre of circle
    ycen: float
        Y coordinate of centre of circle
    radius: float
        Radius of circle

    :Returns:

    True or False
        
    """
    # First determine the closest distance from the centre of the
    # circle to any point along the line
    distance = distance_from_point_to_line(xcen, ycen, x1, y1, x2, y2)
    
    # If the distance is less than the radius of the circle, the circle
    # intersects with the line.
    if distance < radius:
#         print("Line intersects circle (%f < %f)" % (distance, radius))
        return True
    else:
        return False

def distance_from_point_to_line(xpoint, ypoint, x1, y1, x2, y2):
    """
    
    Helper function which calculates the closest approach distance
    between a point and the specified line segment.
    
    """
    # Determine line equation from the end points
    (slope, intercept) = line_from_two_points(x1, y1, x2, y2)
#     print("Line slope=%s, intercept=%f" % (str(slope),intercept))
    if slope is not None:
    
        # Evaluate the line at the X coordinate of the point and
        # calculate the closest approach distance for an infinite
        # line. 
        yline = intercept + slope * xpoint
        ydiff = ypoint - yline       
        theta = math.atan(slope)
        distance = abs(ydiff) * math.cos(theta)
        # Calculate the coordinates of this closest point.
        xclosest = xpoint - (distance * math.sin(theta))
#         print("Closest approach is %f at x=%f" % (distance,xclosest) )
        
        # The infinite line distance is valid only if this closest
        # point actually lies on the line segment. If it doesn't,
        # the distance needs to be recalculated.
        if xclosest < min(x1,x2):
            # The closest point is off the left hand end. The real
            # distance is between the point and the left hand end.
            if x1 < x2:
                distsq = distance_squared(x1, y1, xpoint, ypoint)
            else:
                distsq = distance_squared(x2, y2, xpoint, ypoint)
            distance = math.sqrt(distsq)
#             print("Distance to left hand end is %f", distance)
        elif xclosest > max(x1,x2):
            # The closest point is off the right hand end. The real
            # distance is between the point and the right hand end.
            if x1 > x2:
                distsq = distance_squared(x1, y1, xpoint, ypoint)
            else:
                distsq = distance_squared(x2, y2, xpoint, ypoint)
            distance = math.sqrt(distsq)
#             print("Distance to right hand end is %f", distance)
    else:
        # Infinite slope
        if ypoint < min(y1,y2):
            # Point is below the bottom of the line
            if y1 < y2:
                distsq = distance_squared(x1, y1, xpoint, ypoint)
            else:
                distsq = distance_squared(x2, y2, xpoint, ypoint)
            distance = math.sqrt(distsq)
#             print("Distance to bottom end is %f", distance)
        elif ypoint > max(y1,y2):
            # Point is above the top of the line
            if y1 > y2:
                distsq = distance_squared(x1, y1, xpoint, ypoint)
            else:
                distsq = distance_squared(x2, y2, xpoint, ypoint)
            distance = math.sqrt(distsq)
#             print("Distance to top end is %f", distance)
        else:
            # Point is within the Y range of the line
            distance = abs(x1 - xpoint)
#             print("Closest approach is %f at y=%f" % (distance,ypoint) )
    return distance

def quadrangles_intersect(quad1, quad2):
    """
        
    Helper function to determine whether two quadrangles
    intersect.
        
    :Parameters:

    quad1: tuple of ((x1,y1),(x2,y2),(x3,y3),(x4,y4))
        The corner coordinates of the first quadrangle
    quad2: tuple of ((x1,y1),(x2,y2),(x3,y3),(x4,y4))
        The corner coordinates of the second quadrangle

    :Returns:

    True or False
        
    """
    assert isinstance(quad1, (list,tuple))
    assert isinstance(quad2, (list,tuple))
    assert len(quad1) > 3
    assert len(quad2) > 3
    assert isinstance(quad1[0], (list,tuple))
    assert isinstance(quad2[0], (list,tuple))
    assert len(quad1[0]) > 1
    assert len(quad2[0]) > 1
        
    # Two quadrangles intersect either if any of the lines bordering
    # those quadrangles intersect.
    line11 = (quad1[0][0], quad1[0][1], quad1[1][0], quad1[1][1])
    line12 = (quad1[1][0], quad1[1][1], quad1[2][0], quad1[2][1])
    line13 = (quad1[2][0], quad1[2][1], quad1[3][0], quad1[3][1])
    line14 = (quad1[3][0], quad1[3][1], quad1[0][0], quad1[0][1])

    line21 = (quad2[0][0], quad2[0][1], quad2[1][0], quad2[1][1])
    line22 = (quad2[1][0], quad2[1][1], quad2[2][0], quad2[2][1])
    line23 = (quad2[2][0], quad2[2][1], quad2[3][0], quad2[3][1])
    line24 = (quad2[3][0], quad2[3][1], quad2[0][0], quad2[0][1])

    for line1 in (line11, line12, line13, line14):
        for line2 in (line21, line22, line23, line24):
            if lines_intersect(line1, line2):
                # Two edges intersect
                return True

    # None of the edges intersect, but the quadrangles might still
    # intersect if one of them is entirely inside the other.
    xmin1 = min( quad1[0][0], quad1[1][0], quad1[2][0], quad1[3][0])
    xmax1 = max( quad1[0][0], quad1[1][0], quad1[2][0], quad1[3][0])
    ymin1 = min( quad1[0][1], quad1[1][1], quad1[2][1], quad1[3][1])
    ymax1 = max( quad1[0][1], quad1[1][1], quad1[2][1], quad1[3][1])
    xmin2 = min( quad2[0][0], quad2[1][0], quad2[2][0], quad2[3][0])
    xmax2 = max( quad2[0][0], quad2[1][0], quad2[2][0], quad2[3][0])
    ymin2 = min( quad2[0][1], quad2[1][1], quad2[2][1], quad2[3][1])
    ymax2 = max( quad2[0][1], quad2[1][1], quad2[2][1], quad2[3][1])

    if (xmin1 < xmin2) and (xmax1 > xmax2) and \
       (ymin1 < ymin2) and (ymax1 > ymax2):
        # Quadrangle 2 entirely within quadrangle 1
        return True

    if (xmin2 < xmin1) and (xmax2 > xmax1) and \
       (ymin2 < ymin1) and (ymax2 > ymax1):
        # Quadrangle 1 entirely within quadrangle 2
        return True
    
    # If all the above tests have failed the quadrangles do not intersect.
    return False

def triangles_intersect(tri1, tri2):
    """
        
    Helper function to determine whether two triangles
    intersect.
        
    :Parameters:

    tri1: tuple of ((x1,y1),(x2,y2),(x3,y3))
        The corner coordinates of the first triangle
    tri2: tuple of ((x1,y1),(x2,y2),(x3,y3))
        The corner coordinates of the second triangle

    :Returns:

    True or False
        
    """
    assert isinstance(tri1, (list,tuple))
    assert isinstance(tri2, (list,tuple))
    assert len(tri1) > 2
    assert len(tri2) > 2
    assert isinstance(tri1[0], (list,tuple))
    assert isinstance(tri2[0], (list,tuple))
    assert len(tri1[0]) > 1
    assert len(tri2[0]) > 1
        
    # Two triangles intersect if any point from one triangle lies within
    # the other.
    for xpoint,ypoint in tri1:
        if point_inside_triangle(xpoint, ypoint, tri2[0][0], tri2[0][1],
                tri2[1][0], tri2[1][1], tri2[2][0], tri2[2][1]):
            return True
    for xpoint,ypoint in tri2:
        if point_inside_triangle(xpoint, ypoint, tri1[0][0], tri1[0][1],
                tri1[1][0], tri1[1][1], tri1[2][0], tri1[2][1]):
            return True
    return False

def triangle_intersects_circle(triang, xcen, ycen, radius):
    """
        
    Helper function to determine whether a triangle and a
    circle intersect.
        
    :Parameters:

    triang: tuple of ((x1,y1),(x2,y2),(x3,y3))
        The corner coordinates of the triangle
    xcen: float
        X coordinate of centre of circle
    ycen: float
        Y coordinate of centre of circle
    radius: float
        Radius of circle

    :Returns:

    True or False
        
    """
    assert isinstance(triang, (list,tuple))
    assert isinstance(triang, (list,tuple))
    assert len(triang) > 2
    assert len(triang) > 2
    
    # A circle cannot intersect a triangle if it lies entirely outside
    # the rectangle bounding the triangle.
    xmin = min(triang[0][0], triang[1][0], triang[2][0]) - radius
    xmax = max(triang[0][0], triang[1][0], triang[2][0]) + radius
    ymin = min(triang[0][1], triang[1][1], triang[2][1]) - radius
    ymax = max(triang[0][1], triang[1][1], triang[2][1]) + radius
#     print("Circle centre is", xcen, ycen)
#     print("Bounds of triangle are", xmin, xmax, ymin, ymax)
    if xcen < xmin or xcen > xmax or ycen < ymin or ycen > ymax:
#         print("Circle entirely outside bounds of triangle")
        return False
   
    # A circle intersects a triangle if it intersects any of the
    # straight lines making up the edges of the triangle.
#     print("Triangle side 1-2", triang[0][0], triang[0][1],
#                               triang[1][0], triang[1][1])
    if line_intersects_circle(triang[0][0], triang[0][1],
                              triang[1][0], triang[1][1],
                              xcen, ycen, radius):
#         print("Triangle side 1-2 intersects circle")
        return True
#     print("Triangle side 2-3", triang[1][0], triang[1][1],
#                               triang[2][0], triang[2][1])
    if line_intersects_circle(triang[1][0], triang[1][1],
                              triang[2][0], triang[2][1],
                              xcen, ycen, radius):
#         print("Triangle side 2-3 intersects circle")
        return True
#     print("Triangle side 1-3", triang[0][0], triang[0][1],
#                               triang[2][0], triang[2][1])
    if line_intersects_circle(triang[0][0], triang[0][1],
                              triang[2][0], triang[2][1],
                              xcen, ycen, radius):
#         print("Triangle side 1-3 intersects circle")
        return True
    return False

def polygons_intersect(poly1, poly2):
    """
    
    Helper function to determine whether two polygons
    intersect.
        
    :Parameters:

    poly1: tuple of ((x1,y1),(x2,y2),...)
        The coordinates of the vertices of the first polygon.
        There must be at least 3 vertices.
    poly2: tuple of ((x1,y1),(x2,y2),(x3,y3),...)
        The coordinates of the vertices of the second polygon
        There must be at least 3 vertices.

    :Returns:

    True or False
    
    """
    assert isinstance(poly1, (list,tuple))
    assert isinstance(poly2, (list,tuple))
    assert len(poly1) > 2
    assert len(poly2) > 2
    assert isinstance(poly1[0], (list,tuple))
    assert isinstance(poly2[0], (list,tuple))
    assert len(poly1[0]) > 1
    assert len(poly2[0]) > 1
    
    # TO BE IMPLEMENTED
    raise NotImplementedError("Polygon intersection not implemented yet.")
    return True

def line_from_two_points(x1, y1, x2, y2):
    """
    
    Helper function to return the equation of a line
    passing through any two points.

    :Parameters:

    x1: float
        X coordinate of first point
    y1: float
        Y coordinate of first point
    x2: float
        X coordinate of second point
    y2: float
        Y coordinate of second point
    
    :Returns:
    
    (slope, intercept) or (None, xposition) if the slope is infinite.
    
    """
    xdiff = (x2-x1)
    if abs(xdiff) > 0.0:
        ydiff = (y2-y1)
        slope = ydiff / xdiff
        intercept = y1 - slope * x1
        return (slope, intercept)
    else:
        return (None, x1)

def y_coord_on_circle(xcoord, xcen, ycen, radius):
    """
    
    Return the Y coordinates of the intersection point between the
    vertical line at xcoord and the circle at centre (xcen, ycen)
    and given radius.
    
    """
    xdiff = xcoord - xcen
    # The problem is not solveable if the line is outside the circle.
    if xdiff < -radius or xdiff > radius:
        return None
    
    rsqxsq = (radius * radius) - (xdiff * xdiff)
    ydiff = math.sqrt(rsqxsq)
    return (ycen-ydiff, ycen+ydiff)

def point_inside_boundary(xpoint, ypoint, xmin, xmax, ymin, ymax):
    """
    
    Helper function to determine if a point (xpoint,ypoint) is
    inside a rectangle bounded by (xmin,xmax,ymin,ymax).
    
    """
    if xpoint > xmin and xpoint < xmax and ypoint > ymin and ypoint < ymax:
        return True
    else:
        return False

def point_inside_circle(xpoint, ypoint, xcen, ycen, radius, left_only=False,
                        right_only=False, bottom_only=False, top_only=False):
    """
    
    Helper function to determine if a point (xpoint,ypoint) is
    inside a circle described by (xcen,ycen,radius).

    Flags left_only, right_only, bottom_only and top_only can be used to
    tailor the function to check for the point lying within one half or one
    quarter of the circle.

    """
    # If the point is outside the bounding rectangle, it can't be
    # inside the circle.
    if right_only:
        xmin = xcen
    else:
        xmin = xcen - radius
    if left_only:
        xmax = xcen
    else:
        xmax = xcen + radius
    if top_only:
        ymin = ycen
    else:
        ymin = ycen - radius
    if bottom_only:
        ymax = ycen
    else:
        ymax = ycen + radius
    if point_inside_boundary(xpoint, ypoint, xmin, xmax, ymin, ymax):
        # Point is inside rectangle. Now determine if it is also
        # inside the circle.
        radiussq = radius * radius
        xdiff = xpoint - xcen
        ydiff = ypoint - ycen 
        distsq = (xdiff * xdiff) + (ydiff * ydiff)
        if distsq < radiussq:
            # Inside circle
            return True
        else:
            return False
    else:
        return False

def point_inside_ellipse(xpoint, ypoint, xcen, ycen, xsemimajor, ysemiminor,
                         orient,
                         left_only=False, right_only=False, bottom_only=False,
                         top_only=False):
    """
    
    Helper function to determine if a point (xpoint,ypoint) is
    inside an ellipse described by (xcen,ycen,xsemimajor,ysemiminor,orient).
    
    Flags left_only, right_only, bottom_only and top_only can be used to
    tailor the function to check for the point lying with one half or one
    quarter of the ellipse.
    
    """
    # FIXME: This function doesn't work properly when the ellipse is tilted.
    # If the point is outside the bounding rectangle, it can't be
    # inside the ellipse.
    if right_only:
        xmin = xcen
    else:
        xmin = xcen - xsemimajor
    if left_only:
        xmax = xcen
    else:
        xmax = xcen + xsemimajor
    if top_only:
        ymin = ycen
    else:
        ymin = ycen - ysemiminor
    if bottom_only:
        ymax = ycen
    else:
        ymax = ycen + ysemiminor
    if point_inside_boundary(xpoint, ypoint, xmin, xmax, ymin, ymax):
        # Point is inside rectangle. Now determine if it is also
        # inside the ellipse.
        sxmajsq = xsemimajor * xsemimajor
        syminsq = ysemiminor * ysemiminor
        xdiff = xpoint - xcen
        ydiff = ypoint - ycen
        (xnew, ynew) = rotate_coordinates(xdiff, ydiff, -orient)
        xratio = (xnew * xnew) / sxmajsq
        yratio = (ynew * ynew) / syminsq
        if (xratio + yratio) < 1.0:
            # Inside ellipse
            return True
        else:
            return False
    else:
        return False

def point_inside_triangle(xpoint, ypoint, ax, ay, bx, by, cx, cy):
    """

    Helper function to determine if a point (xpoint,ypoint)
    is inside the triangle bounded by the points A (ax,ay),
    B (bx,by) and C (cx,cy).

    :Parameters:

    xpoint: float
        X coordinate of point
    ypoint: float
        Y coordinate of point
    ax: float
        X coordinate of first corner of triangle
    ay: float
        Y coordinate of first corner of triangle
    bx: float
        X coordinate of second corner of triangle
    by: float
        Y coordinate of second corner of triangle
    cx: float
        X coordinate of third corner of triangle
    cy: float
        Y coordinate of third corner of triangle

    :Returns:

    True or False
    
    """
#     print("Testing for", xpoint, ypoint, "inside", ax, ay, bx, by, cx, cy)
    # This function is reused 3 times
    def test_point(slope, intercept, xtriangle, ytriangle, xpoint, ypoint):
        if slope is not None:
            # Is the third point above or below this line?
            testty = slope * xtriangle + intercept
            if ytriangle < testty:
                # The third point is below the AB line. The given point must
                # also be below this line to be inside the triangle.
                testpy = slope * xpoint + intercept
                if ypoint > testpy:
                    # Point outside triangle
#                     print("Point above triangle (%f > %f)" % (ypoint, testpy))
                    return False
            else:
                # The third point is above the AB line. The given point must
                # also be above this line to be inside the triangle.
                testpy = slope * xpoint + intercept
                if ypoint < testpy:
                    # Point outside triangle
#                     print("Point below triangle (%f < %f)" % (ypoint, testpy))
                    return False
        else:
            # The side of the triangle is vertical.
            # Is the third point to the left or right of this line?
            if xtriangle < intercept:
                # The third point is to the left of the AB line. The given
                # point must also be on the left of this line to be inside
                # the triangle.
                if xpoint > intercept:
                    # Point outside triangle
#                     print("Point to right of triangle (%f > %f)" % (xpoint, intercept))
                    return False
            else:
                # The third point is to the right of the AB line. The given
                # point must also be on the right of this line to be inside
                # the triangle.
                if xpoint < intercept:
                    # Point outside triangle
#                     print("Point to left of triangle (%f < %f)" % (xpoint, intercept))
                    return False
#         print("Point inside triangle side")
        # Ok so far
        return True
    
    # Determine the equation of the line A-B and compare with C
    (slope, intercept) = line_from_two_points(ax, ay, bx, by)
    inside = test_point(slope, intercept, cx, cy, xpoint, ypoint)
    if not inside:
        return False
    # Determine the equation of the line B-C and compare with A
    (slope, intercept) = line_from_two_points(bx, by, cx, cy)
    inside = test_point(slope, intercept, ax, ay, xpoint, ypoint)
    if not inside:
        return False
    # Determine the equation of the line A-C and compare with B
    (slope, intercept) = line_from_two_points(ax, ay, cx, cy)
    inside = test_point(slope, intercept, bx, by, xpoint, ypoint)
    if not inside:
        return False
#     print("Point is within all 3 sides, so is inside triangle")
    return True

def point_inside_avoidance_zone(xpoint, ypoint, xfibre, yfibre, favoid,
                                length1, length2, avoid2, halfwidth2):
    """

    Helper function to determine if a point (xpoint,ypoint)
    is inside the avoidance zone for a positioner arm assigned
    to a target (xfibre,yfibre). All coordinates are referenced
    to an origin at the positioner centre (same as the alpha axis).
    
    OBSOLETE FUNCTION.

    :Parameters:

    xpoint: float
        X coordinate of point
    ypoint: float
        Y coordinate of point
    xfibre: float
        X coordinate of fibre
    yfibre: float
        Y coordinate of fibre.
    favoid: float
        Radius of fibre avoidance circle. No two targets can come closer
        than this limit.
    length1: float
        Length of the alpha arm (determining the location of the elbow joint).
    length2: float
        Length of the beta arm (determining the distance of the fibre
        from the elbow joint).
    avoid2: float
        Length of the outer portion of the beta arm which can collide with
        another beta arm.
    halfwidth2: float
        Width of the outer portion of the beta arm, which can collide with
        another beta arm.
        
    """
    # If xfibre is negative, reflect everything around the X axis.
    if xfibre < 0.0:
        xfibre = -xfibre
        xpoint = -xpoint

    # If yfibre is not zero, rotate all the coordinates so they are in
    # a frame of reference where yfibre is zero.
    if abs(yfibre) > EPS:
        angle = math.atan2(yfibre, xfibre)
        xfibre = math.sqrt( (xfibre * yfibre) + (yfibre * yfibre) )
        
        newxpoint = xpoint * math.cos(angle) + ypoint * math.sin(angle)
        newypoint = xpoint * math.sin(angle) + ypoint * math.cos(angle)
        xpoint = newxpoint
        ypoint = newypoint
        
    # Calculate some significant locations
    # Left-most extent of beta arm when tucked in, which is also the
    # left-most extent of ellipse swept out by beta arm.
    xleft = (length2 - length1) - avoid2
    
    # Right-most extent of ellipse swept out by beta arm
    xright = length2 + length1 - avoid2
    
    # Centre of avoidance ellipse
    xboundary = (xleft + xright)/2.0
    
    # Bottom of avoidance ellipse, determined by the Y extent of the
    # beta arm when the alpha arm is vertical.
    ybottom = -length1 * avoid2 / length2
    
    # The avoidance zone is oversized either by the fibre avoidance
    # or by half the width of the beta arm.
    xfarleft = xleft - favoid
    xfarright = xfibre + favoid
    yfarbottom = ybottom - halfwidth2
    yfartop = halfwidth2
    
    # If the test point is outside this extreme bounding rectangle then it
    # can't be inside the avoidance zone.
    if not point_inside_boundary(xpoint, ypoint, xfarleft, xfarright,
                                  yfarbottom, yfartop):
        return False
    
    # A more detailed check is needed. If the test point is inside any
    # of the following zones then it is inside the avoidance zone. The
    # easiest zones are tested first for efficiency
    
    # (1) Rectangle
    if point_inside_boundary(xpoint, ypoint, xleft, xfibre, -halfwidth2,
                              halfwidth2):
        return True
    
    # (2) Circles at either end of the travel
    if point_inside_circle(xpoint, ypoint, xleft, 0.0, favoid):
        return True
    if point_inside_circle(xpoint, ypoint, xfibre, 0.0, favoid):  # favoid or halfwidth2?
        return True
    
    # (3) Quarter ellipse
    semi_major = favoid + (xright - xleft)/2.0
    semi_minor = abs(ybottom) + halfwidth2
    if point_inside_ellipse(xpoint, ypoint, xboundary, 0.0, semi_major, 
                            semi_minor, 0.0, left_only=True, bottom_only=True):
        return True
    
    # (4) Triangle
    if point_inside_triangle(xpoint, ypoint, xboundary, yfarbottom, 
                             xboundary, -halfwidth2, xfibre, -halfwidth2):
        return True
    
    return False

def generate_avoidance_perimeter(xfibre, length1, length2, mlength, mwidth,
                                 favoid, padding=0.0, inc1=1.0, inc2=0.5):
    """
    
    Generate a list of X,Y coordinates defining the perimeter of the
    avoidance zone swept out by the metrology avoidance zone and the
    fibre avoidance zone when the fibre holder is moved radially from
    xfibre to the minimum radius dicated by the arm lengths.

    OBSOLETE FUNCTION.
    
    :Parameters:
    
    xfibre: float
        The fibre X coordinate for which to generate an avoidance
        perimeter.
    length1: float
        Length of inner (alpha) arm
    length2: float
        Length of outer (beta) arm
    mlength: float
        Length of metrology zone at end of beta arm.
    mwidth: float
        Width of metrology zone at end of beta arm.
    favoid: float
        Radius of the avoidance zone surrounding the fibre.
    padding: float (optional)
        An optional padding to be applied to the edges of the
        avoidance zone, increasing the clearance.
        By default the padding is zero.
    inc1: float (optional)
        The X coordinate increment for the avoidance zone.
        By default this is 2.0.
    inc2: float (optional)
        The X coordinate increment for the xfibre tests.
        By default this is 0.5.

    :Returns:
    
    (xpoints, ypoints): tuple of tuples
        A list of x,y coordinates defining the perimeter of the
        avoidance zone at the given xfibre.
    
    """
    # First generate a list of X coordinates to be tested.
    maxlength = length2 + length1
    minlength = length2 - length1
    xmin = 0.0
    xmax = maxlength + padding + favoid
    xcoords = np.arange(xmin, xmax+inc1, inc1)
    # Fill the Y extremities with dummay values
    yminlist = 1000.0 * np.ones_like(xcoords)
    ymaxlist = -1000.0 * np.ones_like(xcoords)
    
    # Generate the avoidance zone swept out by the metrology zone
    # as the fibre moves inwards.
    if xfibre is None:
        xfibre = maxlength
    xtest = xfibre
    while (xtest >= minlength):
#         print "Fibre moves to=", xtest
        (corner_bottom, corner_left, corner_top, corner_right) = \
            solve_metrology_edge(xtest, length1, length2, mlength, mwidth,
                                 padding=padding)
        bottom_line = line_from_two_points(corner_bottom[0], corner_bottom[1],
                                           corner_right[0], corner_right[1])
        left_line = line_from_two_points(corner_left[0], corner_left[1],
                                         corner_bottom[0], corner_bottom[1])
        top_line = line_from_two_points(corner_left[0], corner_left[1],
                                        corner_top[0], corner_top[1])
        right_line = line_from_two_points(corner_top[0], corner_top[1],
                                          corner_right[0], corner_right[1])

        for ii in range(0, len(xcoords)):
#             print "Testing", xcoords[ii]
            if bottom_line[0] is not None:
                if xcoords[ii] >= corner_bottom[0] and \
                   xcoords[ii] <= corner_right[0]:
#                     print "Inside bottom line range"
                    ybottom = bottom_line[1] + bottom_line[0] * xcoords[ii]
                    if ybottom < yminlist[ii]:
                        yminlist[ii] = ybottom
            if left_line[0] is not None:
                if xcoords[ii] >= corner_left[0] and \
                   xcoords[ii] <= corner_bottom[0]:
#                     print "Inside left line range"
                    ybottom = left_line[1] + left_line[0] * xcoords[ii]
                    if ybottom < yminlist[ii]:
                        yminlist[ii] = ybottom
            if top_line[0] is not None:
#                 print "Inside top line range"
                if xcoords[ii] >= corner_left[0] and \
                   xcoords[ii] <= corner_top[0]:
                    ytop = top_line[1] + top_line[0] * xcoords[ii]
                    if ytop > ymaxlist[ii]:
                        ymaxlist[ii] = ytop
            if right_line[0] is not None:
                if xcoords[ii] >= corner_top[0] and \
                   xcoords[ii] <= corner_right[0]:
#                     print "Inside right line range"
                    ytop = right_line[1] + right_line[0] * xcoords[ii]
                    if ytop > ymaxlist[ii]:
                        ymaxlist[ii] = ytop
        xtest -= inc2
        
    # Add the fibre avoidance zone if needed.
    if favoid > 0.0:
        for ii in range(0, len(xcoords)):
            result = y_coord_on_circle(xcoords[ii], xfibre, 0.0, favoid+padding)
            if result is not None:
                if result[0] < yminlist[ii]:
                    yminlist[ii] = result[0]
                if result[1] > ymaxlist[ii]:
                    ymaxlist[ii]  = result[1]
                
    xpoints = []
    ypoints = []
    for ii in range(0, len(xcoords)):
        if ymaxlist[ii] > -999.0:
            xpoints.append(xcoords[ii])
            ypoints.append(ymaxlist[ii])
    for ii in range(0, len(xcoords)):
        jj = -ii - 1
        if yminlist[jj] < 999.0:
            xpoints.append(xcoords[jj])
            ypoints.append(yminlist[jj])
    xpoints.append(xpoints[0])
    ypoints.append(ypoints[0])
    return (xpoints,ypoints)


def generate_avoidance_perimeter2(xstart, ystart, xfinish, yfinish, 
                                  length1, length2, mlength, mwidth,
                                  favoid, padding=0.0, inc1=1.0, inc2=0.5):
    """
    
    Generate a list of X,Y coordinates defining the perimeter of the
    avoidance zone swept out by the metrology avoidance zone and the
    fibre avoidance zone when the fibre holder moves from
    (xstart, ystart) to (xfinish, yfinish).

    OBSOLETE FUNCTION.
    
    :Parameters:
    
    xstart: float
        The fibre X coordinate at the beginning of travel.
    ystart: float
        The fibre Y coordinate at the beginning of travel.
    xfinish: float
        The fibre X coordinate at the end of travel.
    yfinish: float
        The fibre Y coordinate at the end of travel.
    length1: float
        Length of inner (alpha) arm
    length2: float
        Length of outer (beta) arm
    mlength: float
        Length of metrology zone at end of beta arm.
    mwidth: float
        Width of metrology zone at end of beta arm.
    favoid: float
        Radius of the avoidance zone surrounding the fibre.
    padding: float (optional)
        An optional padding to be applied to the edges of the
        avoidance zone, increasing the clearance.
        By default the padding is zero.
    inc1: float (optional)
        The X coordinate increment for the avoidance zone.
        By default this is 2.0.
    inc2: float (optional)
        The X coordinate increment for the xfibre tests.
        By default this is 0.5.

    :Returns:
    
    (xpoints, ypoints): tuple of tuples
        A list of x,y coordinates defining the perimeter of the
        avoidance zone at the given xfibre.
    
    """
    # First generate a list of X coordinates to be tested.
    maxlength = length2 + length1
    minlength = length2 - length1
    xmin = min( minlength, xstart, xfinish ) - mlength - padding - favoid
    xmax = max(maxlength, xstart, xfinish) + mlength + padding + favoid
    xcoords = np.arange(xmin, xmax+inc1, inc1)
    # Fill the Y extremities with dummay values
    yminlist = 1000.0 * np.ones_like(xcoords)
    ymaxlist = -1000.0 * np.ones_like(xcoords)
    
    # Generate the avoidance zone swept out by the metrology zone
    # as the fibre moves from start to finish.
    xtest = xstart
    ytest = ystart
    npoints = int(0.5 + abs(xfinish - xstart) / inc2)
    if xstart > xfinish:
        xinc = -inc2
    else:
        xinc = inc2
    yinc = (yfinish - ystart) / float(npoints)  
    for point in range(0,npoints):
        # Convert to polar coordinates
        (rtest, thtest) = cartesian_to_polar( xtest, ytest )
#         print "Fibre moves to=", xtest, ytest, rtest, thtest
        # Solve the metrology edge for a nominal location at x=rtest, theta=0
        # and then rotate the coordinates by thtest
        (cb, cl, ct, cr) = \
            solve_metrology_edge(rtest, length1, length2, mlength, mwidth,
                                 padding=padding)
        corner_bottom = rotate_coordinates( cb[0], cb[1], thtest)
        corner_left = rotate_coordinates( cl[0], cl[1], thtest)
        corner_top = rotate_coordinates( ct[0], ct[1], thtest)
        corner_right = rotate_coordinates( cr[0], cr[1], thtest)
        bottom_line = line_from_two_points(corner_bottom[0], corner_bottom[1],
                                           corner_right[0], corner_right[1])
        left_line = line_from_two_points(corner_left[0], corner_left[1],
                                         corner_bottom[0], corner_bottom[1])
        top_line = line_from_two_points(corner_left[0], corner_left[1],
                                        corner_top[0], corner_top[1])
        right_line = line_from_two_points(corner_top[0], corner_top[1],
                                          corner_right[0], corner_right[1])

        for ii in range(0, len(xcoords)):
#             print "Testing", xcoords[ii]
            if bottom_line[0] is not None:
                if xcoords[ii] >= corner_bottom[0] and \
                   xcoords[ii] <= corner_right[0]:
#                     print "Inside bottom line range"
                    ybottom = bottom_line[1] + bottom_line[0] * xcoords[ii]
                    if ybottom < yminlist[ii]:
                        yminlist[ii] = ybottom
            if left_line[0] is not None:
                if xcoords[ii] >= corner_left[0] and \
                   xcoords[ii] <= corner_bottom[0]:
#                     print "Inside left line range"
                    ybottom = left_line[1] + left_line[0] * xcoords[ii]
                    if ybottom < yminlist[ii]:
                        yminlist[ii] = ybottom
            if top_line[0] is not None:
#                 print "Inside top line range"
                if xcoords[ii] >= corner_left[0] and \
                   xcoords[ii] <= corner_top[0]:
                    ytop = top_line[1] + top_line[0] * xcoords[ii]
                    if ytop > ymaxlist[ii]:
                        ymaxlist[ii] = ytop
            if right_line[0] is not None:
                if xcoords[ii] >= corner_top[0] and \
                   xcoords[ii] <= corner_right[0]:
#                     print "Inside right line range"
                    ytop = right_line[1] + right_line[0] * xcoords[ii]
                    if ytop > ymaxlist[ii]:
                        ymaxlist[ii] = ytop
        xtest += xinc
        ytest += yinc
        
    # Add the fibre avoidance zone if needed.
    if favoid > 0.0:
        xtest = xstart
        ytest = ystart
        for point in range(0,npoints):
            for ii in range(0, len(xcoords)):
                result = y_coord_on_circle(xcoords[ii], xtest, ytest, favoid+padding)
                if result is not None:
                    if result[0] < yminlist[ii]:
                        yminlist[ii] = result[0]
                    if result[1] > ymaxlist[ii]:
                        ymaxlist[ii]  = result[1]
            xtest += xinc
            ytest += yinc
        for ii in range(0, len(xcoords)):
            result = y_coord_on_circle(xcoords[ii], xfinish, yfinish, favoid+padding)
            if result is not None:
                if result[0] < yminlist[ii]:
                    yminlist[ii] = result[0]
                if result[1] > ymaxlist[ii]:
                    ymaxlist[ii]  = result[1]
               
    xpoints = []
    ypoints = []
    for ii in range(0, len(xcoords)):
        if ymaxlist[ii] > -999.0:
            xpoints.append(xcoords[ii])
            ypoints.append(ymaxlist[ii])
    for ii in range(0, len(xcoords)):
        jj = -ii - 1
        if yminlist[jj] < 999.0:
            xpoints.append(xcoords[jj])
            ypoints.append(yminlist[jj])
    xpoints.append(xpoints[0])
    ypoints.append(ypoints[0])
    return (xpoints,ypoints)

# def xfibre_to_avoidance(xfibre, length1, length2, mlength, mwidth, padding=0.0):
#     """
#     
#     Given a particular X fibre coordinate, return the coodinates of
#     the left-most, right-most, bottom-most and top-most corners of
#     the metrology avoidance zone.
#     
#     :Parameters:
#     
#     xfibre: float
#         X coordinate of the centre location of the fibre.
#         The Y coordinate is assumed to be zero (i.e. the avoidance
#         zone is predicted for fibre movement along the X axis - other
#         zones may be calculated by rotating this one).
#     length1: float
#         Length of inner (alpha) arm
#     length2: float
#         Length of outer (beta) arm
#     mlength: float
#         Length of metrology zone at end of beta arm.
#     mwidth: float
#         Width of metrology zone at end of beta arm.
#     padding: float (optional)
#         An optional padding to be applied to the edges of the
#         avoidance zone, increasing the clearance.
#         By default the padding is zero.
# 
#     :Returns:
#     
#     quad: tuple of ((x1,y1),(x2,y2),(x3,y3),(x4,y4))
#         The corner coordinates of the metrology zone, sorted
#         into left-most, right-most, bottom-most and top-most.
#     
#     """
#     # First solve the rectangle
#     quad = solve_metrology_edge(xfibre, length1, length2, mlength, mwidth,
#                                 padding=padding)
#     
#     # Sort the corners into order.
#     corner_left = quad[0]
#     corner_right = quad[0]
#     corner_bottom = quad[0]
#     corner_top = quad[0]
#     for corner in quad:
#         if corner[0] < corner_left[0]:
#             corner_left = corner
#         elif corner[0] > corner_right[0]:
#             corner_right = corner
#         if corner[1] < corner_bottom[0]:
#             corner_bottom = corner
#         elif corner[1] > corner_top[0]:
#             corner_top = corner
#             
#     return (corner_left, corner_right, corner_bottom, corner_top)

def solve_metrology_edge(xfibre, length1, length2, mlength, mwidth,
                         padding=0.0):
    """
    
    Give a fibre location (assumed to be at the end of the metrology zone)
    and the positioner arm parameters, find the location of the corners
    of a metrology rectangle of given length and width.
    
    :Parameters:
    
    xfibre: float
        X coordinate of the centre location of the fibre.
        The Y coordinate is assumed to be zero (i.e. the avoidance
        zone is predicted for fibre movement along the X axis - other
        zones may be calculated by rotating this one).
    length1: float
        Length of inner (alpha) arm
    length2: float
        Length of outer (beta) arm
    mlength: float
        Length of metrology zone at end of beta arm.
    mwidth: float
        Width of metrology zone at end of beta arm.
    padding: float (optional)
        An optional padding to be applied to the edges of the
        avoidance zone, increasing the clearance.
        By default the padding is zero.
    
    :Returns:
    
    quad: tuple of ((x1,y1),(x2,y2),(x3,y3),(x4,y4))
        The corner coordinates of the metrology zone
    
    """
    # First solve the triangle rule to determine the angle of orientation
    # of the beta arm.
    length2sq = length2 * length2
    xfibresq = xfibre * xfibre
    length1sq = length1 * length1
    angle = solve_triangle(length2, length2sq, xfibre, xfibresq, length1sq)
    
    if angle is not None:
    
        # Now solve the metrology rectangle to determine the location of the
        # corners
        quad = solve_tilted_rectangle(xfibre, 0.0, mlength, mwidth, angle,
                                      padding=padding, pad_upper=False)
    else:
        strg = "Triangle equation cannot be solved"
        raise ValueError(strg)
    return quad

def solve_tilted_rectangle(xend, yend, length, width, angle, padding=0.0,
                           pad_upper=True):
    """
    
    Given a rectangle of a certain length, width and orientation,
    knowing the coordinates of the centre of one end of the
    rectangle, return the coordinates of the corners.

    :Parameters:
    
    xend: float
        X coordinate of the centre of the upper edge (at the extremity
        of the length) of the rectangle.
    yend: float
        Y coordinate of the centre of the same edge of the rectangle.
    length: float
        Length of the rectangle
    width: float
        Width of the rectangle
    angle: float
        Angle of the rectangle (radians).
    padding: float (optional)
        An optional padding to be applied to the edges of the
        rectangle, increasing the length and the width by
        2 * padding. This parameter can be used to determine the
        corners of a new rectangle which avoids the edges of the
        original rectangle by at least this amount.
        By default the padding is zero and the corners of the
        original rectangle are returned.
    pad_upper: boolean (optional)
        Set True (the default) to pad the upper edge of the rectangle.
        Setting this to False allows one end of the rectangle to
        have a much smaller padding.

    :Returns:
    
    quad: tuple of ((x1,y1),(x2,y2),(x3,y3),(x4,y4))
        The corner coordinates of the rectangle
   
    """
    assert float(length) > 0.0
    assert float(width) > 0.0
    
    # The coordinates of the other edge of the rectangle can be calculated
    # from the length and orientation.
    xlength = length * math.cos(angle)
    ylength = length * math.sin(angle)
    xother = xend - xlength
    yother = yend - ylength
    
    # The X and Y increments of the corners from these ends depend on
    # the width and orientation
    xwidth2 = width * math.sin(angle) / 2.0
    ywidth2 = width * math.cos(angle) / 2.0
        
    x1 = xother + xwidth2
    y1 = yother - ywidth2
    x2 = xother - xwidth2
    y2 = yother + ywidth2
    x3 = xend - xwidth2
    y3 = yend + ywidth2
    x4 = xend + xwidth2
    y4 = yend - ywidth2
    
    # If required, apply a padding to the corner coordinates.
    if padding > 0.0:
        xlength_pad = padding * math.cos(angle)
        ylength_pad = padding * math.sin(angle)
        xwidth_pad = padding * math.sin(angle)
        ywidth_pad = padding * math.cos(angle)
        
        x1 = x1 - xlength_pad + xwidth_pad
        y1 = y1 - ylength_pad - ywidth_pad
        x2 = x2 - xlength_pad - xwidth_pad
        y2 = y2 - ylength_pad + ywidth_pad
        if pad_upper:
            x3 = x3 + xlength_pad - xwidth_pad
            y3 = y3 + ylength_pad + ywidth_pad
            x4 = x4 + xlength_pad + xwidth_pad
            y4 = y4 + ylength_pad - ywidth_pad
        else:
            # Only pad the width at the upper end of the rectangle
            x3 = x3 - xwidth_pad
            y3 = y3 + ywidth_pad
            x4 = x4 + xwidth_pad
            y4 = y4 - ywidth_pad
     
    quad = [(x1,y1), (x2,y2), (x3,y3), (x4,y4)]
    return quad

def solve_tangent_angle(distance, radius):
    """
    
    Helper function to calculate the angle between the
    centre of a circle and the tangent point, as seen from
    a point a certain distance from the circle.
    
    :Parameters:
    
    distance: float
        Distance of point from centre of circle.
    radius: float
        Radius of circle

    :Returns:
    
    tangent_angle: float
        The tangent angle in radians, or None if there is
        no solution.
    
    """
    sinangle = float(radius) / float(distance)
    if abs(sinangle) <= 1.0:
        angle = math.asin(sinangle)
    else:
        angle = None
    return angle

def solve_triangle(side1, side1sq, side2, side2sq, side3sq):
    """
        
    Helper function to solve the triangular cosine rule.
    Calculate an angle given all three sides.
        
    c**2 = a**2 + b**2 - 2 * a * b cos( C )

    :Parameters:
    
    side1: float
        Length of first side of triangle
    side1sq: float
        Square of the length of the first side of the triangle
        (given as well to save computation time).
    side2: float
        Length of the second side of the triangle
    side2sq: float
        Square of the length of the second side of the triangle
        (given as well to save computation time).
    side3sq: float
        Square of the length of the third side of the triangle
        
    :Returns:
    
    angle3: float
        Angle subtended by the third side of the triangle (radians)
        
    """
    if abs(side1) < EPS or abs(side2) < EPS:
        # Equation cannot be solved
        strg = "Triangle equation with side1=%f, side1sq=%f, " \
            "side2=%f, side2sq=%f, side3sq=%f " % \
            (side1, side1sq, side2, side2sq, side3sq)
        strg += "cannot be solved. Divide by zero!"
        warnings.warn(strg)
        return None
    
    cosangle = (side1sq + side2sq - side3sq) / \
            (2.0 * side1 * side2)
            
    if cosangle >= -1.0 and cosangle <= 1.0:
        angle = math.acos(cosangle)
        return angle
    else:
        # Equation cannot be solved
        strg = "Triangle equation with side1=%f, side1sq=%f, " \
            "side2=%f, side2sq=%f, side3sq=%f " % \
            (side1, side1sq, side2, side2sq, side3sq)
        strg += "cannot be solved. cos(angle) = %f" % cosangle
        warnings.warn(strg)
        return None

def solve_shoulder(theta_local, shoulder_fibre, parity):
    """
        
    Helper function to calculate shoulder angle from fibre angle
    (theta_local) and fibre angle at the shoulder joint (shoulder_fibre).
    In the polar coordinates, theta is an angle measured clockwise from
    the Y axis.

    :Parameters:

    theta_local: float
        Angular coordinate of fibre holder, clockwise from Y axis (radians)
    shoulder_fibre: float
        Shoulder to fibre angle (radians)
    parity: int (optional)
        The elbow option to be calculated:
            
        * 1 means elbow right armed
        * -1 means elbow left armed
       
    :Returns:
    
    angle1: float
        The shoulder angle (or alpha motor angle) (radians)
        
    """
    # There are two possible parities for the shoulder angle,
    # depending on which quadrant the elbow joint is contained in.
    if parity == PARITY_RIGHT:
        # Right armed orientation
        angle1 = (math.radians(90.0) - theta_local) - shoulder_fibre
    else:
        # Left armed orientation
        angle1 = (math.radians(90.0) - theta_local) + shoulder_fibre

    return angle1

def solve_shoulder_xy(xpos, ypos, shoulder_fibre, parity):
    """
        
    Helper function to calculate shoulder angle from fibre location
    (xpos, ypos) and fibre angle at the shoulder joint (shoulder_fibre).
    In the polar coordinates, theta is an angle measured clockwise from
    the Y axis.

    :Parameters:

    xpos: float
        X coordinate of fibre holder
    ypos: float
        Y coordinate of fibre holder
    shoulder_fibre: float
        Shoulder to fibre angle (or beta motor angle) (radians)
    parity: int (optional)
        The elbow parity to be calculated:
            
        * 1 means elbow right armed
        * -1 means elbow left armed
       
    :Returns:
    
    angle1: float
        The shoulder angle (or alpha motor angle) (radians)
        
    """
    # Determine the fibre angle.
    fibre_angle = math.atan2(ypos, xpos)

    # There are two possible parities for the shoulder angle,
    # depending on which quadrant the elbow joint is contained in.
    if parity == PARITY_RIGHT:
        # Right armed parity
        angle1 = fibre_angle - shoulder_fibre
    else:
        # Left armed parity
        angle1 = fibre_angle + shoulder_fibre

    return angle1

def solve_elbow(r_local, theta_local, parity, length1, length1sq, length2, length2sq):
    """
    
    Helper function to solve the elbow location given a target position
    (R, theta), parity choice (-1 or 1) and arm lengths.
    In the polar coordinates, theta is an angle measured clockwise from
    the Y axis.

    :Parameters:

    r_local: float
        Radial coordinate of fibre holder
    theta_local: float
        Angular coordinate of fibre holder (radians)
    parity: int (optional)
        The elbow option to be adopted:
            
        * 1 means elbow right armed
        * -1 means elbow left armed
    length1: float
        Length of inner (alpha) arm
    length1sq: float
        Square of the length of the inner (alpha) arm
        (given as well to save computation time).
    length2: float
        Length of outer (beta) arm
    length2sq: float
        Square of the length of the outer (beta) arm
        (given as well to save computation time).
        
    :Returns:
    
    (r_elbow,theta_elbow) - the location of the elbow joint,
    or (None,None) if a the equation cannot be solved.

    """
    # Solve the triangle cosine rule to determine the angle between
    # the shoulder and the fibre target.
    r_squared = r_local * r_local
    shoulder_fibre = solve_triangle(r_local, r_squared, length1, length1sq,
                                    length2sq )
    if shoulder_fibre is None:
        # Equation cannot be solved
        return (None, None)

    # Convert shoulder to fibre angle into shoulder angle
    # and return the elbow location in (R,theta).
    angle1 = solve_shoulder(theta_local, shoulder_fibre, parity)       
    return (length1, angle1)

def solve_elbow_xy(xpos, ypos, parity, length1, length1sq, length2, length2sq):
    """
    
    Helper function to solve the elbow location given a target position
    (xpos, ypos), parity choice (-1 or 1) and arm radii.
    In the polar coordinates, theta is an angle measured clockwise from
    the Y axis.

    :Parameters:

    xpos: float
        X coordinate of fibre holder
    ypos: float
        Y coordinate of fibre holder
    parity: int (optional)
        The elbow parity to be adopted:
            
        * 1 means elbow right armed
        * -1 means elbow left armed
    length1: float
        Length of inner (alpha) arm
    length1sq: float
        Square of the length of the inner (alpha) arm
        (given as well to save computation time).
    length2: float
        Length of outer (beta) arm
    length2sq: float
        Square of the length of the outer (beta) arm
        (given as well to save computation time).
        
    :Returns:
    
    (xelbow,yelbow) - the location of the elbow joint,
    or (None,None) if a the equation cannot be solved.

    """
    # Determine the third length of the triangle that needs
    # to be made by the arm's shoulder and elbow to reach
    # the fibre position.
    reachsq = (xpos * xpos) + (ypos * ypos)
    reach = math.sqrt(reachsq)

    # Solve the triangle cosine rule to determine the angle between
    # the shoulder and the fibre target.
    shoulder_fibre = solve_triangle(reach, reachsq, length1, length1sq,
                                    length2sq )
    if shoulder_fibre is None:
        # Equation cannot be solved
        return (None, None)

    # Convert shoulder to fibre angle into shoulder angle
    # and solve the elbow position.
    angle1 = solve_shoulder_xy(xpos, ypos, shoulder_fibre, parity)        
    xelbow = length1 * math.cos(angle1)
    yelbow = length1 * math.sin(angle1)
    return (xelbow,yelbow)

def read_grid_configuration(filename, pitch, ignore_header=True):
    """
    
    Read a CSV file and return a configuration list describing the locations
    of fibre positioners on a hexagonal grid.
    
    :Parameters:
        
    filename: string
        Name of the CSV file to be read.
    pitch: float
        Distance between neighbouring positioners (in mm)
    ignore_header: bool, optional
        If True (default), ignore the first line of the file.

    :Returns:
        
    config_list: list of (index, xcen, ycen, orient, column, row)
        List of configuration parameters describing the IDs,
        locations and orientations of the fibre positioners.
    
    """
    import csv    
    config_list = []
    mincol = 0
    minrow = 0
    
    with open(filename, 'r') as csvfile:
        cfile = csv.reader(csvfile)
        if ignore_header:
            cfile.next()
        for row in cfile:
            ident = int(row[0])
            xcen = float(row[4])
            ycen = float(row[5])
            orient = 0.0
            (column, row) = cartesian_to_hexagonal(xcen, ycen, pitch)
            mincol = min(column, mincol)
            minrow = min(row, minrow)
            config = [ident, xcen, ycen, orient, column, row]
            config_list.append(config)
    csvfile.close()
    
    # Adjust to make the minimum row and column zero.
    for config in config_list:
        config[4] -= mincol
        config[5] -= minrow
    return config_list

def remove_grid_cells( config_list, skip_cells ):
    """
    
    Remove the given list of cells from a configuration list
    
    :Parameters:
    
    config_list: list of (index, xcen, ycen, orient, column, row)
        List of configuration parameters describing the IDs,
        locations and orientations of the fibre positioners.
    skip_cells: list of tuples of 2 ints (optional)
        A list of (column,row) combinations to be removed. This allows
        gaps to be left in the grid to simulate the location of
        acquisition cameras. If empty, the list is unchanged.
    
    :Returns:
        
    new_config_list: list of (index, xcen, ycen, orient, column, row)
        List of configuration parameters describing the IDs,
        locations and orientations of the fibre positioners, with
        the skipped cells removed.
    
    """
    if skip_cells:
        new_config_list = []
        for config in config_list:
            skip_this = False
            column = config[4]
            row = config[5]
            for testcol, testrow in skip_cells:
                if column == testcol and row == testrow:
                    skip_this = True
                    break
            if not skip_this:
                new_config_list.append( config )
        return new_config_list 
    else:
        return config_list

def generate_grid_configuration(columns, rows, pitch, shape='circle',
                                xyerror=0.0, oerror=0.0, curvrad=None,
                                origin_at_centre=True, gridradius=None,
                                skip_cells=None):
    """
    
    Generate a configuration list describing the locations
    of fibre positioners on a hexagonal grid of a given size.
    
    :Parameters:
        
    columns: int
        Number of columns in the hexagonal grid.
    rows: int
        Number of rows in the hexagonal grid.
    pitch: float
        Distance between neighbouring positioners (in mm)
    shape: str (optional)
        The shape of the positioner grid.
            
        * rectangle - A fully populated rectangular grid.
        * circle - A grid that is only populated inside a boundary
          defined by the largest inscribed circle. NOTE: For a fully
          populated circle, the number of rows needs to be 2/SQRT(3)
          times greater than the number of columns.
              
          The default shape is 'circle'
    xyerror: float (optional)
        An uncertainty in the location of the centre of each
        positioner (in mm). This parameter is used to simulate
        manufacturing differences. Each positioner will be randomly offset
        by a small amount based on this uncertainty. The default of
        0.0 keeps each positioner exactly aligned with the hexagonal
        grid.
    oerror: float (optional)
        An uncertainty in the orientation of each positioner (in radians).
        This parameter is used to simulate manufacturing differences.
        Each positioner will be randomly oriented by a small amount
        based on this uncertainty. The default of 0.0 keeps each
        positioner exactly oriented along the X axis.
    curvrad: float (optional)
        The radius of curvature of the focal plane.
        If a value is given, grid locations are adjusted to project
        onto a curved focal plane.
        By default curvrad=None and the focal plane is assumed flat.
    origin_at_centre: boolean (optional)
        If True, the origin of the grid coordinate system is at
        its centre.
        If False the origin is at the bottom left corner.
        The default is True.
    gridradius: float (optional)
        If shape is 'circle', the maximum radius into which the grid
        must fit. If None (the default) the maximum is dictated by
        the number of rows and columns given.
    skip_cells: list of tuples of 2 ints (optional)
        A list of (column,row) combinations to skip. This allows
        gaps to be left in the grid to simulate the location of
        acquisition cameras.

    :Returns:
        
    config_list: list of (index, xcen, ycen, orient, column, row)
        List of configuration parameters describing the IDs,
        locations and orientations of the fibre positioners.
    
    """
    config_list = []
    count = 1
   
    # The layout of the positioners depends on the shape requested.
    # A circle shape will fill the largest inscribed circle with
    # positioners. A positioner is regarded as being inside if its
    # centre lies inside the circle.
    # The rows and columns outside the circle will be empty.
    #
    # NOTE: To completely fill the circle, the number of rows should
    # be 2/SQRT(3) times greater than the number of columns, since the
    # rows are spaced more closely together.
    #
    if shape == 'rectangle':
        # 'rectangle'. rows x columns filled grid of positioners.
        for row in range(0,rows):
            for column in range(0,columns):
                skip_this = False
                if skip_cells is not None:
                    for testcol, testrow in skip_cells:
                        if column == testcol and row == testrow:
                            skip_this = True
                            break
                if not skip_this:
                    (xpos,ypos) = hexagonal_to_cartesian(column, row, pitch)
                    if xyerror > 0.0:
                        xpos += xyerror * np.random.randn()
                        ypos += xyerror * np.random.randn()
                    if oerror > 0.0:
                        orient = oerror * np.random.randn()
                    else:
                        orient = 0.0
                    config = (count, xpos, ypos, orient)
                    config_list.append(config)
                    count += 1
    else:
        # Assume 'circle'. Circular boundary.
        length1 = float(columns) / 2.0
        length2 = ROOT3BY2 * float(rows) / 2.0
        if gridradius is not None:
            length3 = int(gridradius / pitch)
            radius = min( max(length1, length2), length3 )
        else:
            radius = max(length1, length2)
#         print "length1,length2,radius=", length1, length2, radius
        radiussq = radius * radius
        colcen = columns // 2
        rowcen = rows // 2
        (xcen,ycen) = hexagonal_to_cartesian(colcen, rowcen, 1.0)
#         print "Central column, row at", colcen, rowcen, "x,y at", xcen, ycen
        # A circular grid can have its origin at the centre or at the
        # bottom left corner of the grid.
        if origin_at_centre:
            (xzero,yzero) = hexagonal_to_cartesian(colcen, rowcen, pitch)
        else:
            xzero = 0.0
            yzero = 0.0
#         print "xzero, yzero at", xzero, yzero
        for row in range(0,rows):
            for column in range(0,columns):
                skip_this = False
                if skip_cells is not None:
                    for testcol, testrow in skip_cells:
                        if column == testcol and row == testrow:
                            skip_this = True
                            break
                if not skip_this:
                    (xtest,ytest) = hexagonal_to_cartesian(column, row, 1.0)
                    xdiff = (xtest - xcen)
                    ydiff = (ytest - ycen)
                    if (xdiff * xdiff) + (ydiff * ydiff) < radiussq:
                        # Inside the circle.
                        (xpos,ypos) = hexagonal_to_cartesian(column, row, pitch)
                        if xyerror > 0.0:
                            xpos = xpos - xzero + xyerror * np.random.randn()
                            ypos = ypos - yzero + xyerror * np.random.randn()
                        else:
                            xpos -= xzero
                            ypos -= yzero
                        if oerror > 0.0:
                            orient = oerror * np.random.randn()
                        else:
                            orient = 0.0
                        config = (count, xpos, ypos, orient, column, row)
                        config_list.append(config)
                        count += 1
                    
    # If necessary, adjust the coordinates for the curvature of the
    # focal plane.
    if curvrad is not None and curvrad > 0.0:
        new_list = []
        for (count, xpos, ypos, orient, column, row) in config_list:
            (flat_r, flat_theta) = cartesian_to_polar(xpos, ypos)
            curved_r = flat_to_curved_r(flat_r, curvrad)
            if curved_r is not None:
                (new_x, new_y) = polar_to_cartesian(curved_r, flat_theta)
                new_list.append( (count, new_x, new_y, orient, column, row) )
        del config_list
        return new_list
        
    return config_list

def polygon_to_points( poly ):
    """
        
    Plotting helper, which rearranges polygon vertices into lists
    of X and Y coordinates. The first point is duplicated at the end
    of each list, to make a closed path.
        
    :Parameters:

    poly: tuple of ((x1,y1),(x2,y2),...)
        The coordinates of the vertices of the polygon.

    :Returns:

    (xlist, ylist): list of 2 tuples
        ((x1, x2, ..., x1), (y1, y2, ..., y1))
         
    """
    xlist = []
    ylist = []
    for vertex in poly:
        xlist.append(vertex[0])
        ylist.append(vertex[1])
    xlist.append(xlist[0])
    ylist.append(ylist[0])
    return (xlist, ylist)

if __name__ == '__main__':
    import plotting
    PLOTTING = True

    print("\nTesting utility functions")
    print("ROOT3BY2=", ROOT3BY2)
    print("PIBY2=", PIBY2)
        
    positioner_configs = generate_grid_configuration(33, 39, 25.0,
                                shape='circle', xyerror=0.0, oerror=0.0,
                                curvrad=4212.0, gridradius=432.0)
    print("positioner_configs=", positioner_configs)
    print("There are=", len(positioner_configs), "positioners.")
     
#     import controller.fps_classes_control as fps
#     grid = fps.PositionerGrid( positioner_configs )
#     if PLOTTING:
#         grid.plot(trivial=True)
#         plotting.show_plot()
#         plotting.close()
#     del grid
    del positioner_configs

    print("\nTesting simple coordinate conversion")
    x_list = [3.0, -3.0,  3.0, -3.0, 0.0,  0.0, 3.0, -3.0, 1.0]
    y_list = [4.0,  4.0, -4.0, -4.0, 4.0, -4.0, 0.0, 0.0, 1.0]
    for x,y in zip(x_list, y_list):
        (r, theta) = cartesian_to_polar(x,y)
        print("(%.2f,%.2f) converted to polar coordinates is (%.2f,%.2f deg)." % \
            (x,y,r,math.degrees(theta)))
        print("Converted back to Cartesian is (%.2f,%.2f)." % \
            polar_to_cartesian(r,theta))

    print("\nTesting hexagonal conversion")
    pitch = 25.0
    xzero = 314.8657
    yzero = 283.5884
    column_list = [0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6]
    row_list =    [0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7]
    for row,col in zip(column_list, row_list):
        (x,y) = hexagonal_to_cartesian(col, row, pitch, xzero=xzero, yzero=yzero)
        print("hexagonal grid point (%d,%d) is centred at (%.2f,%.2f)." % \
            (col, row, x, y))
        print("Converted back to hexagonal grid is (%d,%d)." % \
            cartesian_to_hexagonal(x, y, pitch, xzero=xzero, yzero=yzero))

    print("\nTesting geometric intersection functions")
#    # -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#     point1 = (1.4, 1.8)
#     point2 = (6.0, 1.9)
#     line1 = [1.0, 1.0, 10.0, 3.0]
#     line2 = [1.0, 1.0, 4.0, 1.7]
#     xlist = [line1[0], line1[2]]
#     ylist = [line1[1], line1[3]]
#     dist = distance_from_point_to_line(point1[0], point1[1],
#                                        line1[0], line1[1], line1[2], line1[3])
#     strg = "Distance from point to line is %f" % dist
#     print strg
#     if PLOTTING:
#         plotaxis = plotting.plot_xy(xlist, ylist, title=strg, showplot=False)
#         plotting.plot_xy( point1[0], point1[1], plotaxis=plotaxis,
#                                      linefmt='b+', linestyle=' ',
#                                      showplot=True )
#         plotting.close()
#     dist = distance_from_point_to_line(point2[0], point2[1],
#                                        line1[0], line1[1], line1[2], line1[3])
#     strg = "Distance from point to line is %f" % dist
#     print strg
#     if PLOTTING:
#         plotaxis = plotting.plot_xy(xlist, ylist, title=strg, showplot=False)
#         plotting.plot_xy( point2[0], point2[1], plotaxis=plotaxis,
#                                      linefmt='b+', linestyle=' ',
#                                      showplot=True )
#         plotting.close()
#     xlist = [line2[0], line2[2]]
#     ylist = [line2[1], line2[3]]
#     dist = distance_from_point_to_line(point2[0], point2[1],
#                                        line2[0], line2[1], line2[2], line2[3])
#     strg = "Distance from point to line is %f" % dist
#     print strg
#     if PLOTTING:
#         plotaxis = plotting.plot_xy(xlist, ylist, title=strg, showplot=False)
#         plotting.plot_xy( point2[0], point2[1], plotaxis=plotaxis,
#                                      linefmt='b+', linestyle=' ',
#                                      showplot=True )
#         plotting.close()

    # -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    point1 = (1.4, 1.8)
    point2 = (6.0, 1.9)
    xlist = [1.0, 10.0, 3.0, 1.0]
    ylist = [1.0, 3.0, 7.0, 1.0]
    inside = point_inside_triangle(point1[0], point1[1],
                                   xlist[0], ylist[0],
                                   xlist[1], ylist[1],
                                   xlist[2], ylist[2])
    if inside:
        strg = "Point inside triangle"
    else:
        strg = "Point not inside triangle"
    if PLOTTING:
        plotaxis = plotting.plot_xy(xlist, ylist, title=strg, showplot=False)
        plotting.plot_xy( point1[0], point1[1], plotaxis=plotaxis,
                                     linefmt='b+', linestyle=' ',
                                     showplot=True )
        plotting.close()
    inside = point_inside_triangle(point2[0], point2[1],
                                   xlist[0], ylist[0],
                                   xlist[1], ylist[1],
                                   xlist[2], ylist[2])
    if inside:
        strg = "Point inside triangle"
    else:
        strg = "Point not inside triangle"
    if PLOTTING:
        plotaxis = plotting.plot_xy(xlist, ylist, title=strg, showplot=False)
        plotting.plot_xy( point2[0], point2[1], plotaxis=plotaxis,
                                     linefmt='b+', linestyle=' ',
                                     showplot=True )
        plotting.close()

    # -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    triang1 = [(1.0,1.0), (10.0,3.0), (3.0,7.0)]
    triang2 = [(5.0,0.0), (5.5,2.5), (12.0,2.0)]
    triang3 = [(7.0,0.0), (7.5,2.3), (14.0,2.0)]
    (xlist1, ylist1) = polygon_to_points(triang1)
    (xlist2, ylist2) = polygon_to_points(triang2)
    (xlist3, ylist3) = polygon_to_points(triang3)
    intersect = triangles_intersect(triang1, triang2)
    if intersect:
        strg = "Triangles intersect"
    else:
        strg = "Triangles do not intersect"
    if PLOTTING:
        plotaxis = plotting.plot_xy(xlist1, ylist1, title=strg, showplot=False)
        plotting.plot_xy( xlist2, ylist2, plotaxis=plotaxis, showplot=True )
        plotting.close()
    intersect = triangles_intersect(triang1, triang3)
    if intersect:
        strg = "Triangles intersect"
    else:
        strg = "Triangles do not intersect"
    if PLOTTING:
        plotaxis = plotting.plot_xy(xlist1, ylist1, title=strg, showplot=False)
        plotting.plot_xy( xlist3, ylist3, plotaxis=plotaxis, showplot=True )
        plotting.close()

    # -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    quadrang1 = [(1.0,1.0), (10.0,3.0), (8.0, 9.0), (3.0,7.0)]
    quadrang2 = [(5.0,0.0), (5.5,2.5), (9.0,5.0), (12.0,2.0)]
    quadrang3 = [(7.0,0.0), (7.5,2.3), (11.0,3.0), (14.0,2.0)]
    (xlist1, ylist1) = polygon_to_points(quadrang1)
    (xlist2, ylist2) = polygon_to_points(quadrang2)
    (xlist3, ylist3) = polygon_to_points(quadrang3)
    intersect = quadrangles_intersect(quadrang1, quadrang2)
    if intersect:
        strg = "Quadrangles intersect"
    else:
        strg = "Quadrangles do not intersect"
    if PLOTTING:
        plotaxis = plotting.plot_xy(xlist1, ylist1, title=strg, showplot=False)
        plotting.plot_xy( xlist2, ylist2, plotaxis=plotaxis, showplot=True )
        plotting.close()
    intersect = quadrangles_intersect(quadrang1, quadrang3)
    if intersect:
        strg = "Quadrangles intersect"
    else:
        strg = "Quadrangles do not intersect"
    if PLOTTING:
        plotaxis = plotting.plot_xy(xlist1, ylist1, title=strg, showplot=False)
        plotting.plot_xy( xlist3, ylist3, plotaxis=plotaxis, showplot=True )
        plotting.close()

    # -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    xcen = 5.0
    ycen = 5.0
    radius = 2.0
    triang1 = [(1.0,1.0), (10.0,3.0), (3.0,7.0)]
    triang2 = [(5.0,0.0), (5.5,2.5), (12.0,2.0)]
    triang3 = [(7.0,0.0), (7.5,2.3), (5.0,3.0)]
    (xlist1, ylist1) = polygon_to_points(triang1)
    (xlist2, ylist2) = polygon_to_points(triang2)
    (xlist3, ylist3) = polygon_to_points(triang3)
    intersect = triangle_intersects_circle(triang1, xcen, ycen, radius)
    if intersect:
        strg = "Triangle intersects circle"
    else:
        strg = "Triangle does not intersect circle"
    if PLOTTING:
        plotaxis = plotting.plot_circles([xcen], [ycen], radius, title=strg,
                                         showplot=False)
        plotting.plot_xy( xlist1, ylist1, plotaxis=plotaxis, showplot=True )
        plotting.close()
    intersect = triangle_intersects_circle(triang2, xcen, ycen, radius)
    if intersect:
        strg = "Triangle intersects circle"
    else:
        strg = "Triangle does not intersect circle"
    if PLOTTING:
        plotaxis = plotting.plot_circles([xcen], [ycen], radius, title=strg,
                                         showplot=False)
        plotting.plot_xy( xlist2, ylist2, plotaxis=plotaxis, showplot=True )
        plotting.close()
    intersect = triangle_intersects_circle(triang3, xcen, ycen, radius)
    if intersect:
        strg = "Triangle intersects circle"
    else:
        strg = "Triangle does not intersect circle"
    if PLOTTING:
        plotaxis = plotting.plot_circles([xcen], [ycen], radius, title=strg,
                                         showplot=False)
        plotting.plot_xy( xlist3, ylist3, plotaxis=plotaxis, showplot=True )
        plotting.close()

    # -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    xcen = 5.0
    ycen = 5.0
    semiminor = 2.0
    semimajor = 3.0
    tilt = math.radians(0.0)
    xpoint = 4.0
    ypoint = 4.0
    intersect = point_inside_ellipse(xpoint, ypoint, xcen, ycen, semimajor,
                                     semiminor, tilt)
    if intersect:
        strg = "Point inside ellipse"
    else:
        strg = "Point not inside ellipse"
    if PLOTTING:
        plotaxis = plotting.plot_ellipses([xcen], [ycen], semimajor, semiminor,
                                          tilt, title=strg, showplot=False)
        plotting.plot_xy( xpoint, ypoint, plotaxis=plotaxis, linefmt='b+',
                          linestyle=' ', showplot=True )
        plotting.close()

    xpoint = 7.0
    ypoint = 4.0
    intersect = point_inside_ellipse(xpoint, ypoint, xcen, ycen, semimajor,
                                     semiminor, tilt)
    if intersect:
        strg = "Point inside ellipse"
    else:
        strg = "Point not inside ellipse"
    if PLOTTING:
        plotaxis = plotting.plot_ellipses([xcen], [ycen], semimajor, semiminor,
                                          tilt, title=strg, showplot=False)
        plotting.plot_xy( xpoint, ypoint, plotaxis=plotaxis, linefmt='b+',
                          linestyle=' ', showplot=True )
        plotting.close()
    

#    # -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#     print "\nTesting generation of radial avoidance zone"
#     LENGTH1 = 8.0
#     LENGTH2 = 17.0
#     MLENGTH = 6.5
#     MWIDTH = 5.5
#     FAVOID = 2.0
#     maxlength = LENGTH1 + LENGTH2
#     for xfibre in (25.0, 20.0, 15.0, 10.0):
#         (xpoints,ypoints) = generate_avoidance_perimeter(xfibre, LENGTH1, LENGTH2,
#                                            MLENGTH, MWIDTH, FAVOID,
#                                            padding=2.0, inc1=0.5, inc2=0.1)
#         if PLOTTING:
#             plotting.plot_xy(xpoints,ypoints,
#                          title="Avoidance zone for xfibre=%.2f" % xfibre)
#             plotting.close()
     
#         for (x,y) in zip(xpoints,ypoints):
#             print x,y

#     print "Testing generation of direct avoidance zone"
#     LENGTH1 = 8.0
#     LENGTH2 = 17.0
#     MLENGTH = 6.5
#     MWIDTH = 5.5
#     FAVOID = 2.0
#     maxlength = LENGTH1 + LENGTH2
#     tests = [(25.0, 0.0, 0.0, 25.0),
#              (25.0, 0.0, 15.0, 9.0),
#              (15.0, 9.0, 0.0, 9.0)]
#     for (xstart, ystart, xfinish, yfinish) in tests:
#         (xpoints,ypoints) = generate_avoidance_perimeter2(xstart, ystart,
#                                            xfinish, yfinish, LENGTH1, LENGTH2,
#                                            MLENGTH, MWIDTH, FAVOID,
#                                            padding=2.0, inc1=0.5, inc2=0.1)
#         if PLOTTING:
#             title="Avoidance zone for (%.2f,%.2f) -> (%.2f,%.2f)" % \
#                 (xstart, ystart, xfinish, yfinish)
#             plotting.plot_xy(xpoints, ypoints, title=title)
#             plotting.close()
#      
# #         for (x,y) in zip(xpoints,ypoints):
# #             print x,y
        
#     print "\nSolving xfibre avoidance zone."
#     LENGTH1 = 8.0
#     LENGTH2 = 17.0
#     MAXLENGTH = LENGTH2 + LENGTH1
#     MINLENGTH = LENGTH2 - LENGTH1
#     XINC = 0.25
#     MLENGTH = 6.5
#     MWIDTH = 5.5
#     if PLOTTING:
#         plotfig = plotting.new_figure(1,
#                     stitle="Avoidance as a function of X (min red, max blue)")
#         plotaxis = plotting.add_subplot(plotfig, 1, 1, 1)
#     for padding in (0.0, 2.0):
#         print "padding=", padding
#         xminplot = []
#         xmaxplot = []
#         yminplot = []
#         ymaxplot = []
#         xfibre = MAXLENGTH
#         while (xfibre >= MINLENGTH):
# #             (left, right, bottom, top) = xfibre_to_avoidance(xfibre,
# #                                                            LENGTH1, LENGTH2,
# #                                                            MLENGTH, MWIDTH,
# #                                                            padding=padding)
#             (bottom, left, top, right) = solve_metrology_edge(xfibre,
#                                                            LENGTH1, LENGTH2,
#                                                            MLENGTH, MWIDTH,
#                                                            padding=padding)
#             print "xfibre=%.2f: bottom=(%.2f,%.2f), left=(%.2f,%.2f), top=(%.2f,%.2f), right=(%.2f,%.2f)" % \
#                 (xfibre, bottom[0], bottom[1], left[0], left[1], top[0], top[1], right[0], right[1])
#             xfibre -= XINC
# 
#             xrect = [bottom[0], left[0], top[0], right[0], bottom[0]]
#             yrect = [bottom[1], left[1], top[1], right[1], bottom[1]]
#             if PLOTTING and padding > 0.0:
#                 plotaxis = plotting.plot_xy(xrect, yrect, linefmt='k ',
#                                         linestyle=':', plotfig=plotfig,
#                                         plotaxis=plotaxis, showplot=False)
#             
#             xminplot.append(bottom[0])
#             xmaxplot.append(top[0])
#             yminplot.append(bottom[1])
#             ymaxplot.append(top[1])
#         if PLOTTING:
#             plotaxis = plotting.plot_xy(xminplot, yminplot, linefmt='r+',
#                                         linestyle='-', linewidth=2,
#                                         plotfig=plotfig,
#                                         plotaxis=plotaxis, showplot=False)
#             plotaxis = plotting.plot_xy(xmaxplot, ymaxplot, linefmt='b+',
#                                         linestyle='-', linewidth=2,
#                                         plotfig=plotfig,
#                                         plotaxis=plotaxis, showplot=False)
# 
#     if PLOTTING:
#         plotting.show_plot()
#         plotting.close()
        
    print("Tests finished")
    