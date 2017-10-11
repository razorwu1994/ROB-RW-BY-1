

import math
import itertools
import numpy as np

#helper function
def GCD(a, b):
    if b == 0:
        return a
    else:
        return GCD(b, a % b)

class plane(object):
    #plane equation : Ax+By+Cz+D=0
    def __init__(self, A, B, C,D):
        self.A = A
        self.B = B
        self.C = C
        self.D = D

class point(object):
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z


class sphere(object):
    def __init__(self, point, radius):
        self.center = point
        self.radius = radius

def get_two_points_distance(p1, p2):
    return math.sqrt(pow((p1.x - p2.x), 2) + pow((p1.y - p2.y), 2) + pow((p1.z - p2.z), 2))


def roundFunction(number):
    if abs(number)==number:number=abs(number)

    if abs(number-round(number,5))<math.pow(10,-8):
        return round(number,5)
    else:
        return number

def get_two_spheres_intersecting_plane(c1, c2):
    p1 = c1.center
    p2 = c2.center
    r1 = c1.radius
    r2 = c2.radius


    d = get_two_points_distance(p1, p2)
    # print("distance between ("+str(p1.x)+","+str(p1.y)+","+str(p1.z)+") and ("
    #       +str(p2.x)+","+str(p2.y)+","+str(p2.z) +") is "+str(d))

    # if to far away, or self contained - can't be done
    if d >= (r1 + r2) or d <= math.fabs(r1 - r2):
        return None

    # a = math.acos((pow(r1, 2)+ pow(d, 2) - pow(r2, 2) ) / (2 * d * r1))
    # h = r1*math.sin(a)

    A= 2*(p2.x - p1.x)
    B= 2*(p2.y - p1.y)
    C= 2*(p2.z - p1.z)
    D= math.pow(p1.x,2)-math.pow(p2.x,2)+math.pow(p1.y,2)-math.pow(p2.y,2)+math.pow(p1.z,2)-math.pow(p2.z,2)-math.pow(r1,2)+math.pow(r2,2)

    tempA= A
    tempB= B
    tempC= C
    tempD= D

    gcd = reduce(GCD, (tempA, tempB, tempC, tempD))
    if abs(gcd) > 1:
        tempA=  tempA / gcd
        tempB = tempB / gcd
        tempC = tempC / gcd
        tempD = tempD / gcd

    # print "gcd is " + str(gcd)
    TempPlane = plane(tempA,tempB,tempC,tempD)

    print("intersecting plane equation : "+str(tempA)+"X+"+str(tempB)+"Y+"+str(tempC)+"Z+"+str(tempD)+"=0")

    return [TempPlane]


def get_all_intersecting_planes(spheres):
    planes = []
    num = len(spheres)
    for i in range(num):
        j = i + 1
        for k in range(j, num):
            res = get_two_spheres_intersecting_plane(spheres[i], spheres[k])
            if res:
                planes.extend(res)
    return planes


def is_contained_in_spheres(point, spheres):
    # print str(point.x)+","+str(point.y)+","+str(point.z)
    for i in range(len(spheres)):
        if (get_two_points_distance(point, spheres[i].center) > (spheres[i].radius)):
            print "distance"+str(get_two_points_distance(point, spheres[i].center))
            print str(spheres[i].radius)
            return False
    return True



def check_duplicate(points,testPoint):
    if len(points)==0:return False

    completeDupFlag=True
    for p in points:
        if p.x != testPoint.x:
            completeDupFlag=False
            break
        if p.y != testPoint.y:
            completeDupFlag = False
            break
        if p.z != testPoint.z:
            completeDupFlag = False
            break
    return completeDupFlag

if __name__ == '__main__':

    p1 = point(1.0, 0,0)
    p2 = point(0, 1.0,0)
    p3 = point(0, 0,1.0)
    p4 = point(-1.0, 0,0)

    c1 = sphere(p1, 1.0)
    c2 = sphere(p2, 1.0)
    c3 = sphere(p3, 1.0)
    c4 = sphere(p4, 1.0)

    sphere_list = [c1, c2,c3]
    #, c3,c4

    #print(get_all_intersecting_points(circle_list)[0].z)

    # inner_points = []
    counter=0
    planes = get_all_intersecting_planes(sphere_list)

    points=[]
    for singleCombinations in itertools.combinations(planes, 3):
        # planeCounter=0
        planesMatrix=[]
        resultXdetMatrix= []
        resultYdetMatrix = []
        resultZdetMatrix = []
        resultX_point=0.0
        resultY_point=0.0
        resultZ_point=0.0
        for plane in singleCombinations:
            planesMatrix.append([plane.A,plane.B,plane.C])
            resultXdetMatrix.append([plane.D,plane.B,plane.C])
            resultYdetMatrix.append([plane.A,plane.D,plane.C])
            resultZdetMatrix.append([plane.A,plane.B,plane.D])

            # print "plane "+str(planeCounter)+" equation : "+str(plane.A)+"X+"+str(plane.B)+"Y+"+str(plane.C)+"Z+"+str(plane.D)+"=0"
            # planeCounter+=1

        a = np.array(planesMatrix)
        result_det = float(np.linalg.det(a))
        #if det is 0 , that means these 3 planes has no single point intersection
        if result_det == 0: pass
        #we got the point
        else:
            # print "found the point"
            #according to cramer's rule, we found the point and remove difference between -0.0 and 0.0
            resultX_point=roundFunction(float(np.linalg.det(np.array(resultXdetMatrix)))/result_det)
            resultY_point=roundFunction(float(np.linalg.det(np.array(resultYdetMatrix)))/result_det)
            resultZ_point=roundFunction(float(np.linalg.det(np.array(resultZdetMatrix)))/result_det)
            tempPoint = point(resultX_point,resultY_point,resultZ_point)
            if not check_duplicate(points,tempPoint):
                points.append(tempPoint)
                print "at "+"("+str(resultX_point)+","+str(resultY_point)+","+str(resultZ_point)+")"
            else:
                print "duplicate point"
        print "\n"

    if len(points)==0: print "no point found reason could be : intersecting at a plane or a line, or no intersection at all"
    for p in points:
        counter+=1
        print str(counter)+" points"


