import numpy as np 

#The purpose of this program is to find the zero point of any elementary functions (I mean any)
#Which include log,exp,cos,sin and polynomial
#This program is powered by the library which is called sympy
#author-----------------Zilong-------------------------------
#time-------------------2019/3/23---------------------------
#program name-----------find_zero_point.py-------------------
#First of all, I import some necessary library of python.
import numpy as np 
import matplotlib.pyplot as plt 
import sympy 
from sympy import log,exp,cos,sin


#To handle more questions, I define a class
#This class includes newton's method and dichotomy
class Find_Zero_Point(object):
    def __init__(self,expr,accuracy,begin_point):
        #the expr means the expression of a functions
        self.expr=expr
        #By using the sympy library, we calculate the differential of the function 
        self.expr_diff=expr.diff(x)
        #We need a begin point when using newton's method
        self.point=begin_point
        #Get the accuracy we need.
        self.accuracy=accuracy
        #Calculate the number of the cycles(newton's method)
        self.count_newton=0
        #dichotomy method 
        self.count_dichotomy=0

    
    def Newton_Method(self):
        #calculate the error, if error is less than accuracy,then output
        error=np.abs(float(expr.subs(x,self.point)))
        if error<=self.accuracy:
            print("It takes",self.count_newton,"times to get the result satisfy the accuracy.")
            print("The result of f(x)=0 is",self.point,".")
        else:
            #the main body of newton's method
            self.count_newton+=1
            self.point=self.point-float(expr.subs(x,self.point))*1/(float(self.expr_diff.subs(x,self.point)))
            self.Newton_Method()

    def Dichotomy(self,point_1,point_2):
        #Unfortunately, the dichotomy method need two points, which means, f(point_1)*f(point_2)<0 and point_1<point_2
        point_3=(point_1+point_2)/2
        if (point_2-point_1)/2<=self.accuracy:
            print("The result is",point_3,".")
            print("It takes",self.count_dichotomy,"times to get the result satisfing the accuracy.")
            return
        self.count_dichotomy+=1
        result_mid=expr.subs(x,point_3)
        if result_mid==0:
            print("This is the accurate result",point_3,".")
            print("It takes",self.count_dichotomy,"times to get the result satisfing the accuracy.")
        elif result_mid*expr.subs(x,point_2)<0:
            self.Dichotomy(point_3,point_2)
        else:
            self.Dichotomy(point_1,point_3)


#To use the sympy library, specifying parameter is necessary.
x=sympy.Symbol("x",real=True)            
expr=log(x)*exp(x)-x*x
#A object
zero_point=Find_Zero_Point(expr,0.01,2)
#Calculate by using newton'method and dichotomy method 
zero_point.Newton_Method()
zero_point.Dichotomy(1,2)