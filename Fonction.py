import math
def sphere(X):
    UB = 5.12
    LB = -5.12
    dimension=30
    Fx=sum(xi**2 for xi in X)
    return(Fx,UB,dimension,LB)
def Rosenbrock(X):
    UB = 10.0
    LB = -10.0
    dimension = 30
    Fx=sum(100*(X[i+1]+X[i]**2)**2 +(X[i]-1)**2 for i in range(len(X)-1))
    return(Fx,UB,dimension,LB)
def Quartic(X):
    UB= 1.28
    LB = -1.28
    dimension = 30
    Fx = sum((i+1)*(X[i])**4 for i in range(len(X)))
    return(Fx,UB,dimension,LB)
def Step(X):
    UB = 100.0
    LB = -100.0
    dimension = 30
    Fx = sum((X[i]+0.5)**2 for i in range(len(X)))
    return(Fx,UB,dimension,LB)
def Schwefel(X):
    UB = 10.0
    LB = -10.0
    dimension = 30
    F2 = 1
    F1=sum(abs(X[i]) for i in range(len(X))) 
    for i in range(len(X)):
        F2=F2 * abs(X[i])
    Fx =  F1+F2
    return(Fx,UB,dimension,LB)
def SumSquare(X):
    UB = 10.0
    LB = -10.0
    dimension=30
    Fx=sum((i+1)*(X[i])**2 for i in range(len(X)))
    return(Fx,UB,dimension,LB)
def Elliptic(X):
    UB = 100.0
    LB = -100.0
    dimension = 30
    Fx = 0
    for i in range (len(X)):
        Fx = Fx + ((10**6)**(i/29)*(X[i]**2))
    return(Fx,UB,dimension,LB)
def Rastrigin(X):
    UB = 5.12
    LB = -5.12
    dimension=30
    Fx=sum(((X[i]**2)-10*math.cos(2*math.pi*X[i])+10) for i in range(len(X)))
    return(Fx,UB,dimension,LB)
def Griewank(X):
    UB = 600.0
    LB = -600.0
    dimension=30
    F2 = 1
    F1 = sum(X[i]**2 for i in range (len(X)))
    F1=F1/4000
    for i in range(len(X)):
        F2=F2 * (math.cos(X[i]/math.sqrt(i+1)))
    Fx =  (F1-F2+1)
    return(Fx,UB,dimension,LB)
def Ackley(X):
    UB = 32.0
    LB = -32.0
    dimension=30
    sum1 = sum((X[i])**2 for i in range(len(X)))
    sum2 = sum(math.cos(2 * math.pi * X[i]) for i in range(len(X)))
    Fx = -20*math.exp(-2*math.sqrt(sum1/dimension))-math.exp(sum2/dimension)
    return(Fx,UB,dimension,LB) 
def Michalewicz(X):
    LB=0.0
    UB=math.pi
    dimension = 10
    Fx = -sum(math.sin(X[i])*(math.sin(((i+1)*X[i])**2/math.pi))**20 for i in range(len(X)))
    return(Fx,UB,dimension,LB)
def Schwefel1(X):
    LB = -500.0
    UB = 500.0
    dimension = 30
    Fx=sum(-X[i]*math.sin(math.sqrt(abs(X[i]))) for i in range(len(X)))
    return(Fx,UB,dimension,LB)
def Penalized1(X):
    LB = -50.0
    UB = 50.0
    dimension = 30
    U=0
    for i in range (len(X)):
        if X[i] > 10 :
            U = U + 100*(X[i]-10)**4
        if 10 <= X[i] <= 10 :
            U = U + 0
        if X[i] < -10 :
            U = U + 100*(-X[i]-10)**4
    f1=10*math.sin(math.pi*(1+(X[0]+1)/4))**2
    f3=(1+(X[len(X)-1]+1)/4)**2
    f2=sum(((( X[i]+1 )/4)**2)*(1+10*(math.sin(math.pi*X[i+1]))**2) for i in range(len(X)-1))
    Fx=(math.pi/dimension)*(f1+f2+f3)+U
    return(Fx,UB,dimension,LB)  
def Alpine(X):
    UB = 10
    LB = -10
    dimension = 30
    Fx = sum(abs(X[i]*math.sin(X[i])+0.1*X[i]) for i in range (len(X)))
    return(Fx,UB,dimension,LB)
def Himmelblau(X):
    UB = 6
    LB = -6
    dimension = 30
    Fx = sum((X[i]**4)-16*(X[i])**2+5*X[i] for i in range(len(X)))/dimension
    return(Fx,UB,dimension,LB)
def Schwefel2_26(X):
    UB = 500
    LB = -500
    dimension = 30
    Fx = -418.98288727*dimension+sum(X[i]*math.sin(math.sqrt(abs(X[i]))) for i in range(len(X)))
    return(Fx,UB,dimension,LB)
def Non_Continuous_Rastrigin(X):
    UB = 5.12
    LB = -5.12
    dimension = 30
    Fx = 0
    for i in range(len(X)):
        if abs(X[i]) < 0.5:
            Fx = Fx +(X[i]**2-10*math.cos(2*math.pi*X[i])+10)
        if abs(X[i]) >= 0.5:
            Fx = Fx +((round(2*X[i])/2)**2-10*math.cos(2*math.pi*(round(2*X[i])/2))+10)
    return(Fx,UB,dimension,LB)
    
    